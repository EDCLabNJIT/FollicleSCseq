library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(EnhancedVolcano)
library(readxl)
library(ggpubr)
library(grid)

#' Creates Volcanos
#' 
#' @return nothing
#' 
#' @param data input data (list of `data.frames`s)
#' @param DEFolder folder for DE data if `df` is `NULL`
#' @param plotLocation folder to plot images
#' @param level level(s) to search at
#' @param pCutoff cutoff adjusted p for volcano
#' @param FCutoff cutoff log2fold change for volcano
makevolcano <- function(data = NULL, DEFolder = "./data/",
                        ylab = "FDR", ycol = "FDR",
                        plotLocation = "./images/",
                        level = 0, pCutoff = 0.05, FCutoff = 1,
                        xpadding = 1.5, ypadding = 5, col = c("grey30", "forestgreen", "royalblue", "red2")) {

if (is.null(data)){
filename <- paste0(DEFolder,"DE_level",level,".xlsx")
sheetnames <- excel_sheets(filename)
} else {
  sheetnames = names(data)
}
plts <- c()

for (sheetname in sheetnames) {
  df <- read_xlsx(filename, sheetname, n_max = 1000)
  
  plts[[length(plts)+1]] <- EnhancedVolcano(df ,x ="avg_log2FC", y=ycol, lab=df$symbol, pCutoff = pCutoff, drawConnectors = TRUE,
                                            title = sheetname, subtitle = NULL,xlab = NULL,ylab = NULL, FCcutoff = FCutoff,
                                            caption = NULL, colAlpha = 1, col = col, pointSize = 4, labSize = 5,
                                            xlim = c(min(df[["avg_log2FC"]], na.rm = TRUE) - xpadding, max(df[["avg_log2FC"]], na.rm = TRUE) + xpadding),
                                            ylim = c(0, max(-log10(df[[ycol]]), na.rm = TRUE) + ypadding)) + theme(legend.position = "none")
}

figure <- ggarrange(plotlist = plts)

ann_figure <- annotate_figure(figure, left = textGrob(bquote(~-Log[10] ~ italic(.(ylab))), rot = 90, vjust = .75, gp = gpar(cex = 3)),
                bottom = textGrob(bquote(~Log[2] ~ "fold change"), gp = gpar(cex = 3))) + theme_bw() + plot_layout(guides = "collect")

out_file = paste0(plotLocation,"volcano_",level,"_",pCutoff,"_",FCutoff,".png")
if (level == 0) {
  ggsave(out_file, plot = ann_figure, width = 12, height = 8, dpi = 300, units = "in")
} else {
  ggsave(out_file, plot = ann_figure, width = 16, height = 16, dpi = 300, units = "in")
}

}