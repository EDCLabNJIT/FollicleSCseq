library(readxl)
library(ggplot2)

#' Creates dot plot from DAVID results given a specified pathway file
#' 
#' @return nothing
#' 
#' @param davidresults results text file from DAVID
#' @param keptids file (or vector) of ids to plot
#' @param plotname name of plot
#' @param height height of plot (in inches)
davidplot <- function(davidresults = "DAVIDresults/gran chart defaults.txt", keptids = "DAVIDresults/gran ids.txt",
                        plotname = "./images/david_gran.png", height = 4.5) {

t <- read.delim(davidresults)
colnames(t)[4] <- "Percent"
if (typeof(keptids) == "character") { 
  kept <- read.table(keptids)[,1]
} else {
  kept <- keptids
}

t <- t %>% filter(rowSums(  sapply(kept ,function(x) grepl(x, Term))  ) > 0)

p <- ggplot(t, aes(x = -log10(FDR), 
                   y = factor(Term, levels = rev(Term)),
                   shape = factor(Category), color = Percent, size = Count)) + 
  geom_point() + scale_x_continuous(n.breaks = 6) + theme_bw(base_size = 12) + 
  labs(y = "Term", shape = "Category", color = "% of Genes", size = "# of genes") + 
  scale_color_gradient(low="blue",high="red") + 
  scale_shape_manual(values = c(16,15)) + guides(shape = guide_legend(order=1), color = guide_colorbar(order=3), size=guide_legend(order=2, reverse = TRUE))

ggsave(plotname, p, width = 10,height = 4.5)
}