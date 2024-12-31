library(ggplot2)
library(readxl)

#' creates plots of up/downreugulated counts
#' 
#' @return nothing
#' 
#' @param DEFile1 DE level 0
#' @param DEFile2 DE level 1
#' @param order1 celltype order level 0
#' @param order2 celltype order level 1
#' @param padding horizontal padding in graphs
#' @param out1 output image level 0
#' @param out2 output image level 1
regulationplots <- function(DEFile1, DEFile2,
                            order1, order2,
                            padding,
                            out1, out2) {

filename <- DEFile1
sheetnames <- order1

bardata <- data.frame(Category=character(),
                      Celltype=character(), 
                      values=numeric(), 
                      values2=numeric(),
                      stringsAsFactors = FALSE) 

for (sheetname in sheetnames) {
  df <- read_xlsx(filename, sheetname, n_max = 1000)
  s <- df[df[,"FDR"] < .05 & !is.na(df[,"FDR"]),"avg_log2FC"] > 0
  if (sum(s) > 0){
    bardata <- bardata %>% add_row(Category = "upregulated", Celltype = sheetname, values = sum(s), values2 = sum(s))
    bardata <- bardata %>% add_row(Category = "downregulated", Celltype = sheetname, values = sum(1-s), values2 = -sum(1-s))
  }
}

p <- ggplot(data = bardata) + geom_bar(aes(x=factor(Celltype, levels = sheetnames),y=values2,fill=Category),stat="identity",position="identity") +
  geom_text(aes(x=Celltype,y=values2,label=abs(values2)),hjust = ifelse(bardata$values2 >= 0, -.5, 1.5)) +
  scale_y_continuous(labels=abs) + theme_bw() + xlab("") + ylab("Differentially Expressed Genes") +
  coord_flip() + geom_hline(yintercept = 0, linewidth=1) + ylim(c(min(bardata$values2)-padding,max(bardata$values2)+padding))
ggsave(out1, p, width = 6, height = 4)


filename <- DEFile2
sheetnames <- order2

bardata2 <- data.frame(Category=character(),
                       Celltype=character(),
                       Maintype=character(),
                       values=numeric(), 
                       values2=numeric(),
                       stringsAsFactors = FALSE) 

for (sheetname in sheetnames) {
  subst = substr(sheetname,1,1)
  if (subst == "G"){
    main <- "Granulosa"
  } else if (subst == "M") {
    main <- "Mesenchyme"
  } else {
    main <- sheetname
  }
  df <- read_xlsx(filename, sheetname, n_max = 1000)
  s <- df[df[,"FDR"] < .05 & !is.na(df[,"FDR"]),"avg_log2FC"] > 0
  if (length(s) > 0){
    bardata2 <- bardata2 %>% add_row(Category = "upregulated", Celltype = sheetname, Maintype = main, values = sum(s), values2 = sum(s))
    bardata2 <- bardata2 %>% add_row(Category = "downregulated", Celltype = sheetname, Maintype = main, values = sum(1-s), values2 = -sum(1-s))
  }
}

p2 <- ggplot(data = bardata2) + geom_bar(aes(x=factor(Celltype, levels = sheetnames),y=values2,fill=Category),stat="identity",position="identity") +
  geom_text(aes(x=Celltype,y=values2,label=abs(values2)),hjust = ifelse(bardata2$values2 > 0, -.5, 1.5)) +
  scale_y_continuous(labels=abs) + theme_bw() + xlab("") + ylab("Differentially Expressed Genes") +
  coord_flip() + geom_hline(yintercept = 0, linewidth=1) + geom_vline(xintercept = c(2.5,8.5),linewidth=2) +
  scale_fill_manual(values = c(upregulated = "#F8766D",downregulated = "#00BFC4")) +
  guides(color = guide_legend(reverse = TRUE), ) + ylim(c(min(bardata2$values2)-padding,max(bardata2$values2)+padding))
ggsave(out2, p2, width = 6, height = 4)

}