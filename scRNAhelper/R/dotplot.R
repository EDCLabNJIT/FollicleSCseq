library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(stringr)

#' Creates DotPlots
#' 
#' @return nothing
#' 
#' @param seuratObj input seruat object
#' @param seuratFile file for seurat rds if `seuratObj` is `NULL`
#' @param plotLocation folder to plot images
#' @param morris.gran.markers csv where morris top granulosa markers are stored
#' @param morris.mes.merkers csv where morris top mesenchyme markers are stored
#' @param nGenes number of top genes for each cluster to plot
makedot <- function(seuratObj = NULL, seuratFile = "./saveddata/mapped_cells.rds",
                    plotLocation = "./images/", morris.gran.markers = "./morrismarkers/morris_gran_markers.csv",
                    morris.mes.merkers = "./morrismarkers/morris_mes_markers.csv",
                    nGenes = 5) {

# load data
if (is.null(seuratObj)){
  sc.cells <- readRDS(seuratFile)
} else {
  sc.cells <- seuratObj
}

Idents(sc.cells) <- sc.cells$predicted.Level1

sc.cells@active.ident <- factor(sc.cells@active.ident, Idents(sc.cells) %>% table %>% sort(decreasing = FALSE) %>% names)

gran.markers <- read.csv(morris.gran.markers, nrows = nGenes)
gran.types <- gran.markers %>% colnames

mes.markers <- read.csv(morris.mes.merkers, nrows = nGenes)
mes.types <- mes.markers %>% colnames


sc.cell.types <- sc.cells$predicted.Level1 %>% table %>% sort(decreasing = TRUE) %>% names
sc.mes.types <- str_subset(sc.cell.types, "^M+")
sc.gran.types <- str_subset(sc.cell.types, "^GC+")


p1 <- DotPlot(sc.cells, gran.markers, idents=sc.gran.types, cols= "RdBu") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 40, hjust=1), axis.text = element_text(color = "black"),
        panel.grid.major = element_blank(), strip.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold"))

p2 <- DotPlot(sc.cells, mes.markers, idents=sc.mes.types, cols= "RdBu") +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 40, hjust=1), axis.text = element_text(color = "black"),
        panel.grid.major = element_blank(), strip.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold"))

ggsave(paste0(plotLocation,"gran_dot.png"), plot = p1, width = 16, height = 4, dpi = 300, units = "in")
ggsave(paste0(plotLocation,"mes_dot.png"), plot = p2, width = 16, height = 4, dpi = 300, units = "in")


}