library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#' Performs Quality Control
#' 
#' @return a Seurat Object
#' 
#' @param seuratObj input seruat object
#' @param seuratFile file for seurat rds if `seuratObj` is `NULL`
#' @param make.plots Make plots?
#' @param plotLocation Folder to put plots if `make.plots`
#' @param nFeature.lower removes cells with `nFeature < nFeature.lower`
#' @param nFeature.upper removes cells with `nFeature > nFeature.upper`
#' @param percent.mt.upper removes cells with `percemt_mt > percent.mt.upper`
#' @param nCount.rna.upper removes cells with `nCount_RNA > nCount.rna.upper`
#' @param saverds saves as an rds if `TRUE`.
#' @param write.rdsfile path to save rds file to if `saverds == TRUE`
QC <- function(seuratObj = NULL, seuratFile = "./saveddata/merged_cells.rds", 
               make.plots = TRUE, plotLocation = "./images/",
               nFeature.lower = 3500, nFeature.upper = 7000, percent.mt.upper = 5, nCount.rna.upper = 60000,
               saverds = TRUE, write.rdsfile = "./saveddata/filtered_cells.rds") {

# load merged cells
if (is.null(seuratObj)){
  sc.ser <- readRDS(seuratFile)
} else {
  sc.ser <- seuratObj
}

# find filter range
sc.ser[["percent_mt"]] <- PercentageFeatureSet(sc.ser, pattern = "^mt-")

if (make.plots) {
  # Supplementary image preQCViolin
  p <- VlnPlot(sc.ser, features = c("nCount_RNA","nFeature_RNA","percent_mt"), alpha = .25) 
  ggsave(paste0(plotLocation,"preQCViolin.png"), plot = p, width = 10, height = 4, dpi = 300, units = "in")
}

# filter data
sc.filter <- subset(sc.ser, subset = nFeature_RNA > nFeature.lower & nFeature_RNA < nFeature.upper & percent_mt < percent.mt.upper & nCount_RNA < nCount.rna.upper)

if (make.plots) {
  # Supplementary image postQCViolin
  p <- VlnPlot(sc.filter, features = c("nCount_RNA","nFeature_RNA","percent_mt"), alpha = .25)
  ggsave(paste0(plotLocation,"postQCViolin.png"), plot = p, width = 10, height = 4, dpi = 300, units = "in")
}

if (saverds) {
  # save filtered cells
  saveRDS(sc.filter, write.rdsfile)
}

sc.filter
}