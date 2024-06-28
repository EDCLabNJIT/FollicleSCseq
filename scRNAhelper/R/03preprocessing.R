library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


#' Performs Preprocessing
#' 
#' @return a Seurat Object
#' 
#' @param seuratObj input seruat object
#' @param seuratFile file for seurat rds if `seuratObj` is `NULL`
#' @param make.plots Make plots?
#' @param plotLocation Folder to put plots if `make.plots`
#' @param nfeatures features to select
#' @param regress.mt regress out mitochondrial percentage when scaling?
#' @param saverds saves as an rds if `TRUE`.
#' @param write.rdsfile path to save rds file to if `saverds == TRUE`
preprocess <- function(seuratObj = NULL, seuratFile = "./saveddata/filtered_cells.rds", 
               make.plots = TRUE, plotLocation = "./images/",
               nfeatures = 2000, regress.mt = TRUE,
               saverds = TRUE, write.rdsfile = "./saveddata/preprocessed_cells.rds") {

# load filtered cells
if (is.null(seuratObj)){
  sc.filtered <- readRDS(seuratFile)
} else {
  sc.filtered <- seuratObj
}

# normalize data
sc.filtered <- NormalizeData(sc.filtered, verbose = FALSE)

# feature selecting
sc.filtered <- FindVariableFeatures(sc.filtered, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)

if (make.plots) {
  top10 <- head(VariableFeatures(sc.filtered), 10)
  plot1 <- VariableFeaturePlot(sc.filtered) + theme_classic()
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  ggsave(paste0(plotLocation,"featureSelection.png"), plot = plot2, width = 8, height = 5, dpi = 300, units = "in")
}

# scaling
if (regress.mt) {
  sc.filtered <- ScaleData(sc.filtered, vars.to.regress = "percent_mt")
} else {
  sc.filtered <- ScaleData(sc.filtered)
}

# PCA
sc.filtered <- RunPCA(sc.filtered, features = VariableFeatures(object = sc.filtered))

if (make.plots) {
  data <- tibble(PC = 1:30,
                 stdev = sc.filtered@reductions$pca@stdev[1:30],
                 expl = cumsum(sc.filtered@reductions$pca@stdev[1:30] ^2) / sc.filtered@reductions$pca@misc$total.variance)
  p <- ggplot(data, aes(PC, stdev)) +
    geom_point(size = 2) +
    theme_bw() +
    geom_line(aes(y = expl*30), linewidth = 2, col = "blue") +
    scale_y_continuous("Standard Deviation", sec.axis =  sec_axis(~ ./30,name = "Proportion Explained")) +
    theme(axis.title.y.right = element_text(color = "blue"), text = element_text(size=14))
  
  ggsave(paste0(plotLocation,"elbow.png"), plot = p, width = 5, height = 5, dpi = 300, units = "in")
}

# integration
sc.filtered <- IntegrateLayers(object = sc.filtered, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                            verbose = FALSE)
sc.filtered[["RNA"]] <- JoinLayers(sc.filtered[["RNA"]])

# save data
if (saverds) {
  saveRDS(sc.filtered, write.rdsfile)
}

sc.filtered
}