library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(gridExtra)
library(readxl)

#' Maps the cells onto the Morris atlas
#' 
#' @return a Seurat Object
#' 
#' @param seuratObj input seruat object
#' @param seuratFile file for seurat rds if `seuratObj` is `NULL`
#' @param atlasFile file for Morris atlas
#' @param make.plots Make plots?
#' @param plotLocation Folder to put plots if `make.plots`
#' @param atlas.rename .
#' @param umapPC number of PC to use in UMAP. Morris uses 20
#' @param PC Number of principal components to use
#' @param anchor.n.trees `n.trees` paramter for `FindTransferAnchors`
#' @param anchor.k.anchor `k.anchor` paramter for `FindTransferAnchors`
#' @param anchor.k.score `k.score` paramter for `FindTransferAnchors`
#' @param feature.folder folder with `.csv` files that contain markers to use
#' @param saverds saves as an rds if `TRUE`.
#' @param write.rdsfile path to save rds file to if `saverds == TRUE`
mapcells <- function(seuratObj = NULL, seuratFile = "./saveddata/preprocessed_cells.rds",
                       atlasFile = "./morrisdata/ovary_0.rds",
                       make.plots = TRUE, plotLocation = "./images/",
                       atlas.rename = c(), umapPC = 20,
                       PC = 15, anchor.n.trees = 250, anchor.k.anchor = 100, anchor.k.score = 200,
                       feature.folder = "./morrismarkers/",
                       saverds = TRUE, write.rdsfile = "./saveddata/mapped_cells.rds", projumapargs = list()) {

#import data
if (is.null(seuratObj)){
  sc.cells <- readRDS(seuratFile)
} else {
  sc.cells <- seuratObj
}
morris.atlas <- readRDS(atlasFile)

# renames identities
if (length(atlas.rename) > 0) {
  rename.df <- stack(atlas.rename)
  names(rename.df) <- c("new", "old")
  for (i in 1:nrow(rename.df)) {
    morris.atlas$Level1[morris.atlas$Level1 == rename.df[i,"old"]  ] <- rename.df[i,"new"]
  }
}

morris.atlas <- FindVariableFeatures(morris.atlas, selection.method = "vst", nfeatures = 5000, verbose = FALSE)

if (is.null(feature.folder)) {
  all.features <- NULL
} else {
  df.list = list()
  for (filename in list.files(feature.folder, pattern = ".csv")) {
    df.list[[filename]] <- read.csv(paste0(feature.folder,filename), nrows = 10)
  }
  all.features <- Reduce(cbind,df.list)
}

# Calculate UMAP on Morris Atlas
# It already contains a umap, but didnt return the model, which is necessary to map
morris.atlas <- RunUMAP(morris.atlas, dims = 1:20, reduction = "pca",reduction.name = "umap",
                        return.model = TRUE)

#calculate anchors and map cells
sc.anchors <- FindTransferAnchors(reference = morris.atlas, query = sc.cells, dims = 1:PC,
                                  features = unlist(all.features),
                                  reference.reduction = "pca", n.trees = anchor.n.trees,
                                  k.anchor = anchor.k.anchor, k.score = anchor.k.score)

sc.cells <- MapQuery(anchorset = sc.anchors, query = sc.cells, reference = morris.atlas,
                        reference.reduction = "pca", reduction.model = "umap", refdata = list(Level0 = "Level0", Level1 = "Level1"),
                     projectumap.args = projumapargs)

Idents(sc.cells) <- sc.cells$predicted.Level0


# plotting
if (make.plots) {
  col_palette <- c("dodgerblue2","#E31A1C","green4","#6A3D9A","#FF7F00","gold1","skyblue2",
                         "#FB9A99","palegreen2","#CAB2D6","#FDBF6F","gray70","khaki2","maroon","orchid1",
                         "deeppink1","steelblue4","darkturquoise","green1","yellow4","yellow3","darkorange4","brown" )
  
  celltypes.level0 <- morris.atlas$Level0 %>% table %>% sort(decreasing = TRUE) %>% names
  celltypes.level1 <- morris.atlas$Level1 %>% table %>% sort(decreasing = TRUE) %>% names
  colmap.level0 <- setNames(col_palette[1:length(celltypes.level0)], celltypes.level0)
  colmap.level1 <- setNames(col_palette[1:length(celltypes.level1)], celltypes.level1)
  
  p0 <- DimPlot(morris.atlas, reduction = "umap", group.by = "Level0",
                label = TRUE, repel=TRUE, label.box = TRUE, label.size = 3, cols=colmap.level0) + NoLegend()
  p1 <- DimPlot(morris.atlas, reduction = "umap", group.by = "Level1", shuffle = TRUE,
                label = TRUE, repel=TRUE, label.box = TRUE, label.size = 3, cols=colmap.level1) + NoLegend()
  p3 <- DimPlot(sc.cells, reduction = "ref.umap", group.by = "predicted.Level0",
                label = TRUE, repel=TRUE, label.box = TRUE, label.size = 3, cols=colmap.level0) + NoLegend()
  p4 <- DimPlot(sc.cells, reduction = "ref.umap", group.by = "predicted.Level1", shuffle = TRUE, seed = 2,
                label = TRUE, repel=TRUE, label.box = TRUE, label.size = 3, cols=colmap.level1) + NoLegend()
  p <- ggarrange(p0, p3, p1, p4, nrow=2, ncol=2)
  ggsave(paste0(plotLocation,"mapping_nolegend.png"), plot = p, width = 12, height = 12, dpi = 300, units = "in")
  
  q0 <- DimPlot(morris.atlas, reduction = "umap", group.by = "Level0",
                cols=colmap.level0) + guides(col = guide_legend(ncol = 1, override.aes = c(size = 3)))
  q1 <- DimPlot(morris.atlas, reduction = "umap", group.by = "Level1",
                cols=colmap.level1) + guides(col = guide_legend(ncol = 1, override.aes = c(size = 3)))
  legend1 <- cowplot::get_legend(q0)
  legend2 <- cowplot::get_legend(q1)
  q0 <- q0 + NoLegend()
  q1 <- q1 + NoLegend()
  q3 <- DimPlot(sc.cells, reduction = "ref.umap", group.by = "predicted.Level0",
                cols=colmap.level0) + NoLegend()
  q4 <- DimPlot(sc.cells, reduction = "ref.umap", group.by = "predicted.Level1",
                cols=colmap.level1) + NoLegend()
  
  q <- ggarrange(q0,q3,legend1, q1, q4, legend2, nrow = 2, ncol = 3, widths = c(1,1,.4)) +
          theme(panel.background = element_rect(fill="white"))
  ggsave(paste0(plotLocation,"mapping.png"), plot = q, width = 14.5, height = 12, dpi = 300, units = "in")
}

# save objects
saveRDS(sc.cells, write.rdsfile)

sc.cells
}