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
#' @return a Splitted Seurat Object
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
mapcells <- function(seuratObj = NULL, seuratFile = "./saveddata/clustered_cells.rds",
                       atlasFile = "./morrisdata/ovary_0.rds",
                       make.plots = TRUE, plotLocation = "./images/",
                       atlas.rename = c(), umapPC = 20,
                       PC = 15, anchor.n.trees = 50, anchor.k.anchor = 15, anchor.k.score = 100,
                       marker.folder = "./morrismarkers/",
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

# gets markers
markers <- list()
for (filename in list.files("morrismarkers/", pattern = ".csv")) {
  markers[[filename]] <- read.csv(paste0("morrismarkers/",filename), nrows = 10)
}

# Calculate UMAP on Morris Atlas
# It already contains a umap, but didnt return the model, which is necessary to map
morris.atlas <- RunUMAP(morris.atlas, dims = 1:umapPC, reduction = "pca",reduction.name = "umap",
                        return.model = TRUE)


# projects umap onto sc.cells
sc.cells <- ProjectUMAP(reference = morris.atlas, query = sc.cells, query.dims = 1:umapPC, reference.dims = 1:umapPC,
                        query.reduction = "pca", reference.reduction = "pca", reduction.model = "umap", reduction.name = "umap2")




# split sc and atlas by level 0
sc.split = SplitObject(sc.cells, split.by = "Level0")
atlas.split = SplitObject(morris.atlas, split.by = "Level0")

sc.mes = sc.split$Mesenchyme
sc.gran = sc.split$Granulosa

atlas.mes <- atlas.split$Mesenchyme
atlas.gran <- atlas.split$Granulosa


# granlosa level 1
gran.anchors <- FindTransferAnchors(reference = atlas.gran, query = sc.gran, dims = 1:PC,
                                    reference.reduction = "pca", k.anchor = anchor.k.anchor, k.score = anchor.k.score, n.trees = anchor.n.trees)

sc.gran <- MapQuery(anchorset = gran.anchors, query = sc.gran, reference = atlas.gran,
                    reference.reduction = "pca", reduction.model = "umap", refdata = list(Level1 = "Level1"))
sc.gran$Level1 <- sc.gran$predicted.Level1
Idents(sc.gran) <- "Level1"

# mesenchyme level 1
mes.anchors <- FindTransferAnchors(reference =atlas.mes, query = sc.mes, dims = 1:PC,
                                    reference.reduction = "pca", k.anchor = anchor.k.anchor, k.score = anchor.k.score, n.trees = anchor.n.trees)

sc.mes <- MapQuery(anchorset = mes.anchors, query = sc.mes, reference = atlas.mes,
                               reference.reduction = "pca", reduction.model = "umap", refdata = list(Level1 = "Level1"))
sc.mes$Level1 <- sc.mes$predicted.Level1
Idents(sc.mes) <- "Level1"

# moves Level1 objects up to original object
Idents(sc.cells) <- "Level0"
cells <- Idents(sc.cells)
levels(cells) <- c(levels(cells), levels(sc.gran), levels(sc.mes))
cells[cells == "Granulosa"] <- Idents(sc.gran)
cells[cells == "Mesenchyme"] <- Idents(sc.mes)
sc.cells$Level1 <- cells
Idents(sc.cells) <- "Level1"

# dot plots
if (make.plots) {
  dotgran <- DotPlot(sc.gran, features = markers$morris_gran_markers.csv, group.by = "Level1", cols= "RdBu") +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 40, hjust=1), axis.text = element_text(color = "black"),
          panel.grid.major = element_blank(), strip.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold"))
  
  dotmes <- DotPlot(sc.mes, features = markers$morris_mes_markers.csv, group.by = "Level1", cols= "RdBu") +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 40, hjust=1), axis.text = element_text(color = "black"),
          panel.grid.major = element_blank(), strip.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold"))
  
  ggsave(paste0(plotLocation,"gran_dot.png"), plot = dotgran, width = 18, height = 4, dpi = 300, units = "in")
  ggsave(paste0(plotLocation,"mes_dot.png"), plot = dotmes, width = 18, height = 4, dpi = 300, units = "in")
}


# plotting
if (make.plots) {
  col_palette <- c("dodgerblue2","#E31A1C","green4","#6A3D9A","#FF7F00","gold1","skyblue2",
                         "#FB9A99","palegreen2","#CAB2D6","#FDBF6F","gray70","khaki2","maroon","orchid1",
                         "deeppink1","steelblue4","darkturquoise","green1","yellow4","yellow3","darkorange4","brown" )
  
  celltypes.level0 <- sc.cells$Level0 %>% table %>% sort(decreasing = TRUE) %>% names
  celltypes.level1 <- sc.cells$Level1 %>% table %>% sort(decreasing = TRUE) %>% names
  colmap.level0 <- setNames(col_palette[1:length(celltypes.level0)], celltypes.level0)
  colmap.level1 <- setNames(col_palette[1:length(celltypes.level1)], celltypes.level1)
  
  #p0 <- DimPlot(morris.atlas, reduction = "umap", group.by = "Level0",
  #              label = TRUE, repel=TRUE, label.box = TRUE, label.size = 3, cols=colmap.level0) + NoLegend()
  #p1 <- DimPlot(morris.atlas, reduction = "umap", group.by = "Level1", shuffle = TRUE,
  #              label = TRUE, repel=TRUE, label.box = TRUE, label.size = 3, cols=colmap.level1) + NoLegend()
  p3 <- DimPlot(sc.cells, reduction = "umap", group.by = "Level0",
                label = TRUE, repel=TRUE, label.box = TRUE, label.size = 3, cols=colmap.level0) + NoLegend()
  p4 <- DimPlot(sc.cells, reduction = "umap", group.by = "Level1", shuffle = TRUE, seed = 2,
                label = TRUE, repel=TRUE, label.box = TRUE, label.size = 3, cols=colmap.level1) + NoLegend()
  p <- ggarrange(p3, p4, nrow=2, ncol=1)
  ggsave(paste0(plotLocation,"mapping_nolegend.png"), plot = p, width = 6, height = 12, dpi = 300, units = "in")
  
  #q0 <- DimPlot(morris.atlas, reduction = "umap", group.by = "Level0",
  #              cols=colmap.level0) + guides(col = guide_legend(ncol = 1, override.aes = c(size = 3)))
  #q1 <- DimPlot(morris.atlas, reduction = "umap", group.by = "Level1",
  #              cols=colmap.level1) + guides(col = guide_legend(ncol = 1, override.aes = c(size = 3)))
  #legend1 <- cowplot::get_legend(q0)
  #legend2 <- cowplot::get_legend(q1)
  #q0 <- q0 + NoLegend()
  #q1 <- q1 + NoLegend()
  q3 <- DimPlot(sc.cells, reduction = "umap", group.by = "Level0",
                cols=colmap.level0)# + NoLegend()
  q4 <- DimPlot(sc.cells, reduction = "umap", group.by = "Level1",
                cols=colmap.level1)# + NoLegend()
  legend1 <- cowplot::get_legend(q3)
  legend2 <- cowplot::get_legend(q4)
  q3 <- q3 + NoLegend()
  q4 <- q4 + NoLegend()
  q <- ggarrange(q3, legend1, q4, legend2, nrow = 2, ncol = 2, widths = c(1,.4)) +
          theme(panel.background = element_rect(fill="white"))
  ggsave(paste0(plotLocation,"mapping.png"), plot = q, width = 8, height = 12, dpi = 300, units = "in")
}
  
# save objects
if (saverds) {
  saveRDS(sc.cells, write.rdsfile)
}

sc.cells
}
