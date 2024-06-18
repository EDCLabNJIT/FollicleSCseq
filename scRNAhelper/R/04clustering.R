

#' Performs level 0 clustering
#' 
#' @return a Seurat Object
#' 
#' @param seuratObj input seruat object
#' @param seuratFile file for seurat rds if `seuratObj` is `NULL` 
#' @param resolution resolution for clustering
#' @param identities cluster renames
#' @param pc Number of principal components to use
#' @param umappc number of PC to use in UMAP. Morris uses 20
#' @param make.plots Make plots?
#' @param plotLocation Folder to put plots if `make.plots`
#' @param dotmarkersFile csv file with markers for dotplot
#' @param saverds saves as an rds if `TRUE`.
#' @param write.rdsfile path to save rds file to if `saverds == TRUE`
cluster <- function(seuratObj = NULL, seuratFile = "./saveddata/preprocessed_cells.rds",
                    resolution = 0.01, identities = c(),
                    PC = 15, umapPC = 20,
                    makeplots = TRUE, plotLocation = "./images/", dotmarkersFile = "./morrismarkers/morris_markers.csv",
                    saverds = TRUE, write.rdsfile = "./saveddata/clustered_cells.rds") {

#import data
if (is.null(seuratObj)){
  sc.cells <- readRDS(seuratFile)
} else {
  sc.cells <- seuratObj
}

sc.cells <- FindNeighbors(sc.cells, reduction = "pca", dims = 1:PC)

sc.cells <- RunUMAP(sc.cells, dims = 1:umapPC)
  
sc.cells <- FindClusters(sc.cells, resolution = resolution)

if (makeplots) {
  markers <- read.csv(dotmarkersFile, nrows = 10)
  p1 <- DotPlot(sc.cells, features = markers, cols= "RdBu") +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 40, hjust=1), axis.text = element_text(color = "black"),
          panel.grid.major = element_blank(), strip.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold"))
  ggsave(paste0(plotLocation,"precluster_dot.png"), plot = p1, width = 18, height = 4, dpi = 300, units = "in")
}

sc.cells$Level0 <- Idents(sc.cells)
Idents(sc.cells) <- "Level0"
sc.cells <- RenameIdents(sc.cells, identities)
sc.cells$Level0 <- Idents(sc.cells)

if (makeplots) {
  p2 <- DotPlot(sc.cells, features = markers, cols= "RdBu") +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = 40, hjust=1), axis.text = element_text(color = "black"),
          panel.grid.major = element_blank(), strip.text.x = element_text(face="bold"), axis.text.y = element_text(face="bold"))
  ggsave(paste0(plotLocation,"postcluster_dot.png"), plot = p2, width = 18, height = 4, dpi = 300, units = "in")
}

if (saverds) {
  saveRDS(sc.cells, write.rdsfile)
}

sc.cells
}