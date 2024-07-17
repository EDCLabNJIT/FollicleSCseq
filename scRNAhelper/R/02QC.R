library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scDblFinder)

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
#' @param findDbls use `scDblFinder`?
#' @param saverds saves as an rds if `TRUE`.
#' @param write.rdsfile path to save rds file to if `saverds == TRUE`
QC <- function(seuratObj = NULL, seuratFile = "./saveddata/merged_cells.rds", 
               make.plots = TRUE, plotLocation = "./images/", dataLoc = "./data/",
               nFeature.lower = 3500, nFeature.upper = 7000, percent.mt.upper = 5, nCount.rna.upper = 60000,
               findDbls = FALSE,
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

if (findDbls) {
  keptcells = c()
  for (i in 2262:2267) {
    count_matrix = sc.filter@assays$RNA@layers[[paste0("counts.",i)]]
    is_singlet <- scDblFinder(count_matrix, returnType = "table")$class == "singlet"
    keptcells <- append(keptcells, names(Idents(sc.filter))[Idents(sc.filter) == i][is_singlet])
  }
  subset(sc.filter, cells = keptcells)
}

if (make.plots) {
  # Supplementary image postQCViolin
  p <- VlnPlot(sc.filter, features = c("nCount_RNA","nFeature_RNA","percent_mt"), alpha = .25)
  ggsave(paste0(plotLocation,"postQCViolin.png"), plot = p, width = 10, height = 4, dpi = 300, units = "in")
  
  # QC_test
  testloc = paste0(dataLoc,"QC_test.txt")
  if (file.exists(testloc)) {
    file.remove(testloc)
  }
  precount <- sc.ser$orig.ident %>% table
  postcount <- sc.filter$orig.ident %>% table
  proportions <- postcount / precount
  
  
  
  sink(testloc)
  print("cells before\n")
  print(precount)
  print("\ncells after\n")
  print(postcount)
  print("\nproportions\n")
  print(proportions)
  print("\n")
  print(prop.test(cbind(postcount, precount-postcount)))
  print(t.test(proportions[1:3], proportions[4:6]))
  sink()
}

if (saverds) {
  # save filtered cells
  saveRDS(sc.filter, write.rdsfile)
}

sc.filter
}