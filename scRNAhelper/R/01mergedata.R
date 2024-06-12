library(dplyr)
library(Seurat)
library(patchwork)

#' Grabs 10X data from file, merges into single Seurat object
#' 
#' @details
#' 10X files must be located in `TenX.folder`/num
#' 
#' @return a Seurat Object
#' 
#' @param min.features min.features for `CreateSeuratObject` function
#' @param min.cells min.cells for `CreateSeuratObject` function
#' @param saverds saves as an rds if `TRUE`.
#' @param write.rdsfile path to save rds file to if `saverds == TRUE`
mergedata <- function(min.features = 1000, min.cells = 3,
                      TenX.folder = "./10X/",
                      saverds = TRUE, write.rdsfile = "./saveddata/merged_cells.rds") {

# Grab datasets from file
data.in <- list()
for (num in as.character(2262:2267)){
  data.in[[num]] <- Read10X(paste0(TenX.folder,num))
}

# Convert datasets to Seurat objects
# note: want 'Gene Expression'
# min features: 1000
# min cells: 3
ser.objs <- list()
for (num in as.character(2262:2267)){
  ser.objs[[num]] <- CreateSeuratObject(data.in[[num]][["Gene Expression"]],
                                        min.features = min.features, min.cells = min.cells, project = num)
}

# merge cells together
merge.cells <- merge(ser.objs[[1]], y =ser.objs[-1], add.cell.ids = as.character(2262:2267), proj="scSeq")

# label STIM and CTRL
merge.cells$stim <- ifelse(merge.cells$orig.ident %in% as.character(2262:2264),"CTRL","MM")

# save merged cells
if (saverds) {
  saveRDS(merge.cells, write.rdsfile)
}

merge.cells
}