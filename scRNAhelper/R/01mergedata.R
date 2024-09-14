library(dplyr)
library(Seurat)
library(patchwork)
library(hdf5r)

#' Grabs 10X data from file, merges into single Seurat object
#' 
#' 
#' @return a Seurat Object
#' 
#' @param min.features min.features for `CreateSeuratObject` function
#' @param min.cells min.cells for `CreateSeuratObject` function
#' @param saverds saves as an rds if `TRUE`.
#' @param write.rdsfile path to save rds file to if `saverds == TRUE`
mergedata <- function(min.features = 1000, min.cells = 3,
                      raw_matrix = "./raw_feature_bc_matrix.h5",
                      saverds = TRUE, write.rdsfile = "./saveddata/merged_cells.rds") {

# Grab datasets from file
#data.in <- list()
#for (num in as.character(2262:2267)){
#  data.in[[num]] <- Read10X(paste0(TenX.folder,num))
#}
data <- Read10X_h5(raw_matrix)

# Convert datasets to Seurat objects
# note: want 'Gene Expression'
idents <- substring(data$`Gene Expression`@Dimnames[[2]],18,18) %>% factor
idents <- (as.numeric(levels(idents))[idents] + 2261) %>% as.character
df <- data.frame(orig.ident = idents)
rownames(df) <- data$`Gene Expression`@Dimnames[[2]]

ser.obj <- CreateSeuratObject(data[["Gene Expression"]], meta.data =df,
                              min.features = min.features, min.cells = min.cells, proj="scSeq")

ser.objs <- SplitObject(ser.obj, split.by = "orig.ident")

#ser.objs <- list()
#for (num in as.character(2262:2267)){
#  ser.objs[[num]] <- CreateSeuratObject(data.in[[num]][["Gene Expression"]],
#                                        min.features = min.features, min.cells = min.cells, project = num)
#}

# merge cells together
merge.cells <- merge(ser.objs[[1]], y =ser.objs[-1], add.cell.ids = as.character(2262:2267), proj="scSeq")
Layers(merge.cells) <- as.character(2262:2267)

# label STIM and CTRL
merge.cells$stim <- ifelse(ser.obj$orig.ident %in% as.character(2262:2264),"CTRL","MM")

# save merged cells
if (saverds) {
  saveRDS(merge.cells, write.rdsfile)
}

merge.cells
}
