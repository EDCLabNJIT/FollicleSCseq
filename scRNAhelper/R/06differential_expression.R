library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readxl)
library(writexl)

#' Performs differential expression
#' 
#' @return list of lists of `data.frame`s [level]
#' 
#' @param seuratObj input seruat object
#' @param seuratFile file for seurat rds if `seuratObj` is `NULL`
#' @param logfc.threshhold minimum logfc to be considered for DE.
#' @param savedata saves as an xslx if `TRUE`.
#' @param DEFolder path to DE file if `savedata == TRUE
calcDE <- function(seuratObj = NULL, seuratFile = "./saveddata/mapped_cells.rds",
                   logfc.threshhold = 0.1,
                   savedata = TRUE, DEFolder = "./data/") {

# import data
if (is.null(seuratObj)){
  sc.cells <- readRDS(seuratFile)
} else {
  sc.cells <- seuratObj
}
  
out = c()

for (level in 0:1) {
  # pseudobulks each identity
  pseudo.sc <- AggregateExpression(sc.cells, assays = "RNA", return.seurat = TRUE, group.by = c("stim","orig.ident",paste0("Level",level)))
  
  # creates the identities of the bulks
  preds <- pseudo.sc[[paste0("Level",level)]] %>% as.vector
  stim <- pseudo.sc$stim %>% as.vector
  
  #combines identities into single identifier
  pseudo.sc$celltype.stim <- paste(preds[[1]], stim, sep = "_")
  Idents(pseudo.sc) <- "celltype.stim"
  
  # removes file to write to if it exists
  options(warn=2)
  xlsxfile <- paste0(DEFolder,"DE_level",level,".xlsx")
  if (file.exists(xlsxfile)) {
    file.remove(xlsxfile)
  }
  options(warn=0)
  
  # Only choosing cell types that occur in every batch
  cell_types.count <- preds[[1]] %>% table
  cell_types <- names(cell_types.count)[cell_types.count == 6]
  
  # use DESeq2 to calculate the bulk differential expression and write it to file
  sheets <- list()
  for (cell_type in cell_types){
    bulk.de <- FindMarkers(object = pseudo.sc, 
                           ident.1 = paste0(cell_type,"_MM"), ident.2 = paste0(cell_type,"_CTRL"),
                           test.use = "DESeq2", verbose = FALSE, features = VariableFeatures(pseudo.sc))
    
    bulk.de$FDR <- p.adjust(bulk.de$p_val, method ="BH")
    sheets[[cell_type]] <- tibble::rownames_to_column(bulk.de,"gene")
  }
  out[[level %>% as.character]] <- sheets
  write_xlsx(sheets, xlsxfile)
}

out
}