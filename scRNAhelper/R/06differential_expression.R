library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(clusterProfiler)
library(readxl)
library(writexl)

#' Performs differential expression
#' 
#' @return list of lists of `data.frame`s [level]
#' 
#' @param seuratObj input seruat object
#' @param seuratFile file for seurat rds if `seuratObj` is `NULL`
#' @param logfc.threshhold minimum logfc to be considered for DE.
#' @param pparFolder folder for ppar information
#' @param savedata saves as an xslx if `TRUE`.
#' @param DEFolder path to DE file if `savedata == TRUE
calcDE <- function(seuratObj = NULL, seuratFile = "./saveddata/mapped_cells.rds",
                   logfc.threshhold = 0.1, pparFolder = "./ppar/",
                   savedata = TRUE, DEFolder = "./data/") {

# import data
if (is.null(seuratObj)){
  sc.cells <- readRDS(seuratFile)
} else {
  sc.cells <- seuratObj
}

#ppar.pred <- read_xls(paste0(pparFolder,"Predicted PPAR target genes.xls")) %>% as.data.frame
#rownames(ppar.pred) <-  ppar.pred[,"Gene Symbol"]
#ppar.pred <- ppar.pred[,-1]
#colnames(ppar.pred) <- c("PPAR-P value", "Confidence Level")

#ppar.alpha <- read_xls(paste0(pparFolder,"PPAR alpha target genes.xls"))
#ppar.delta <- read_xls(paste0(pparFolder,"PPAR delta target genes.xls"))
#ppar.gamma <- read_xls(paste0(pparFolder,"PPAR gamma target genes.xls")) 
#ppar.all <- bind_rows(ppar.alpha, ppar.delta, ppar.gamma)
#ppar.all <- ppar.all[ppar.all$Species == "mouse",]

#ppar.all <- ppar.all %>% group_by(`Gene Symbol`) %>% summarise(across(c(`Tissue / Cell`,Regulation,`PPAR subtype`,`References (PMID)`),function(x) paste(x, collapse="/")))
#ppar.all <- as.data.frame(ppar.all)

#row.names(ppar.all) <- ppar.all[,"Gene Symbol"]
#ppar.all <- ppar.all[,-1]
#ppar.all <- bind_cols(ppar.pred, ppar.all[rownames(ppar.pred),])  

out = c()

for (level in -1:1) {
  # pseudobulks each identity
  if (level != -1){
  pseudo.sc <- AggregateExpression(sc.cells, assays = "RNA", return.seurat = TRUE, group.by = c("stim","orig.ident",paste0("Level",level)))
  } else {pseudo.sc <- AggregateExpression(sc.cells, assays = "RNA", return.seurat = TRUE, group.by = c("stim","orig.ident"))}
  
  # creates the identities of the bulks
  if (level != -1){
  preds <- pseudo.sc[[paste0("Level",level)]] %>% as.vector
  }
  stim <- pseudo.sc$stim %>% as.vector
  
  #combines identities into single identifier
  if (level != -1){
  pseudo.sc$celltype.stim <- paste(preds[[1]], stim, sep = "_")
  } else {pseudo.sc$celltype.stim <- paste0("_",pseudo.sc$stim)}
  
  Idents(pseudo.sc) <- "celltype.stim"
  
  # removes file to write to if it exists
  options(warn=2)
  xlsxfile <- paste0(DEFolder,"DE_level",level,".xlsx")
  if (file.exists(xlsxfile)) {
    file.remove(xlsxfile)
  }
  options(warn=0)
  
  # Only choosing cell types that occur in every batch
  if (level != -1){
  cell_types.count <- preds[[1]] %>% table
  cell_types <- names(cell_types.count)[cell_types.count == 6]
  } else {cell_types <- c("")}
  
  # use DESeq2 to calculate the bulk differential expression and write it to file
  sheets <- list()
  for (cell_type in cell_types){
    bulk.de <- FindMarkers(object = pseudo.sc, 
                           ident.1 = paste0(cell_type,"_MM"), ident.2 = paste0(cell_type,"_CTRL"),
                           test.use = "DESeq2", verbose = FALSE, features = VariableFeatures(pseudo.sc))
    
    bulk.de$FDR <- p.adjust(bulk.de$p_val, method ="BH")
    bulk.de <- tibble::rownames_to_column(bulk.de,"symbol")
    
    entrezid <- bitr(bulk.de$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    entrezid <- setNames(entrezid$ENTREZID, entrezid$SYMBOL)
    bulk.de$entrezid <- entrezid[bulk.de$symbol]
    bulk.de <- bulk.de[,c(1,8,2:7)]
    #bulk.de <- bind_cols(bulk.de, ppar.all[bulk.de$symbol,])
    sheets[[cell_type]] <- bulk.de
  }
  out[[level %>% as.character]] <- sheets
  write_xlsx(sheets, xlsxfile)
}

out
}