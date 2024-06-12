library(clusterProfiler)
# If you use clusterProfiler in published research, please cite:
# T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141
library(org.Mm.eg.db)
library(AnnotationDbi)
library(readxl)
library(writexl)
library(ggplot2)

#' Performs enrichment analysis
#' 
#' @return a list of data.frames
#' 
#' @param df input data frame
#' @param DEFolder folder for DE data if `df` is `NULL`
#' @param level level to search at (0 or 1) if `df` is `NULL`
#' @param cell_type cell type to search as labelled by name if `df` is `NULL`
#' @param de.q.val maximum adjusted p value for differentially expressed genes
#' @param log2fold.val minimum log2fold change needed for DE genes (non-negative)
#' @param sign.val -1 for values less than negative. 0 for greater than or less than negative. 1 for values greater
#' @param savedata saves as an xslx if `TRUE`.
#' @param EAFolder path to EA file if `savedata == TRUE`
enrich <- function(df = NULL, DEFolder = "./data/",
                   de.q.val = .05, log2fold.val = 0, sign.val = 0,
                   level = 0, sheetname = "Granulosa",
                   savedata = TRUE, EAFolder = "./data/EA") {

log2fold.val = abs(log2fold.val)

if (is.null(df)) {
  # read data
  filename <- paste0(DEFolder,"DE_level",level,".xlsx")
  df <- read_xlsx(filename, sheetname, n_max = 1000)
}

# filter genes
q.filter <- df$p_val_adj < de.q.val
log2.filter <- q.filter & FALSE # set to all FALSE
if (sign.val <= 0) {
  log2.filter <- log2.filter | df$avg_log2FC <= -log2fold.val
}
if (sign.val >= 0) {
  log2.filter <- log2.filter | df$avg_log2FC >= log2fold.val
}
  
genes <- df[q.filter & log2.filter,"gene"]$gene

# enrich genes
go.bp.en <- enrichGO(gene     = genes,
                  OrgDb    = "org.Mm.eg.db",
                  ont      = "BP",
                  keyType  = "SYMBOL")
go.bp.df <- as.data.frame(go.bp.en)

kegg.en <- enrichKEGG(gene     = bitr(genes, "ALIAS", "UNIPROT", "org.Mm.eg.db")$UNIPROT,
                      organism = "mmu",
                      keyType  = "uniprot")
kegg.df <- as.data.frame(kegg.en)


# write enriched data.
xlsxfile <- paste0(EAFolder,sheetname,"_",
                   de.q.val,"_",log2fold.val,
                   ifelse(sign.val < 0,"-",""),ifelse(sign.val > 0,"+",""),
                   ".xlsx")
if (file.exists(xlsxfile)) {
  file.remove(xlsxfile)
}

output = list(GO.bp = go.bp.df, KEGG = kegg.df)
write_xlsx(output, xlsxfile)

output
}
