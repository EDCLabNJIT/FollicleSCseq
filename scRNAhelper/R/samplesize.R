library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(writexl)
library(tidyr)

#' Calculates Sample Size
#' 
#' @return a list of `data.frame`s
#' 
#' @param seuratObj input seruat object
#' @param seuratFile file for seurat rds if `seuratObj` is `NULL`
#' @param dataLocation folder to save sample size data
calcsamplesize <- function(seuratObj = NULL, seuratFile = "./saveddata/mapped_cells.rds",
                    dataLocation = "./data/") {

# import data
if (is.null(seuratObj)){
  sc.cells <- readRDS(seuratFile)
} else {
  sc.cells <- seuratObj
}

# split data into groups based on batch
# control is batches 2262 to 2264, treatment is 2265 to 2267
cells.split.by.orig <- SplitObject(sc.cells, split.by="orig.ident")
cells.split.by.stim <- SplitObject(sc.cells, split.by="stim")

# merges the columns of a data frame to a full data frame
merge.all <- function(x, y) {
  merge(x, y, all=TRUE, by="Cell Type")
}

# making and writing of data
sheets = list()
for (level in 0:1) {
  d = list()
  for (ident in as.character(2262:2267)) {
    d[[ident]] <- cells.split.by.orig[[ident]][[paste0("Level",level)]] %>% table %>% as.data.frame
    names(d[[ident]]) <- c("Cell Type",ident)
  }
  for (ident in c("CTRL", "MM")) {
    d[[ident]] <- cells.split.by.stim[[ident]][[paste0("Level",level)]] %>% table %>% as.data.frame
    names(d[[ident]]) <- c("Cell Type",ident)
  }
  
  
  df <- Reduce(merge.all, d) %>% replace(is.na(.), 0)
  rownames(df) <- df[,1]
  df <- df[,-1]
  total = data.frame(Total = colSums(df)) %>% t
  df <- rbind(df, total)
  df$total <- df$CTRL + df$MM
  df = tibble::rownames_to_column(df,var = "Cell Type")
  sheets[[paste0("Level",level)]] <- df
}
write_xlsx(sheets,paste0(dataLocation,"samplesize.xlsx"))

sheets
}