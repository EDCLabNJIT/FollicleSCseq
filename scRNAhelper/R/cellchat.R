library(CellChat)
library(patchwork)
library(ggplot2)
library(writexl)
library(dplyr)
library(keras)
options(stringsAsFactors = FALSE)
future::plan("multisession", workers = 4) # do parallel


prepareCCC <- function() {
CellChatDB.use <- CellChatDB.mouse
seurat_object <- readRDS("saveddata/mapped_cells.rds")
data.input <- seurat_object[["RNA"]]$data # normalized data matrix
labels <- Idents(seurat_object)
meta <- data.frame(labels = labels, samples = seurat_object$orig.ident,row.names = names(labels)) # create a dataframe of the cell labels
ctrlCells <- seurat_object$stim == "CTRL"

meta.ctrl <- meta[ctrlCells,]
data.in.ctrl <- data.input[,ctrlCells]

meta.mm <- meta[!ctrlCells,]
data.in.mm <- data.input[,!ctrlCells]

cellchat.ctrl <- createCellChat(object = data.in.ctrl, meta = meta.ctrl, group.by = "labels")
cellchat.ctrl@DB <- CellChatDB.use
cellchat.ctrl <- subsetData(cellchat.ctrl) # This step is necessary even if using the whole database
cellchat.ctrl <- identifyOverExpressedGenes(cellchat.ctrl)
cellchat.ctrl <- identifyOverExpressedInteractions(cellchat.ctrl)
cellchat.ctrl <- computeCommunProb(cellchat.ctrl, type = "triMean")
cellchat.ctrl <- filterCommunication(cellchat.ctrl, min.cells = 10)
cellchat.ctrl <- computeCommunProbPathway(cellchat.ctrl)
cellchat.ctrl <- aggregateNet(cellchat.ctrl)
cellchat.ctrl <- netAnalysis_computeCentrality(cellchat.ctrl)

cellchat.mm <- createCellChat(object = data.in.mm, meta = meta.mm, group.by = "labels")
cellchat.mm@DB <- CellChatDB.use
cellchat.mm <- subsetData(cellchat.mm) # This step is necessary even if using the whole database
cellchat.mm <- identifyOverExpressedGenes(cellchat.mm)
cellchat.mm <- identifyOverExpressedInteractions(cellchat.mm)
cellchat.mm <- computeCommunProb(cellchat.mm, type = "triMean")
cellchat.mm <- filterCommunication(cellchat.mm, min.cells = 10)
cellchat.mm <- computeCommunProbPathway(cellchat.mm)
cellchat.mm <- aggregateNet(cellchat.mm)
cellchat.mm <- netAnalysis_computeCentrality(cellchat.mm)

object.list <- list(CTRL = cellchat.ctrl, MM = cellchat.mm)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

saveRDS(cellchat.mm, "saveddata/cellchatmm.rds")
saveRDS(cellchat.ctrl, "saveddata/cellchatctrl.rds")
saveRDS(cellchat, "saveddata/cellchat.rds")
}








doCCC <- function() {
cellchat.mm <- readRDS("saveddata/cellchatmm.rds")
cellchat.ctrl <- readRDS("saveddata/cellchatctrl.rds")
cellchat <- readRDS("saveddata/cellchat.rds")

object.list <- list(CTRL = cellchat.ctrl, MM = cellchat.mm)



###############



png("images/interactionslvl1_nolabels.png", width=8, height=4, units="in", res = 500)
par(mfrow = c(1,2), xpd=TRUE)
gg1 <- netVisual_circle(cellchat@net$CTRL$count, weight.scale = TRUE, label.edge= FALSE, title.name = "Number of interactions - CTRL", margin = 0, vertex.label.cex = .Machine$double.xmin)

gg2 <- netVisual_circle(cellchat@net$MM$count, weight.scale = TRUE, label.edge= FALSE, title.name = "Number of interactions - MM", margin = 0, vertex.label.cex = .Machine$double.xmin)
dev.off()



###############



group.cellType <- c("Endothelium", "Immune", "Epithelium", rep("Granulosa",7), rep("Mesenchyme",6))
names(group.cellType) <- levels(cellchat.ctrl@idents)
pathways = c("TGFb","WNT","BMP")
for (path in pathways) {
  png(paste0("images/CCC_",path,".png"), width=16, height=8, units="in", res = 500)
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_chord_cell(cellchat.ctrl, signaling = path, group = group.cellType, big.gap = 8, small.gap = 1.5, title.name = "CONTROL")
  netVisual_chord_cell(cellchat.mm, signaling = path, group = group.cellType, big.gap = 8, small.gap = 1.5, title.name = "MM")
  dev.off()
}



###############



gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight", )
#> Do heatmap based on a merged object
png("images/heatmap.png", width=9, height=5, units="in", res = 500)
gg1 + gg2
dev.off()



###############



# coarse celltype
group.cellType <- c("Endothelium", "Immune", "Epithelium", rep("Granulosa",7), rep("Mesenchyme",6)) %>% as.factor
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))

png("images/interactions.png", width=8, height=4, units="in", res = 500)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = TRUE, label.edge= TRUE, edge.weight.max = weight.max[3], edge.width.max = 12, arrow.size = .5,title.name = paste0("Number of interactions - ", names(object.list)[i]), margin = 0, )
}
dev.off()

png("images/interactions_nolabels.png", width=8, height=4, units="in", res = 500)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = TRUE, label.edge= FALSE, vertex.label.cex = .Machine$double.xmin, edge.weight.max = weight.max[3], edge.width.max = 12, arrow.size = .5,title.name = paste0("Number of interactions - ", names(object.list)[i]), margin = 0, )
}
dev.off()

ctrldat <- sapply(cellchat@net$CTRL[4:7], as.data.frame)
stimdat <- sapply(cellchat@net$MM[4:7], as.data.frame)
diffdat <- mapply(function(x,y) x-y, cellchat@net$MM[4:7], cellchat@net$CTRL[4:7]) %>%
              sapply(as.data.frame)
write_xlsx(ctrldat, "data/CCctrldata.xlsx")
write_xlsx(stimdat, "data/CCmmdata.xlsx")
write_xlsx(diffdat, "data/CCdiffdata.xlsx")

groups <- list()
for (stim in c("CTRL","MM")){
  dn = dimnames(cellchat@net[[stim]]$pval)
  tab = cellchat@net[[stim]]$pval %>% array_reshape(c(16^2,length(dn[[3]])), "F") %>% t
  colgrid <- expand_grid(dn[[1]],dn[[2]])
  rownames(tab) <- dn[[3]]
  colnames(tab) <- sapply(1:256, function(i) paste(colgrid[i,2],colgrid[i,1],sep=" -> "))
  tab <- tab %>% as.data.frame
  tab <- cbind(rownames(tab),tab)
  colnames(tab)[1] <- "Pathway"
  groups[[stim]] <- tab
}
write_xlsx(groups,"data/CCpval.xlsx")
}