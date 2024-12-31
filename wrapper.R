suppressWarnings(devtools::document("./scRNAhelper/"))
sc.cells <- mergedata(min.features = 1000, min.cells = 3,
                      raw_matrix = "./raw_feature_bc_matrix.h5",
                      saverds = FALSE, write.rdsfile = "./saveddata/merged_cells.rds")

sc.cells <- QC(Nsc.cells, make.plots = TRUE, plotLocation = "./images/",
               nFeature.lower = 3500, nFeature.upper = Inf, percent.mt.upper = 5, nCount.rna.upper = 100000, findDbls = TRUE,
               saverds = FALSE, write.rdsfile = "./saveddata/filtered_cells.rds")

sc.cells <- preprocess(sc.cells, make.plots = TRUE, plotLocation = "./images/",
                       nfeatures = 5000,
                       saverds = FALSE, write.rdsfile = "./saveddata/preprocessed_cells.rds")

sc.cells <- cluster(seuratObj = NULL, seuratFile = "./saveddata/preprocessed_cells.rds", PC = 15,
                    resolution = 0.01, identities = c("0" = "Granulosa","1" = "Mesenchyme","2" = "Endothelium","3" = "Immune", "4" = "Epithelium"),
                    makeplots = TRUE, plotLocation = "./images/", dotmarkers = "./morrismarkers/morris_markers.csv",
                    saverds = FALSE, write.rdsfile = "./saveddata/clustered_cells.rds")

#                  old name                  new name
atlas.rename <- c("GC_Estrous"            = "GC_Luteinizing",
                  "GC_CL_Lytic"           = "GC_Regressing CL",
                  "GC_CL_Active"          = "GC_Active CL",
                  "GC_Antral"             = "GC_Mural",
                  "GC_Preantral"          = "GC_Cumulus",
                  "M_Immature Theca"      = "M_Early Theca",
                  "M_Dividing Mesenchyme" = "M_Early Theca", 
                  "M_Cortical Stroma"     = "M_Steroidogenic Stroma",
                  "M_Medullary Stroma"    = "M_Fibroblast-like Stroma")


sc.cells <- mapcells(seuratObj = NULL, seuratFile = "./saveddata/clustered_cells.rds",
                     atlasFile = "./ovary_0.rds",
                     make.plots = TRUE, plotLocation = "./images/",
                     atlas.rename = atlas.rename, PC = 20, umapPC = 20,
                     anchor.n.trees = 50, anchor.k.anchor = 5, anchor.k.score = 30,
                     saverds = TRUE, write.rdsfile = "./saveddata/mapped_cells.rds",
                     marker.folder = "./morrismarkers/", 
                     projumapargs = list(annoy.metric = "cosine"))

DE.data <- calcDE(sc.cells, 
                  savedata = TRUE, DEFolder = "./data/")

# Enrichment done with DAVID
#                   level sheetname  p/qval reg log2fold cutoff
#enrichvals = list(list(1, "GC-Mural" , .05, 0 , 0),
#                  list(1, "GC-Mural" , .05, 1 , 0),
#                  list(1, "GC-Mural" , .05, -1, 0),
#                  list(1, "GC-Mural" , .05, 0 , .3),
#                  list(1, "GC-Mural" , .05, 1 , .3),
#                  list(1, "GC-Mural" , .05, -1, .3),
#                  list(0, "Granulosa", .05, 0 , 0),
#                  list(0, "Granulosa", .05, 1 , 0),
#                  list(0, "Granulosa", .05, -1 , 0))
#for (eval in enrichvals) {
#  EA.data <- enrich(p.val = "FDR", level = eval[[1]], sheetname = eval[[2]], de.q.val = eval[[3]], 
#                    sign.val = eval[[4]], log2fold.val = eval[[5]], DEFolder = "./data/", EAFolder = "./data/")
#}

df <- calcsamplesize(sc.cells)

makevolcano(level=-1, FCutoff = .322, xpadding = .5, ypadding = 1, col = c("grey30", "purple", "cyan", "red"))
makevolcano(level=0, FCutoff = .322, xpadding = .5, ypadding = 1, col = c("grey30", "purple", "cyan", "red"))
makevolcano(level=1, FCutoff = .322, xpadding = .5, ypadding = 1, col = c("grey30", "purple", "cyan", "red"))


davidplot("DAVIDresults/gran chart defaults.txt", "DAVIDresults/gran ids.txt",
          plotname = "images/DAVID_gran.png")

davidplot("DAVIDresults/mural chart defaults.txt", "DAVIDresults/mural ids.txt",
          plotname = "images/DAVID_mural.png")



sc.cells <- readRDS("saveddata/mapped_cells.rds")
geneplot <- VlnPlot(sc.cells, c("Srsf5"), group.by = "Level0", split.by = "stim", split.plot = TRUE) + xlab("")
ggsave("images/Srsf5.png", geneplot, width = 6, height = 4)

regulationplots("data/DE_level0.xlsx", "data/DE_level1.xlsx",
                c("Granulosa", "Mesenchyme", "Epithelium", "Endothelium", "Immune") %>% rev,
                c('Immune', 'Epithelium', 'Endothelium', 'M-Steroidogenic Theca',
                  'M-Steroidogenic Stroma', 'M-Smooth Muscle', 'M-Pericyte', 'M-Fibroblast-like Stroma',
                  'M-Early Theca',
                  'GC-Mural', 'GC-Mitotic', 'GC-Luteinizing', "GC-Regressing CL", "GC-Active CL",
                  'GC-Cumulus', 'GC-Atretic'),
                20,
                "images/regulation0.png", "images/regulation1.png")

if (FALSE) {
  prepareCCC()
}
doCCC()

  