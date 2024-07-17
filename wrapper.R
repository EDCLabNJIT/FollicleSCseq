suppressWarnings(devtools::document("./scRNAhelper/"))
sc.cells <- mergedata(min.features = 1000, min.cells = 3,
                      TenX.folder = "./10X/",
                      saverds = TRUE, write.rdsfile = "./saveddata/merged_cells.rds")

sc.cells <- QC(NULL, make.plots = TRUE, plotLocation = "./images/",
               nFeature.lower = 3500, nFeature.upper = Inf, percent.mt.upper = 5, nCount.rna.upper = 100000, findDbls = TRUE,
               saverds = TRUE, write.rdsfile = "./saveddata/filtered_cells.rds")

sc.cells <- preprocess(make.plots = FALSE, plotLocation = "./images/",
                       nfeatures = 5000,
                       saverds = TRUE, write.rdsfile = "./saveddata/preprocessed_cells.rds")

sc.cells <- cluster(seuratObj = NULL, seuratFile = "./saveddata/preprocessed_cells.rds", PC = 15,
                    resolution = 0.01, identities = c("0" = "Granulosa","1" = "Mesenchyme","2" = "Endothelium","3" = "Immune", "4" = "Epithelium"),
                    makeplots = TRUE, plotLocation = "./images/", dotmarkers = "./morrismarkers/morris_markers.csv",
                    saverds = TRUE, write.rdsfile = "./saveddata/clustered_cells.rds")

#                  old name                  new name
atlas.rename <- c("GC_Estrous"            = "GC_Luteinizing and CL",
                  "GC_CL_Lytic"           = "GC_Luteinizing and CL",
                  "GC_CL_Active"          = "GC_Luteinizing and CL",
                  "GC_Antral"             = "GC_Mural",
                  "GC_Preantral"          = "GC_Cumulus",
                  "M_Immature Theca"      = "M_Early Theca",
                  "M_Dividing Mesenchyme" = "M_Early Theca", 
                  "M_Cortical Stroma"     = "M_Steroidogenic Stroma",
                  "M_Medullary Stroma"    = "M_Fibroblast-like Stroma")

sc.cells <- mapcells(seuratObj = NULL, seuratFile = "./saveddata/clustered_cells.rds",
                     atlasFile = "./morrisdata/ovary_0.rds",
                     make.plots = TRUE, plotLocation = "./images/",
                     atlas.rename = atlas.rename, PC = 20, umapPC = 20,
                     anchor.n.trees = 50, anchor.k.anchor = 5, anchor.k.score = 30,
                     saverds = TRUE, write.rdsfile = "./saveddata/mapped_cells.rds",
                     marker.folder = "./morrismarkers/",
                     projumapargs = list(annoy.metric = "cosine"))

DE.data <- calcDE(NULL, 
                  savedata = TRUE, DEFolder = "./data/")

#                   level sheetname  p/qval reg log2fold cutoff
enrichvals = list(list(1, "GC-Mural" , .05, 0 , 0),
                  list(1, "GC-Mural" , .05, 1 , 0),
                  list(1, "GC-Mural" , .05, -1, 0),
                  list(1, "GC-Mural" , .05, 0 , .3),
                  list(1, "GC-Mural" , .05, 1 , .3),
                  list(1, "GC-Mural" , .05, -1, .3),
                  list(0, "Granulosa", .05, 0 , 0),
                  list(0, "Granulosa", .05, 1 , 0),
                  list(0, "Granulosa", .05, -1 , 0))
for (eval in enrichvals) {
  EA.data <- enrich(p.val = "FDR", level = eval[[1]], sheetname = eval[[2]], de.q.val = eval[[3]], 
                    sign.val = eval[[4]], log2fold.val = eval[[5]], DEFolder = "./data/", EAFolder = "./data/")
}

t <- calcsamplesize(sc.cells)

makevolcano(level=0, FCutoff = .3)
makevolcano(level=1, FCutoff = .3)


