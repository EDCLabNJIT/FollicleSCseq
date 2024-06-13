devtools::document("./scRNAhelper/")
sc.cells <- mergedata(min.features = 1000, min.cells = 3,
                      TenX.folder = "./10X/",
                      saverds = TRUE, write.rdsfile = "./saveddata/merged_cells.rds")

sc.cells <- QC(sc.cells, make.plots = TRUE, plotLocation = "./images/",
               nFeature.lower = 3500, nFeature.upper = 7000, percent.mt.upper = 5, nCount.rna.upper = 50000,
               saverds = TRUE, write.rdsfile = "./saveddata/filtered_cells.rds")

sc.cells <- preprocess(make.plots = FALSE, plotLocation = "./images/",
                       nfeatures = 500,
                       saverds = TRUE, write.rdsfile = "./saveddata/preprocessed_cells.rds")

#                  old name                  new name
atlas.rename <- c("GC_Estrous"            = "GC_Luteinizing",
                  "GC_CL_Lytic"           = "GC_Regressing CL",
                  "GC_CL_Active"          = "GC_Active CL",
                  "M_Immature Theca"      = "M_Early Theca",
                  "M_Dividing Mesenchyme" = "M_Early Theca", 
                  "M_Cortical Stroma"     = "M_Steroidogenic Stroma",
                  "M_Medullary Stroma"    = "M_Fibroblast-like Stroma")

sc.cells <- mapcells(sc.cells, atlasFile = "./morrisdata/ovary_0.rds",
                make.plots = TRUE, plotLocation = "./images/",
                atlas.rename = atlas.rename, PC = 15, umapPC = 20,
                anchor.n.trees = 50, anchor.k.anchor = 5, anchor.k.score = 50,
                saverds = TRUE, write.rdsfile = "./saveddata/mapped_cells.rds",
                feature.folder = NULL,
                projumapargs = list(annoy.metric = "cosine"))

DE.data <- calcDE(sc.cells,
                  savedata = TRUE, DEFolder = "./data/")

enrichvals = list(c(1, "GC-Antral"),
                  c(0, "Granulosa"))
for (eval in enrichvals) {
  EA.data <- enrich(level = eval[1], sheetname = eval[2])
}

makedot(sc.cells,
        plotLocation = "./images/",
        morris.gran.markers = "./morrismarkers/morris_gran_markers.csv",
        morris.mes.merkers = "./morrismarkers/morris_mes_markers.csv",
        nGenes = 5)

t <- calcsamplesize(sc.cells)

makevolcano(level=0)
makevolcano(level=1)
