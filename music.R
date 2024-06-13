library(MuSiC2)

benchmark.eset = readRDS(url("https://xuranw.github.io/MuSiC/data/bulk-eset.rds"))
bulk.control.mtx = exprs(benchmark.eset)[, benchmark.eset$group == 'healthy']
bulk.case.mtx = exprs(benchmark.eset)[, benchmark.eset$group == 't2d']

seger.sce = readRDS(url("https://xuranw.github.io/MuSiC/data/EMATBsce_healthy.rds"))
