sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
library(doParallel)
library(signal, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(dplyr, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(fdapace, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(Rcpp, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
sourceDir("./Funcs/", trace = FALSE)
sourceDir("../fdapace/R/", trace = FALSE)
registerDoParallel(cores = 25)

Kd_dist <- readRDS("../Kd_simulation/Kd_ecdf_25.rds")
names.used = c('pID', 'doy', 'year', 'B1', 'B2', 'B3', 'B4', 'B5', 'B7')
#dsm.test <- readRDS("./DSM_test.rds")
#res.test <- UfchangepointAll(dsm.test, var.names = names.used, bands = c('ndvi', 'mndwi', 'b7'), Kdist =  Kd_dist, detect.method = 'MIE', est.method = 'SCE',
#                 year.range = NULL, FVE = c(99, 85), K = c(NULL, NULL), alpha = 0.05, 
#                 save.results = TRUE, use.results = FALSE, dir = "./fchange_test_results/")
#saveRDS(res.test, file = "./fchange_test_results/test_fchange_res.rds")
for (seg in 1:14) {
cat(seg, "segment of Ankeny dataset", "\n")
file.name <- paste0("./Ankeny", seg, ".rds")
dsm.seg <- readRDS(file.name)
ids <- sort(unique(dsm.seg$pID))
div <- (length(ids)-1)%/%100 + 1
cat("number of subsets: ", div, "\n")
seg.dir <- paste0("./fchange_ankeny_", seg, "/")
if (!dir.exists(seg.dir)) {
    dir.create(seg.dir)
}
res <- NULL
for (ii in 1:div) {
  cat("Doing subset data ", ii, "\n")
  dsm.segii <- subset(dsm.seg, pID %in% ids[((ii-1)*100+1):(ii*100)])
  #dsm.segii <- subset(dsm.seg, name %in% ids[1:10])
  res.seg <- UfchangepointAll(dsm.segii, var.names = names.used, bands = c('ndvi', 'mndwi', 'b7'), Kdist =  Kd_dist, detect.method = 'IME', est.method = 'CSE',
                 year.range = NULL, FVE = c(99.9, 85), K = c(NULL, NULL), alpha = 0.05, 
                 save.results = TRUE, use.results = FALSE, dir = seg.dir) 
  saveRDS(res.seg, file = paste0(seg.dir, "dsm_", seg, "_subset_", ii, "_fchange_res.rds"))
  res <- rbind(res, res.seg)
}
saveRDS(res, file = paste0(seg.dir, "ankeny_", seg, "_all_fchange_res.rds"))
}
