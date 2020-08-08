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
names.used = c('name', 'doy', 'year', 'B1', 'B2', 'B3', 'B4', 'B5', 'B7')
#dsm.test <- readRDS("./DSM_test.rds")
#res.test <- UfchangepointAll(dsm.test, var.names = names.used, bands = c('ndvi', 'mndwi', 'b7'), Kdist =  Kd_dist, detect.method = 'MIE', est.method = 'SCE',
#                 year.range = NULL, FVE = c(99, 85), K = c(NULL, NULL), alpha = 0.05, 
#                 save.results = TRUE, use.results = FALSE, dir = "./fchange_test_results/")
#saveRDS(res.test, file = "./fchange_test_results/test_fchange_res.rds")
validate <- read.csv("./dsm_validation_set.csv")
ids0 <- validate$pointID
for (seg in 1:10) {
cat(seg, "segment of Des Moines dataset", "\n")
file.name <- paste0("./DSM", seg, ".rds")
dsm.seg <- readRDS(file.name)
ids <- sort(unique(dsm.seg$name))
seg.dir <- paste0("./fchange_dsm_", seg, "/")
ids <-intersect(ids, ids0)
if (length(ids) == 0) {next}
dsm.segii <- subset(dsm.seg, name %in% ids)
cat("doing for pointID in ", ids, "\n")

res <- NULL
for (m in c('IME', 'ISE', 'CSE')) {
  for (f in c(55:99, 99.9)) {
    res.seg <- UfchangepointAll(dsm.segii, var.names = names.used, bands = c('ndvi', 'mndwi', 'b7'), Kdist =  Kd_dist, 
                            detect.method = NULL, est.method = m, year.range = NULL, FVE = c(f, f), K = c(NULL, NULL), alpha = 0.05, 
                            save.results = TRUE, use.results =TRUE, dir = seg.dir)
    res.seg$method = m; res.seg$fve = f 
    res <- rbind(res, res.seg)
  }
}
saveRDS(res, file = paste0(seg.dir, "dsm_", seg, "_validate_fchange_res.rds"))
}

