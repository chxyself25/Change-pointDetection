library(doParallel)
library(devtools)
library(fda, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(fdapace, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(Rcpp, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir("../fdapace/R/", trace = FALSE)
devtools::load_all("../SRfpca/.")
library(MASS)
source("../Funcs/urban_detect_correct.R")
registerDoParallel(cores = 50)
# centered functional data in 373 locations with change observed
data <- readRDS("../DataResults/ForTesting/Landsat_All_3Indices_5Bands_Filtered_Sparsesd_centered_MMCE3_fve55.rds")
ids <- unique(data$pointID);
#BSdist1 <- readRDS("../DataResults/ForTesting/bootstrap1000_All_mv_3Indices_1.rds")
#len1 <- length(BSdist1$S)

cat("Start doing the bootstrap sampling...", "\n")

BSdist <- list(S = list(), A = list(), lams = list())
#BSdist <- readRDS("../DataResults/ForTesting/bootstrap1000_All_mv_3Indices.rds")
for (id in ids) {
  #if (id <= len1 && !is.null(BSdist1$S[[id]])) {
  #  next
  #}
  cat("doing location ", id, '\n')
  dati <- subset(data, pointID == id)
  Ly <- lapply(split(dati, f = as.factor(dati$year)), function(x) {
    t(as.matrix(x[, c("fsd_ndvi", "fsd_mndwi", "fsd_b7")]))
  })
  Lt <- lapply(split(dati, f = as.factor(dati$year)), function(x) {x$doy})
  ## bootstrap sampling to get S and A statistics
  btsi <- mfBTSampling(Ly, Lt, K = 24, bs.n = 100, resample.n = 1000)
  BSdist$S[[id]] <- btsi$S
  BSdist$A[[id]] <- btsi$A
  BSdist$lams[[id]] <- btsi$lams
  saveRDS(BSdist, file = "../DataResults/ForTesting/bootstrap1000_All_mv_3Indices_lams.rds")
}
saveRDS(BSdist, file = "../DataResults/ForTesting/bootstrap1000_All_mv_3Indices_lams.rds")
