library(doParallel)
library(devtools)
library(fda, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(fdapace, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir("../fdapace/R/", trace = FALSE)
devtools::load_all("../SRfpca/.")
library(Rcpp, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(MASS)
registerDoParallel(cores = 50)
SNR <- function(diff.func, grid, lams, sigma2, theta) {
  delta <- sum(apply(diff.func, 1, function(y) {trapzRcpp(grid, y^2)}))
  p <- nrow(diff.func)
  num <- theta*(1-theta)*delta
  den <- sum(lams) + sigma2*diff(range(grid))*p
  return(num/den)
}


## read in data
dat <- readRDS("../DataResults/CleanData/Landsat_All_3Indices_5Bands_Filtered_Sparsesd.rds")
## reference data
ref1 <- read.csv(file = "../DataResults/LandsatRaw_ForDrZhu_merged_reference.csv")
ref1$pointID <- ref1$pointID + 120
ref2 <- read.csv(file = "../DataResults/LandsatRaw_ForDrZhu_merged0_reference.csv")
ref2$pointID <- ref2$pointID + 263
# remove ambiguous cases and unchange sites
rm.idx <- c(16,17,94,143) + 120
rm.idx <- c(rm.idx, 34+263)
# test if removing these id could improve results
ref <- subset(rbind(ref1, ref2), !is.na(year) & !(pointID %in% rm.idx) & year > 2000)
ids <- ref$pointID
cat("check number of locations: ", length(ids), "\n")
grid.len <- 200

## calculate SNR for each location
res <- foreach(id = ids, .combine = "rbind") %dopar% {
  dati <- subset(dat, pointID == id)
  # change-point
  P2 <- ref$year[ref$pointID == id]
  # estimate mean before and after the change
  yeari <- unique(dati$year)
  Ly <- lapply(split(dati, f = as.factor(dati$year)), function(x) {
    t(as.matrix(x[,c("fsd_b2", "fsd_b3", "fsd_b4", "fsd_b5", "fsd_b7", "fsd_ndvi", "fsd_mndwi")]))})
  Lt <- lapply(split(dati, f = as.factor(dati$year)), function(x) {x$doy})
  Ly1 <- Ly[yeari<=P2]; Lt1 <- Lt[yeari<=P2]
  Ly2 <- Ly[yeari>P2]; Lt2 <- Lt[yeari>P2]
  day.min <- max(min(unlist(Lt1)), min(unlist(Lt2))); day.max <- min(max(unlist(Lt1)), max(unlist(Lt2)))
  mui1 <- RFPCAMu(Ly1, Lt1, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean", nRegGrid = grid.len,  
                                         kernel = "epan", maxK = 24, methodMuCovEst = "smooth", ToutRange = c(day.min, day.max)))
  mui2 <- RFPCAMu(Ly2, Lt2, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean", nRegGrid = grid.len,
                                         kernel = "epan", maxK = 24, methodMuCovEst = "smooth", ToutRange = c(day.min, day.max)))
  mu.grid <- seq(day.min, day.max, length.out = grid.len)
  cat(range(mui1$workGrid - mui2$workGrid), "\n")
  cat(range(mui1$workGrid - mu.grid), "\n")
  # do FPCA to get lambda and sigma2
  resid1 <- lapply(1:length(Ly1), function(i) {
    Ly1[[i]] - mui1$muObs[, match(Lt1[[i]], mui1$obsGrid)]
  })
  resid2 <- lapply(1:length(Ly2), function(i) {
    Ly2[[i]] - mui2$muObs[, match(Lt2[[i]], mui2$obsGrid)]
  })
  resid <- append(resid1, resid2)
  resid3 <- lapply(resid, function(x) {x[c("fsd_ndvi", "fsd_mndwi", "fsd_b7"), ]})
  pcai3 <- RFPCA(resid3, Lt, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean", nRegGrid = grid.len, 
                                kernel = "epan", maxK = 24, methodMuCovEst = "smooth", methodXi = "CE", ToutRange = c(day.min, day.max)))
  snr3 <-  SNR((mui2$muWork-mui1$muWork)[c("fsd_ndvi", "fsd_mndwi", "fsd_b7"), ], mu.grid, pcai3$lam, pcai3$sigma2, mean(yeari <= P2))
  resid5 <- lapply(resid, function(x) {x[c("fsd_b2", "fsd_b3", "fsd_b4", "fsd_b5", "fsd_b7"), ]})
  pcai5 <- RFPCA(resid5, Lt, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean", nRegGrid = grid.len, 
                                kernel = "epan", maxK = 24, methodMuCovEst = "smooth", methodXi = "CE", ToutRange = c(day.min, day.max)))
  snr5 <- SNR((mui2$muWork-mui1$muWork)[c("fsd_b2", "fsd_b3", "fsd_b4", "fsd_b5", "fsd_b7"), ], mu.grid, pcai5$lam, pcai5$sigma2, mean(yeari <= P2))
  data.frame(pointID = id, SNR = c(snr3, snr5), Dimension = c("3 indicators", "5 bands"), sigma2 = c(pcai3$sigma2, pcai5$sigma2))
}

saveRDS(res, file = "./SNR_all_3_5.rds")
