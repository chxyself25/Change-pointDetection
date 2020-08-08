library(doParallel)
library(devtools)
library(fda, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(Rcpp, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
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
print(sessionInfo())
source("../Funcs/urban_detect_correct.R")
registerDoParallel(cores = 50)
fres <- data.frame(type = rep(c("slow", "fast"), each = 3), sigma = rep(c(0.3,0.6,1.2), 2),
                   fve = c(99.5, 99.9, 99.9, 99, 99.9, 99.9))

## information set up
type <- "slow"; sigma <- 1.2
dir <- paste0("./", type, "/sigma_", sigma, "/")
est.res <- readRDS(paste0(dir, "sim_500_", type, "_", sigma, "_est_res.rds"))
## best FVE for mSUM 
opt.fve <- fres$fve[fres$type == type & fres$sigma == sigma]
est <- subset(est.res, method == "mSUM" & fve == opt.fve)
cs <- c(0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25)

## get centered functional data
yeari <- 1:100; c = cs[8]
dat <- readRDS(paste0(dir, "sim_500_data_", type, "_", c, ".rds"))
BSdist <- list(list(S = list(), A = list()), list(S = list(), A = list()), list(S = list(), A = list()))
mBSdist <- list(S = list(), A = list())
#BSdist <- readRDS(paste0(dir, "sim_500_uv_bootstrap_", type, "_",  c, ".rds"))
#mBSdist <- readRDS(paste0(dir, "sim_500_mv_bootstrap_", type, "_", c, ".rds"))
for (id in 1:500) {
  cat("doing the id ", id, "\n")
  dati <- dat[[id]]
  Ly <- dati$Ly; Lt <- dati$Lt
  P <- est$estimate[est$c == c & est$pointID == id]
  Ly1 <- Ly[yeari<=P]; Lt1 <- Lt[yeari<=P]
  Ly2 <- Ly[yeari>P]; Lt2 <- Lt[yeari>P]
  mui1 <- RFPCAMu(Ly1, Lt1, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean", 
                                         kernel = "epan", maxK = 20, methodMuCovEst = "smooth"))
  mui2 <- RFPCAMu(Ly2, Lt2, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean", 
                                         kernel = "epan", maxK = 20, methodMuCovEst = "smooth"))
  resid1 <- lapply(1:length(Ly1), function(i) {
    Ly1[[i]] - mui1$muObs[, match(Lt1[[i]], mui1$obsGrid)]
  })
  resid2 <- lapply(1:length(Ly2), function(i) {
    Ly2[[i]] - mui2$muObs[, match(Lt2[[i]], mui2$obsGrid)]
  })
  resid <- append(resid1, resid2)
  for (b in 1:3) {
    cat("doing band ", b, "\n")
    Lyb <- lapply(resid, function(x) {x[b,]})
    #b.idx <- sample(1:100, size = 100, replace = TRUE)
    #pcab <- FPCA(Ly[b.idx], Lt[b.idx], optns = list(dataType = "Sparse", kernel = "epan", methodBwMu = "GCV", methodBwCov = "GCV",
    #                                  maxK = 20, methodMuCovEst = "smooth", methodXi = "CE", FVEthreshold = 1))
    #cat(length(pcab), "\n")
    btsib <- fBTSampling(Lyb, Lt, K = 20, bs.n = 100, resample.n = 500)
    BSdist[[b]]$S[[id]] <- btsib$S
    BSdist[[b]]$A[[id]] <- btsib$A
  }
  btsi <- mfBTSampling(resid, Lt, K = 20, bs.n = 100, resample.n =  500)
  mBSdist$S[[id]] <- btsi$S
  mBSdist$A[[id]] <- btsi$A
  cat(id, " is finished", "\n")
  saveRDS(BSdist, file = paste0(dir, "sim_500_uv_bootstrap_", type, "_",  c, ".rds"))
  saveRDS(mBSdist, file = paste0(dir, "sim_500_mv_bootstrap_", type, "_", c, ".rds"))
}
# backup in lss folder
#saveRDS(BSdist, file = paste0("/lss/research/zhuz-lab/xchang/Change-pointDetection/Simulation/", type, "/sigma_", sigma, "/sim_500_uv_bootstrap_", type, "_",  c, ".rds"))
#saveRDS(mBSdist, file = paste0("/lss/research/zhuz-lab/xchang/Change-pointDetection/Simulation/", type, "/sigma_", sigma, "/sim_500_mv_bootstrap_", type, "_",  c, ".rds"))


