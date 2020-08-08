################################################################################
########### simulation performance w.r.t different measurement error ###########
################################################################################
source("../Funcs/urban_detect_correct.R")
source("../Funcs/eval_func.R")
library(doParallel)
registerDoParallel(cores = 50)
Kdist <- readRDS("../Kd_simulation/Kd_ecdf_25.rds")
FVEs <- c(55:99, 99.5, 99.9)
# cs <- c(0,0.05,0.1,0.2,0.4)
cs <- c(0,0.025,0.05,0.075, 0.1, 0.15, 0.2, 0.25)
type = "fast"; sigmas <- c(0.3, 0.6, 1.2)

if (FALSE) {
for (sigma in sigmas) {
dir <- paste0("./", type, "/sigma_", sigma,"/")
res <- NULL
for (c in cs) {
  cat("c = ", c, "\n")
  #BSdist.list <- readRDS(paste0("./SimulationData/fast_simulation/sim_500_BSdist3_fast_", c, ".rds"))
  simres <- readRDS(paste0(dir, "/sim_500_MFPCA_UFPCA_", type, "_", c, ".rds"))
  for (f in FVEs) {
    # MMCE
    resf <- ftest_multi(simres$MFPCA, ids = 1:500, year.range = NULL, FVE = f, beta = NULL, K = NULL,
                        method = "max", type = "mv", BSres = Kdist, Kd = TRUE)
    res <- rbind(res, data.frame(pval = resf$det, method = "mMAX", fve = f, c = c))
    # MSCE
    resf <- ftest_multi(simres$MFPCA, ids = 1:500, year.range = NULL, FVE = f, beta = NULL, K = NULL,
                        method = "sum", type = "mv", BSres = Kdist, Kd = TRUE)
    res <- rbind(res, data.frame(pval = resf$det, method = "mSUM", fve = f, c = c))
    # UMCE
    resf <- ftest_multi(list(simres$UFPCA1, simres$UFPCA2, simres$UFPCA3), ids = 1:500, year.range = NULL, 
                        FVE = f, beta = NULL, K = NULL, method = "max", type = "uv", BSres = Kdist, Kd = TRUE)
    res <- rbind(res, data.frame(pval = resf$det, method = "uMAX", fve = f, c = c))
    # USCE
    # resf <- ftest_multi(list(simres$UFPCA1, simres$UFPCA2, simres$UFPCA3), ids = 1:500, year.range = NULL, 
    #                     FVE = f, beta = NULL, K = NULL, method = "sum", type = "uv", BSres = BSdist.list, Kd = FALSE)
    # res <- rbind(res, data.frame(rej = mean(resf$det), method = "USCE", fve = f, c = c))
  }
}
saveRDS(res, file = paste0(dir, "/sim_500_", type, "_", sigma, "_pvals_res.rds"))
}
}

## use bootstrap distribution to do detetcion
if (TRUE) {
for (sigma in sigmas) {
dir <- paste0("./", type, "/sigma_", sigma,"/")
res <- NULL
for (c in cs) {
  #cat("c = ", c, "\n")
  if (file.exists(paste0(dir, "/sim_500_mv_bootstrap_", type, "_", c, ".rds"))) {
  cat("c = ", c, "\n")
  simres <- readRDS(paste0(dir, "/sim_500_MFPCA_UFPCA_", type, "_", c, ".rds"))
  mbts <- readRDS(paste0(dir, "/sim_500_mv_bootstrap_", type, "_", c, ".rds"))
  ubts <- readRDS(paste0(dir, "/sim_500_uv_bootstrap_", type, "_", c, ".rds"))
  #min(sapply(ubts, function(x) {sapply(x$S, function(y) {nrow(y)})}))
  for (f in FVEs) {
    # MMCE
    resf <- ftest_multi(simres$MFPCA, ids = 1:500, year.range = NULL, FVE = f, beta = NULL, K = NULL,
                        method = "max", type = "mv", BSres = mbts, Kd = FALSE)
    res <- rbind(res, data.frame(pval = resf$det, method = "mMAX", fve = f, c = c))
    # MSCE
    resf <- ftest_multi(simres$MFPCA, ids = 1:500, year.range = NULL, FVE = f, beta = NULL, K = NULL,
                        method = "sum", type = "mv", BSres = mbts, Kd = FALSE)
    res <- rbind(res, data.frame(pval = resf$det, method = "mSUM", fve = f, c = c))
    # UMCE
    resf <- ftest_multi(list(simres$UFPCA1, simres$UFPCA2, simres$UFPCA3), ids = 1:500, year.range = NULL,
                        FVE = f, beta = NULL, K = NULL, method = "max", type = "uv", BSres = ubts, Kd = FALSE)
    res <- rbind(res, data.frame(pval = resf$det, method = "uMAX", fve = f, c = c))
    # USCE
    resf <- ftest_multi(list(simres$UFPCA1, simres$UFPCA2, simres$UFPCA3), ids = 1:500, year.range = NULL, 
                        FVE = f, beta = NULL, K = NULL, method = "sum", type = "uv", BSres = ubts, Kd = FALSE)
    res <- rbind(res, data.frame(pval = resf$det, method = "uSUM", fve = f, c = c))
  }
 }
}
saveRDS(res, file = paste0(dir, "/sim_500_", type, "_", sigma, "_pvals_bts_res.rds"))
}
}

if (FALSE) {
for (sigma in sigmas) {
dir <- paste0("./", type, "/sigma_", sigma,"/")
res <- NULL
for (c in cs) {
  cat("c = ", c, "\n")
  # BSdist.list <- readRDS(paste0("./Testing/simulation/fast_simulation/sim_500_BSdist3_fast_", c, ".rds"))
  simres <- readRDS(paste0(dir, "/sim_500_MFPCA_UFPCA_", type, "_", c, ".rds"))
  for (f in FVEs) {
    # MMCE
    resf <- fchange_multi(simres$MFPCA, ids = 1:500, year.range = NULL, FVE = f, beta = NULL, K = NULL,
                        method = "max", type = "mv")
    res <- rbind(res, data.frame(pointID = resf$pointID, estimate = resf$det, method = "mMAX", fve = f, c = c))
    # MSCE
    resf <- fchange_multi(simres$MFPCA, ids = 1:500, year.range = NULL, FVE = f, beta = NULL, K = NULL,
                        method = "sum", type = "mv")
    res <- rbind(res, data.frame(pointID = resf$pointID, estimate = resf$det, method = "mSUM", fve = f, c = c))
    # UMCE
    resf <- fchange_multi(list(simres$UFPCA1, simres$UFPCA2, simres$UFPCA3), ids = 1:500, year.range = NULL, 
                        FVE = f, beta = NULL, K = NULL, method = "max", type = "uv")
    res <- rbind(res, data.frame(pointID = resf$pointID, estimate = resf$det, method = "uMAX", fve = f, c = c))
    # USCE
    resf <- fchange_multi(list(simres$UFPCA1, simres$UFPCA2, simres$UFPCA3), ids = 1:500, year.range = NULL,
                        FVE = f, beta = NULL, K = NULL, method = "sum", type = "uv")
    res <- rbind(res, data.frame(pointID = resf$pointID, estimate = resf$det, method = "uSUM", fve = f, c = c))
  }
}
saveRDS(res, file = paste0(dir, "/sim_500_", type, "_", sigma, "_est_res.rds"))
}
}

