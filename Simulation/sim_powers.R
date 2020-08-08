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
#source("../fdapace/R/FPCA.R")
library(Rcpp, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
devtools::load_all("../SRfpca/.")
library(MASS)
source("../Funcs/urban_detect.R")
source("./sim_power_prep.R")
registerDoParallel(cores = 50)
mus <- readRDS("./sim_mu1_mu2.rds")
mu1 <- mus$mu1; diff.func <- mus$mu2 - mus$mu1; mu.grid <- mus$grid
cat("range of grid: ", range(mu.grid), "\n")
# function of doing univariate FPCA
GetFPCA <- function(Ly, Lt) {
  pca <- FPCA(Ly, Lt, optns = list(dataType = "Sparse", kernel = "epan", methodBwMu = 'GCV', methodBwCov = 'GCV', FVEthreshold = 1,
                                   maxK = 20, methodMuCovEst = "smooth", methodXi = "CE"))
  #Sig <- pca$fittedCov + diag(pca$rho, nrow = length(pca$workGrid))
  xicov <- lapply(Lt, function(tt) {
    Phi <- ConvertSupport(fromGrid = pca$workGrid, toGrid = tt, phi = pca$phi)
    Cov <- ConvertSupport(fromGrid = pca$workGrid, toGrid = as.numeric(tt), Cov = as.matrix(pca$fittedCov))
    Lam <- diag(pca$lambda, nrow = pca$selectK)
    Lam %*% t(Phi) %*% solve(Cov+diag(pca$sigma2, nrow = length(tt))) %*% Phi %*% Lam
  })
  xicov <- (1/length(Lt))*Reduce('+', xicov)
  res <- list(xi = pca$xiEst, lam = pca$lambda, xicov = xicov, cumFVE = pca$cumFVE, bw = pca$bwMu, years = 1:100)
  return(res)
}

#######################################################################################
cat("Start doing the simulation..", "\n")
cs <- c(0,0.025,0.05,0.075,0.1,0.15,0.2,0.25); sigma <- 1.2
#cs <- c(0.3, 0.35, 0.4, 0.45, 0.5)
# simulate 500 sites with change-point, at scale c
for (c in cs) {
cat("Doing c = ", c, "\n")
simres0 <- foreach(i = 1:500) %dopar% {
  #cat("Doing the ", i, "th simulation...", "\n")
  sim_dati <- fChangeSim(mu1, mu.grid, diff.func, c = c, n = 100, K = 20, theta = 0.5, e.sd = sigma, decay = "slow")
  #saveRDS(sim_dati, file = paste0("debug_sim_dat_", i, ".rds"))
  Lyi <- sim_dati$Ly; Lti <- sim_dati$Lt
  Ly1i <- lapply(Lyi, function(y) {y[1,]})
  Ly2i <- lapply(Lyi, function(y) {y[2,]})
  Ly3i <- lapply(Lyi, function(y) {y[3,]})
  # do MFPCA
  mpca <- RFPCA(Lyi, Lti, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean",
                               kernel = "epan", maxK = 20, methodMuCovEst = "smooth", methodXi = "CE"))
  xicov <- lapply(Lti, function(tt) {
    m <- length(tt)
    Phi <- apply(mpca$phiObsTrunc[match(tt, mpca$obsGrid),,,drop=FALSE], 3, function(x) {c(t(x))})
    Lam <- diag(mpca$lam, nrow = mpca$K)
    Cov0 <- mpca$covObs[match(tt, mpca$obsGrid), match(tt, mpca$obsGrid),,,drop=FALSE]
    Cov <- matrix(NA, nrow = m*3, ncol = m*3)
    for (i in 1:m) {
      indi <- (3*(i-1)+1):(3*i)
        for (j in 1:m) {
          indj <- (3*(j-1)+1):(3*j)
          Cov[indi, indj] <- Cov0[i,j,,] 
      }
    }
    Lam %*% t(Phi) %*% solve(Cov+diag(mpca$sigma2, nrow = m*3)) %*% Phi %*% Lam
  })
  xicov <- (1/length(Lti))*Reduce('+', xicov)
  xi.res <- list(xi = mpca$xi, lam = mpca$lam, xicov = xicov, cumFVE = mpca$cumFVE, bw = mpca$userBwMu, years = 1:100)
  # do FPCA
  #pca1 <- FPCA(Ly1i, Lti, optns = list(dataType = "Sparse", kernel = "epan", methodBwMu = 'GCV', FVEthreshold = 1,
  #                                    maxK = 20, methodMuCovEst = "smooth", methodXi = "CE"))
  xi.res1 <- GetFPCA(Ly1i, Lti)
  #pca2 <- FPCA(Ly2i, Lti, optns = list(dataType = "Sparse", kernel = "epan", methodBwMu = 'GCV', FVEthreshold = 1,
  #                                     maxK = 20, methodMuCovEst = "smooth", methodXi = "CE"))
  xi.res2 <- GetFPCA(Ly2i, Lti)
  #pca3 <- FPCA(Ly3i, Lti, optns = list(dataType = "Sparse", kernel = "epan", methodBwMu = 'GCV', FVEthreshold = 1,
  #                                     maxK = 20, methodMuCovEst = "smooth", methodXi = "CE"))
  xi.res3 <- GetFPCA(Ly3i, Lti)
  #xi.res1 <- list(xi = pca1$xiEst, lam = pca1$lambda, cumFVE = pca1$cumFVE, years = 1:100)
  #xi.res2 <- list(xi = pca2$xiEst, lam = pca2$lambda, cumFVE = pca2$cumFVE, years = 1:100)
  #xi.res3 <- list(xi = pca3$xiEst, lam = pca3$lambda, cumFVE = pca3$cumFVE, years = 1:100)
  
  list(MFPCA = xi.res, UFPCA1 = xi.res1, UFPCA2 = xi.res2, UFPCA3 = xi.res3, data = sim_dati)
  #list(MFPCA = mpca, UFPCA1 = pca1, UFPCA2 = pca2, UFPCA3 = pca3, data = sim_dati)
  #saveRDS(res, file = paste0("./debug_data/res/sim_res_", i, ".rds"))
}
simres <- list()
simres$MFPCA <- lapply(simres0, function(x) {x$MFPCA})
simres$UFPCA1 <- lapply(simres0, function(x) {x$UFPCA1})
simres$UFPCA2 <- lapply(simres0, function(x) {x$UFPCA2})
simres$UFPCA3 <- lapply(simres0, function(x) {x$UFPCA3})
simdata <- lapply(simres0, function(x) {x$data})
saveRDS(simres, file = paste0("./slow/sigma_", sigma, "/sim_500_MFPCA_UFPCA_slow_", c, ".rds"))
saveRDS(simdata, file = paste0("./slow/sigma_", sigma, "/sim_500_data_slow_", c, ".rds"))
}

