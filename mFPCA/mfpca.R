sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
library(devtools)
library(fda, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(fdapace, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(Rcpp, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
sourceDir("../fdapace/R/", trace = FALSE)
devtools::load_all("../SRfpca/.")
library(MASS)
library(doParallel)
registerDoParallel(cores = 50)
data <- readRDS("../DataResults/CleanData/Landsat_All_3Indices_5Bands_Filtered_Sparsesd.rds")
ids <- unique(data$pointID)
cat("3 indices...", "\n")
pcas <- foreach(id = 1:503) %dopar% {
  if (!id %in% ids) {
    return(NULL)
  }
  dati <- subset(data, pointID == id)
  Ly <- lapply(split(dati, f = as.factor(dati$year)), function(x) {
      t(as.matrix(x[,c("fsd_b7", "fsd_ndvi", "fsd_mndwi")]))})
  Lt <- lapply(split(dati, f = as.factor(dati$year)), function(x) {x$doy})
  mpca <- RFPCA(Ly, Lt, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean",
                kernel = "epan", maxK = 24, methodMuCovEst = "smooth", methodXi = "CE"))
  #resi <- GetFPCA(datai, bands = c("NDVI", "MNDWI", "B7"))
  xicov <- lapply(Lt, function(tt) {
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
  xicov <- (1/length(Lt))*Reduce('+', xicov)
  resi <- list(xi = mpca$xi, lam = mpca$lam, xicov = xicov, cumFVE = mpca$cumFVE, bw = mpca$userBwMu, years = unique(dati$year))
  return(resi)
}
#ndvi <- lapply(pcas, function(x) {x[["NDVI"]]})
#mndwi <- lapply(pcas, function(x) {x[["MNDWI"]]})
#b7 <- lapply(pcas, function(x) {x[["B7"]]})
#saveRDS(ndvi, file = "../DataResults/FPCAResults/xilam_est_All_uFPCA_ndvi_GetFPCA.rds")
#saveRDS(mndwi, file = "../DataResults/FPCAResults/xilam_est_All_uFPCA_mndwi_GetFPCA.rds")
#saveRDS(ndvi, file = "../DataResults/FPCAResults/xilam_est_All_uFPCA_b7_GetFPCA.rds")
saveRDS(pcas, file = "../DataResults/FPCAResults/xilam_est_All_mFPCA_3Indices.rds")
cat("5 bands...", "\n")
pcas <- foreach(id = 1:503) %dopar% {
  if (!id %in% ids) {
    return(NULL)
  }
  dati <- subset(data, pointID == id)
  Ly <- lapply(split(dati, f = as.factor(dati$year)), function(x) {
      t(as.matrix(x[,c("fsd_b2", "fsd_b3", "fsd_b4", "fsd_b5", "fsd_b7")]))})
  Lt <- lapply(split(dati, f = as.factor(dati$year)), function(x) {x$doy})
  mpca <- RFPCA(Ly, Lt, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean",
                kernel = "epan", maxK = 24, methodMuCovEst = "smooth", methodXi = "CE"))
  #resi <- GetFPCA(datai, bands = c("NDVI", "MNDWI", "B7"))
  xicov <- lapply(Lt, function(tt) {
    m <- length(tt)
    Phi <- apply(mpca$phiObsTrunc[match(tt, mpca$obsGrid),,,drop=FALSE], 3, function(x) {c(t(x))})
    Lam <- diag(mpca$lam, nrow = mpca$K)
    Cov0 <- mpca$covObs[match(tt, mpca$obsGrid), match(tt, mpca$obsGrid),,,drop=FALSE]
    Cov <- matrix(NA, nrow = m*5, ncol = m*5)
    for (i in 1:m) {
      indi <- (5*(i-1)+1):(5*i)
        for (j in 1:m) {
          indj <- (5*(j-1)+1):(5*j)
          Cov[indi, indj] <- Cov0[i,j,,]
      }
    }
    Lam %*% t(Phi) %*% solve(Cov+diag(mpca$sigma2, nrow = m*5)) %*% Phi %*% Lam
  })
  xicov <- (1/length(Lt))*Reduce('+', xicov)
  resi <- list(xi = mpca$xi, lam = mpca$lam, xicov = xicov, cumFVE = mpca$cumFVE, bw = mpca$userBwMu, years = unique(dati$year))
  return(resi)
}
saveRDS(pcas, file = "../DataResults/FPCAResults/xilam_est_All_mFPCA_5Bands.rds")
