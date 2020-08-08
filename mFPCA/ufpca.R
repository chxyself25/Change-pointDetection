sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
library(fdapace, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(Rcpp, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
sourceDir("../fdapace/R/", trace = FALSE)
library(doParallel)
registerDoParallel(cores = 50)
GetFPCA <- function(dat, bands) {
  res <- list()
  years <- unique(dat$year)
  for (b in bands) {
    Ly <- split(dat[[b]], f = as.factor(dat$year))
    Lt <- split(dat$doy, f = as.factor(dat$year))
    pca <- FPCA(Ly, Lt, optns = list(dataType = 'Sparse', kernel = 'epan', methodBwMu = 'GCV', methodBwCov = 'GCV',
                              methodMuCovEst = 'smooth', methodXi = 'CE', FVEthreshold = 1, maxK = 24))
    #Sig <- pca$fittedCov + diag(pca$rho, nrow = length(pca$workGrid))
    xicov <- lapply(Lt, function(tt) {
      Phi <- ConvertSupport(fromGrid = pca$workGrid, toGrid = tt, phi = pca$phi)
      if (length(tt) == 1) {
        Phi <- t(as.matrix(Phi)) 
      }
      Cov <- ConvertSupport(fromGrid = pca$workGrid, toGrid = as.numeric(tt), Cov = pca$fittedCov)
      Lam <- diag(pca$lambda, nrow = pca$selectK)
      Lam %*% t(Phi) %*% solve(Cov+diag(pca$sigma2, nrow = length(tt))) %*% Phi %*% Lam
    })
    xicov <- (1/length(Lt))*Reduce('+', xicov)
    res[[b]] <- list(xi = pca$xiEst, lam = pca$lambda, xicov = xicov, cumFVE = pca$cumFVE, bw = pca$bwMu, years = years)
  }
  return(res)
}

data <- readRDS("../DataResults/CleanData/Landsat_All_3Indices_5Bands_Filtered_Sparsesd.rds")
ids <- unique(data$pointID)
pcas <- foreach(id = 1:503) %dopar% {
  if (!id %in% ids) {
    return(NULL)
  }
  datai <- subset(data, pointID == id)
  #resi <- GetFPCA(datai, bands = c("B2", "B3", "B4", "B5", "NDVI", "MNDWI", "B7"))
  resi <- GetFPCA(datai, bands = c("fsd_b2", "fsd_b3", "fsd_b4", "fsd_b5", "fsd_b7", "fsd_ndvi", "fsd_mndwi"))
  return(resi)
}
for (b in c("fsd_b2", "fsd_b3", "fsd_b4", "fsd_b5", "fsd_b7", "fsd_ndvi", "fsd_mndwi")) {
  resb <- lapply(pcas, function(x) {x[[b]]})
  saveRDS(resb, file = paste0("../DataResults/FPCAResults/xilam_est_All_uFPCA_", tolower(b), ".rds"))
}
