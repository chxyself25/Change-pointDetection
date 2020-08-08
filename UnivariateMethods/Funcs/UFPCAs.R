# wrap up function for producing PC scores, lambda, etc for change-point statistics
# doing FPCA for each index, save intermediate results if necessary
# prepare the input for doing change-point detection
# require; fdapace

# bands are indices/bands for univariate FPCA
# dat refers to one pixel's observation
# return list of bands, each one is a list of FPCA results
GetFPCA <- function(dat, bands) {
  res <- list()
  years <- unique(dat$year)
  for (b in bands) {
    Ly <- split(dat[[b]], f = as.factor(dat$year))
    Lt <- split(dat$doy, f = as.factor(dat$year))
    pca <- FPCA(Ly, Lt, optns = list(dataType = 'Sparse', kernel = 'epan', methodBwMu = 'GCV', methodBwCov = 'GCV',
                              methodMuCovEst = 'smooth', methodXi = 'CE', FVEthreshold = 1))
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
    res[[b]] <- list(xi = pca$xiEst, lam = pca$lambda, xicov = xicov, cumFVE = pca$cumFVE, years = years)
  }
  return(res)
}

## preprocessing the FPCA results from one band
# year.range is the time restriction for doing detection and estimation
# K is the number of components given 
Prepxi <- function(xi.res, year.range = NULL, FVE = 85, K = NULL) {
  yeari <- xi.res$years
  if (is.null(year.range)) {
    year.range <- range(yeari)
  }
  yr.idx <- which(yeari <=  year.range[2] & yeari >= year.range[1]) # select the range with change-point 
  lambda <- xi.res$lam
  xiest <- xi.res$xi[yr.idx,]
  # select number of principal components
  if (is.null(K)) { # determine the K by FVE
    K <- min(which(xi.res$cumFVE >= FVE)) 
  } else {
    K <- min(K, length(lambda)) # make sure it does not exceed total 
  }
  res <- list(xi = xiest, lam = lambda, xicov = xi.res$xicov, K = K, years = yeari[yr.idx])
  return(res)
}

PrepFPCA <- function(comps.list, year.range = NULL, FVE = 85, K = NULL, Kd.dim = 25) {
  num.comps <- length(comps.list)
  p.list <- rep(list(NULL), num.comps)
  # check if K is given 
  if (!is.null(K)) {
    K <- min(K, Kd.dim*num.comps)
    # determine the K for each index
    Ks <- rep(K%/%num.comps, num.comps)
    Ks <- Ks + sample(c(rep(1, K%%num.comps), rep(0, num.comps - K%%num.comps)), replace = FALSE)
    if (sum(Ks) != K) {
      stop("distribution of K is not correct!")
    }
    for (i in 1:num.comps) {
      p.list[[i]] <- Prepxi(comps.list[[i]], year.range, FVE, K = Ks[i]) 
    }
  } else {
    for (i in 1:num.comps) {
      p.list[[i]] <- Prepxi(comps.list[[i]], year.range, FVE, K = NULL) # determine K by FVE
    }
  }
  return(p.list)
}

## test passed: test sample, pixel 1, 100
# library(signal)
# library(dplyr)
# library(fdapace)
# test_dat <- readRDS("./UnivariateMethods/DSM_test.rds")
# source("./UnivariateMethods/Funcs/Preprocess.R")
# data = Preprocess(test_dat, var.names = c('name', 'doy', 'year', 'B1', 'B2', 'B3', 'B4', 'B5', 'B7'))
# datai <- subset(data, pointID == 1)
# pcai <- GetFPCA(datai, bands = c('ndvi', 'mndwi', 'b7'))
# xc <- PrepFPCA(pcai, year.range = NULL, FVE = 85, K = NULL)



