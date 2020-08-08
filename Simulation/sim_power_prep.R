## simulation: testing power and change size

## create mean functions used in simulation
# dat <- readRDS("./DataResults/CleanData/Landsat_All_3Indices_5Bands_Filtered_Sparsesd.rds")
# dat118 <- subset(dat, pointID == 118+120) # 118 site i has clear change pattern
# grid.len <- 200
# year118 <- unique(dat118$year)
# Ly <- lapply(split(dat118, f = as.factor(dat118$year)), function(x) {
#   t(as.matrix(x[,c("fsd_b7", "fsd_ndvi", "fsd_mndwi")]))})
# Lt <- lapply(split(dat118, f = as.factor(dat118$year)), function(x) {x$doy})
# thetas <- readRDS("./DataResults/ForTesting/change-point_est_All_MSCE5.rds")
# P <- thetas$det[thetas$pointID == 118+120]
# Ly1 <- Ly[year118<=P]; Lt1 <- Lt[year118<=P]
# Ly2 <- Ly[year118>P]; Lt2 <- Lt[year118>P]
# mu1 <- RFPCAMu(Ly1, Lt1, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean",
#                                        kernel = "epan", methodMuCovEst = "smooth", nRegGrid = grid.len))
# mu2 <- RFPCAMu(Ly2, Lt2, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean",
#                                        kernel = "epan", methodMuCovEst = "smooth", nRegGrid = grid.len))
# # make two mean functions on the same grid which is generated from two workGrids
# grid.mu <- seq(max(mu1$workGrid[1], mu2$workGrid[1]), min(mu1$workGrid[grid.len], mu2$workGrid[grid.len]), length.out = grid.len)
# mu1 <- apply(mu1$muWork, 1, function(x) {
#   ConvertSupport(mu1$workGrid, grid.mu, mu = x)
# })
# mu2 <- apply(mu2$muWork, 1, function(x) {
#   ConvertSupport(mu2$workGrid, grid.mu, mu = x)
# })
# saveRDS(list(mu1 = mu1, mu2 = mu2, grid = grid.mu), file = "./Testing/simulation/sim_mu1_mu2.rds")
# 
## create eigen functions used in simulation
# dat <- readRDS("./DataResults/ForTesting/Landsat_497_3Indices_5Bands_Filtered_Sparsesd_centered_MSCE5_fve80.rds")
# dat118 <- subset(dat, pointID == 118+120)
# Ly <- lapply(split(dat118, f = as.factor(dat118$year)), function(x) {
#   t(as.matrix(x[,c("fsd_b7", "fsd_ndvi", "fsd_mndwi")]))})
# Lt <- lapply(split(dat118, f = as.factor(dat118$year)), function(x) {x$doy})
# pca <- RFPCA(Ly, Lt, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean",
#                                    kernel = "epan", maxK = 30, methodMuCovEst = "smooth", methodXi = "CE", nRegGrid = 200))
# saveRDS(list(phi = pca$phi, lambda = pca$lam, grid = pca$workGrid), file = "./Testing/simulation/sim_phis.rds")
# pca$lam[1:3] # 55, 37.5, 3


# VD <- list()
# for (i in 1:10) {
#   VD[[i]] <- randortho(4, type = "orthonormal")/sqrt(10)
# }
# mus <- readRDS("./Testing/simulation/sim_mu1_mu2.rds")
# mu1 <- mus$mu1; diff.func <- mus$mu2 - mus$mu1; mu.grid <- mus$grid
# sample1 <- fChangeSim(mu1, mu.grid, diff.func, c = 0.5, n = 100, K = 20, e.sd = 0.3)
# pca1 <- RFPCA(sample1$Ly, sample1$Lt, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean",
#                                     kernel = "epan", maxK = 24, methodMuCovEst = "smooth", methodXi = "CE"))

## function that simulates functional observation with changes: 18 out of 31
# mu1 is mean function before change-point, f is mu2 - mu1, dimension is 200 by 3
# mu.grid is the corresponding time point of mu1
# phi are three-dimensional eigen functions: 100, 3(variable), K
# first 13 have poisson mean 8, and 18 after
fChangeSim <- function(mu1, mu.grid, diff.func, c, n, K, theta, e.sd = 0.3,
                       decay = "fast", VD = list()) {
  # mean functions
  p <- ncol(mu1); basis.method = ifelse(length(VD) == 0, "split", "random")
  mu2 <- mu1 + c*diff.func
  grid <- seq(1, 365)
  # generate observed points randomly, according to possion distribution
  # the site mimiced has 13 days observed on average
  Lt <- lapply(1:n, function(i) {sort(sample(grid, size = max(2, rpois(1, 13))))})  # in case there is zero obs.
  # convert mean functions to observed points
  cpt <- floor(n*theta)
  mu1obs <- lapply(1:cpt, function(i) {
    apply(mu1, 2, function(x) {ConvertSupport(mu.grid, Lt[[i]], mu = x)})
  })
  mu2obs <- lapply((cpt+1):n, function(i) {
    apply(mu2, 2, function(x) {ConvertSupport(mu.grid, Lt[[i]], mu = x)})
  })
  muobs <- append(mu1obs, mu2obs)
  # eigenfunctions construction
  # fourier basis functions
  if (basis.method == "split") {
    f.basis <- create.fourier.basis(rangeval = c(0, 3*365), nbasis = ifelse(K%%2 != 0, K+2, K+1))
  } else {
    f.basis <- create.fourier.basis(rangeval = c(0, 365), nbasis = ifelse(K%%2 != 0, K+2, K+1))
  }
  # eigenvalues, decay fast
  Ks <- 1:K
  if (decay == "fast") {
    lams <- 150*3^(-Ks)
  } else {
    lams <- 47/(Ks^2)
  }
  # principal component scores
  xis <- mvrnorm(n, mu = rep(0, K), Sigma = diag(lams))
  #xis <- sapply(lams, function(lam) {rexp(n, rate = sqrt(1/lam)) - sqrt(lam)})
  # form the functional time series
  # signv <- sample(c(-1, 1), size = p, prob = c(0.5, 0.5))
  if (basis.method == "split") {
    Ly <- lapply(1:n, function(i) {
      t(muobs[[i]] + do.call("cbind", lapply(1:p, function(j) {
        eval.basis(Lt[[i]]+(j-1)*365, f.basis)[,2:(K+1)] %*% xis[i,1:K]
      })))
    })
  } else {
    Ly <- lapply(1:n, function(i) {
      t(muobs[[i]] + do.call("cbind", lapply(1:p, function(j) {
        eval.basis(Lt[[i]], f.basis)[,2:(K+1)] %*% VD[[j]] %*% xis[i,1:K]
      })))
    })
  }
  Ly <- lapply(Ly, function(y) {
    y + do.call("rbind", lapply(1:p, function(j) {
      rnorm(ncol(y), mean = 0, sd = e.sd)
    }))
  })
  return(list(Ly = Ly, Lt = Lt))
}

### function for doing test for one dataset: MFPCA, calculate Snd
GetSNK <- function(Ly, Lt, FVE) {
  # doing MFPCA
  pca <- RFPCA(Ly, Lt, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean",
                                     kernel = "epan", maxK = 30, methodMuCovEst = "smooth", methodXi = "CE"))
  xiest <- pca$xi
  lam <- pca$lam
  cumFVE <- pca$cumFVE
  # select number of principle components
  K <- min(which(cumFVE > FVE))
  # weighted cusum statistics
  s <- SNd(lam[1:K], xiest[,1:K,drop=FALSE])
  return(list(xiest = xiest, lam = lam, Snk = s, K = K))
}

## function for getting bootstrap distribution
Getboot <- function(Ly, Lt, xiest, lam, K) {
  # estimate change point
  xiest <- as.matrix(xiest) # make sure it is matrix
  N <- nrow(xiest)
  res <- rep(list(NULL), K)
  for (k in 1:K) {
    k.est <- thetaNdc(lambda = lam[k], xi = xiest[,k])
    value <- k.est$TN
    res[[k]] <- c(k.est$P, value)
  }
  d.idx <- which.max(sapply(res, "[", 2))
  cp <- res[[d.idx]][1]
  # estimate two mean functions
  Ly1 <- Ly[1:cp]; Lt1 <- Lt[1:cp]
  Ly2 <- Ly[(cp+1):N]; Lt2 <- Lt[(cp+1):N]
  mui1 <- RFPCAMu(Ly1, Lt1, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean",
                                         kernel = "epan", methodMuCovEst = "smooth"))
  mui2 <- RFPCAMu(Ly2, Lt2, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean",
                                         kernel = "epan", methodMuCovEst = "smooth"))
  # substract mean functions
  resid1 <- lapply(1:length(Ly1), function(i) {
    Ly1[[i]] - mui1$muObs[, match(Lt1[[i]], mui1$obsGrid)]
  })
  resid2 <- lapply(1:length(Ly2), function(i) {
    Ly2[[i]] - mui2$muObs[, match(Lt2[[i]], mui2$obsGrid)]
  })
  resid <- append(resid1, resid2)
  # resampling 500 times
  res <- foreach(i = 1:500, .combine = "rbind") %dopar% {
    b.idx <- sample(1:length(resid), size = 100, replace = TRUE)
    residb <- resid[b.idx]
    Ltb <- Lt[b.idx]
    pcab <- RFPCA(residb, Ltb, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean",
                                         kernel = "epan", maxK = 30, methodMuCovEst = "smooth", methodXi = "CE"))
    xiest <- pcab$xi
    lam <- pcab$lam
    # calculate Snd for k = 1,..., 10
    Sn <- SNd(lam[1:K], xiest[,1:K])
    Sn
  }
  return(res)
}

## function for wrapping up the two procedure: first calculate Snk, and resampling, finally do testing
MWCtest <- function(Ly, Lt, FVE, alpha) {
  s.res <- GetSNK(Ly, Lt, FVE)
  Snk <- s.res$Snk
  dist <- Getboot(Ly, Lt, s.res$xiest, s.res$lam, s.res$K)
  Kdcdf <- ecdf(dist)
  pvalue <- 1-Kdcdf(Snk)
  if (pvalue < alpha) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
