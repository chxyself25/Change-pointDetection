## some old functions, may not used any more
library(fda)
library(Rcpp)
library(doParallel)
registerDoParallel(cores = 8)
source("./ConvertSupport.R") # there has  a IntpFill function
#sourceCpp("./trapzRcpp.cpp")
sourceCpp("./interp2lin.cpp")

#############################################################################
#####################Partially Observed FPCA#################################
#############################################################################
# Bin the data, M is the number of bins
# M <- 100
# functions for doing FPCA on partially observed data
# input mat is the matrix after densing interpolation: has missing values, doy is its observed points
FPCA_po <- function(mat, doy, year, M, FVEthreshold = 0.9999, maxK = 30) {
  binres <- apply(mat, 1, function(x) {GetBinnedCurve(doy, x, M = M)})
  bin.dat <- matrix(0, ncol = M, nrow = length(year)) 
  for (i in 1:length(year)) {
    bin.dat[i, ] <- binres[[i]]$newy
  }
  bin.doy <- binres[[1]]$midpoint
  # mean function estimation
  mu <- colMeans(bin.dat, na.rm = TRUE)
  # covariance matrix estimation, sigma2 estimation
  cov.mat = cov(bin.dat, use = 'pairwise.complete.obs')
  cov.mat = 0.5 * (cov.mat + t(cov.mat))
  cov.mat[is.na(cov.mat)] <- 0
  ord <- 2
  sigma2 <- mean(diff(t(bin.dat), differences=ord)^2, na.rm=TRUE) / choose(2 * ord, ord)
  # eigenfunction and eigenvalue estimation
  eigres <- GetEigenAnalysisResults(cov.mat-diag(sigma2, M), grid = bin.doy, FVEthreshold = FVEthreshold, maxK = maxK, mu = mu)
  phi <- eigres$phi
  lambda <- eigres$lambda
  # principle component score estimation
  # convert mu and phi to observed points
  mu.obs <- IntpFill(fromGrid = bin.doy, toGrid = doy, y = mu)
  phi.obs <- apply(phi, 2, function(x) {IntpFill(fromGrid = bin.doy, toGrid = doy, y = x)})
  xi.est <- GetINScores(mat, doy, mu.obs, lambda, phi.obs, sigma2)$xiEst
  return(list(lambda = lambda, xi = xi.est, mu = mu.obs, phi = phi.obs, cumFVE = eigres$cumFVE))
}

GetBins <-  function(x,y, xx){
  N <- length(xx)
  count = rep(0,N-1);
  newy = count;
  
  for (i in  2:(N-1)){    
    ids = ((x >= xx[i-1]) &  (x < xx[i]));
    if (all(ids == 0)){
      count[i-1] = 0;
      newy[i-1] = NaN;      
    } else {
      count[i-1] =   sum(ids);
      newy[i-1] =  mean(y[ids]); 
    }
  }
  
  # print('GetBins used')
  # for the last bin, include the left and right end point
  ids =  ((x >= xx[i]) &(x <= xx[i+1]));
  if (all(ids == 0)){
    count[i] = 0;
    newy[i] = NaN;
  }else {
    count[i] = sum(ids);
    newy[i] = mean(y[ids]);
  }
  
  zList = list(newy = newy, count = count)
  return( zList );     
}
GetBinnedCurve <- function(x, y, M = 30, limits = c(min(x), max(x))) {
  h <- diff(limits) / M
  xx <- c(limits[1], 
          seq(limits[1] + h / 2, limits[2] - h / 2, length.out=M - 1), 
          limits[2])
  N <- length(xx);
  midpoint <- seq(limits[1], limits[2], length.out=M)
  
  zList = GetBins(x,y,xx)
  newy = zList$newy # average of y in the binned intervals
  count = zList$count# number of points falling into intervals
  # if (nonEmptyOnly){
  #   midpoint = midpoint[count > 0];
  #   newy = newy[count > 0];
  #   count = count[count > 0];
  #   M = length(midpoint);
  # }
  res = list(midpoint = midpoint, newy = newy, count = count, numBin = M, binWidth = h)
  return(res)
}
GetEigenAnalysisResults <- function(smoothCov, grid, FVEthreshold = 0.9999, maxK = 30, mu = NULL) {
  
  gridSize <- grid[2] - grid[1] # grid should be equally discretized
  numGrids <- nrow(smoothCov)
  
  eig <- eigen(smoothCov)
  
  positiveInd <- eig[['values']] >= 0
  if (sum(positiveInd) == 0) {
    stop('All eigenvalues are negative. The covariance estimate is incorrect.')
  }
  d <- eig[['values']][positiveInd]
  eigenV <- eig[['vectors']][, positiveInd, drop=FALSE]
  
  # in case maxK is more than number of principle components from decomposition
  maxK <- min(maxK, length(d))
  
  d <- d[1:maxK]
  eigenV <- eigenV[, 1:maxK, drop=FALSE]
  FVE <- cumsum(d) / sum(d) * 100
  # final number of component chosen based on FVE
  no_opt <- min(which(FVE >= FVEthreshold * 100))
  
  if (is.null(mu)) {
    mu = 1:dim(eigenV)[1]
  }
  phi <- apply(eigenV, 2, function(x) {
    x <- x / sqrt(trapzRcpp(grid, x^2)) # normalize eigenfunction to make sure L2 norm is 1
    # determine the sign of eigenfunction
    if ( 0 <= sum(x*mu) )
      return(x)
    else
      return(-x)
  })
  lambda <- gridSize * d; # eigen value estimation from the dicretized
  
  fittedCov <- phi %*% diag(x=lambda, nrow = length(lambda)) %*% t(phi)
  
  return(list(lambda = lambda[1:no_opt], phi = phi[,1:no_opt, drop=FALSE], 
              cumFVE = FVE, kChoosen=no_opt, fittedCov=fittedCov))
}
# ymat: n by p matrix of dense regular functional observations 
# t: list of observed time grids for the functional observations
# ret: a list of:
#        xiEst: n by length(lambda) matrix of estimated FPC scores
#        fittedY: n by p matrix of fitted/recovered functional observations
GetINScores <- function(ymat, t, mu, lambda, phi, sigma2=NULL){
  if(length(lambda) != ncol(phi)){
    stop('No. of eigenvalues is not the same as the no. of eigenfunctions.')
  }
  
  n = nrow(ymat)
  tau = sort(unique(signif( unlist(t),14 ))) # get observed time grid
  ranget <- diff(range(tau))
  mumat = matrix(rep(mu, n), nrow = n, byrow = TRUE)
  cymat = ymat - mumat
  
  xiEst = matrix(0, nrow = n, ncol = length(lambda)) 
  # Get Scores xiEst
  for(i in 1:length(lambda)){
    tempmat = cymat * matrix(rep(phi[,i],n), nrow = n, byrow = TRUE)
    xiEst[,i] = sapply(1:n, function(j) trapzRcpp(X = tau[!is.na(tempmat[j,])], Y = tempmat[j, !is.na(tempmat[j,])]))
  }
  
  # Get Fitted Y: n by p matrix on observed time grid
  fittedY = mumat + t(phi %*% t(xiEst))
  
  res = list('xiEst' = xiEst, xiVar = NULL, 'fittedY' = fittedY)
  
  return(res)
}

## supposed to estimate CE scores from data from densing
GetCEScores <- function(dat, mu, lambda, phi, cov, doyeval, band) {
  years <- sort(unique(dat$year))
  xi <- foreach(yr = years, .combine= 'rbind') %dopar% {
    datyr <- subset(dat, year == yr)
    # only observations within doyeval are available for computing xi
    obs.idx <- which(datyr$doy <= max(doyeval) & datyr$doy >= min(doyeval))
    # convert mu, phi, cov to valid observed grids
    mu.obs <- IntpFill(fromGrid = doyeval, toGrid = datyr$doy[obs.idx], y = mu)
    phi.obs <- apply(phi, 2, function(x) {IntpFill(fromGrid = doyeval, toGrid = datyr$doy[obs.idx], y = x)})
    cov.obs <- IntpFill(fromGrid = doyeval, toGrid = datyr$doy[obs.idx], Cov = cov)
    cov.obs <- 0.5*(t(cov.obs)+cov.obs)
    # PACE
    t(diag(lambda) %*% t(phi.obs) %*% solve(cov.obs) %*% (datyr[[band]][obs.idx] - mu.obs))
  }
  return(xi)
}

###########################################################################################
#### abondanded: modified regular one, count from back, actually equivalent to TNforward
TNbackward <- function(x, lambda, xi) {
  N <- nrow(xi)
  d <- length(lambda)
  t <- 0
  for (i in 1:d) {
    if (round(x*N)==0) {
      xi.test <- 0 - x*sum(xi[,i])
    } else {
      xi.test <- sum(xi[(ceiling(x*N)):N,i]) - (1-x)*sum(xi[,i])
      #xi.test <- xi[round(x*N), i] - x*sum(xi[,i])
      #xi.est <- xi[round(x*N), i] - x*(0.5*N)
    }
    t <- t + (1/lambda[i])*(xi.test)
  }
  return(t/N)
}

# test statistics for two change-point detection, x1 < x2
TNDt2 <- function(x1, x2, lambda, xi) {
  N <- nrow(xi)
  d <- length(lambda)
  t <- 0
  for (i in 1:d) {
    xi.test <- sum(xi[1:(floor(x1*N)),i])/x1 + sum(xi[(ceiling(x2*N)):N,i])/x2
    t <- t + (1/lambda[i])*(xi.test^2)
  }
  return(t)
}
thetaN2 <- function(lambda, xi) {
  xx <- seq(0.3, 0.9, length.out = 100) 
  N <- nrow(xi)
  xx.idx <- combn(1:100, m = 2)
  tn <- sapply(1:ncol(xx.idx), function(i) {TNDt2(xx[xx.idx[1,i]], xx[xx.idx[2,i]], lambda, xi)})
  max.idx <- which(tn == max(tn))
  if (length(max.idx) == 1) {ceiling(N*xx[xx.idx[,max.idx]])}
  else {
    idx <- c(max(xx.idx[1,max.idx]), max(xx.idx[2,max.idx]))
    ceiling(N*xx[idx])
  }
}

## estimate change-point by doing regression on the first principle component scores
## xi has to be vector
regthetaN <- function(lambda, xi, precision = 10, type = NULL) {
  index <- 1:length(xi)
  #transxi <- 1/(1+exp(-0.5*xi[,1]))
  lmfit <- lm(xi~index)
  P1 <- which.min(lmfit$residuals)
  P2 <- which.max(lmfit$residuals)
  return(list(P1 = P1, P2 = P2, resid1 = min(lmfit$residuals), resid2 = max(lmfit$residuals)))
}

## application level functions: multivariate FPCA, sigmoid transformation, estimation P1 and P2
## ensemble using difference between two xi at P1 and P2
Cfchange <- function(xi.res, id, year.range = NULL, FVE = 85, beta = NULL) {
  prep <- Prepxi(xi.res, id, year.range, FVE, beta)
  xiest2 <- prep$xi; xiest <- prep$xiold
  lambda <- prep$lam; N <- nrow(xiest2)
  K <- prep$K; yeari <- prep$years
  res <- NULL
  for (k in 1:K) {
    #k.est <- thetaNdc(lambda = lambda[k], xi = xiest2[,k], output = "multiple")
    k.est <- thetaN(lambda = lambda[k], xi = xiest2[,k])
    #cat(k.est$P, "\n")
    #value <- abs(xiest[k.est$P1,k] - xiest[k.est$P2,k])
    value <- k.est$TN
    #cat(value, "\n")
    res <- rbind(res, c(yeari[k.est$P2], value))
  }
  dom.pc <- which.max(res[,2])
  return(res[dom.pc,])
}

# Kd critical values, from berke's paper
Kd <- function(x, a, d) {
  dist <- list('10' = c(0.345165,0.606783,0.842567,1.065349,1.279713,1.485200,
                        1.690773,1.897365,2.096615,2.288572,2.496635,2.686238,
                        2.884214,3.066906,3.268958,3.462039,3.650724,3.837678,
                        4.024313,4.214800,4.404677,4.591972,4.778715,4.965613,
                        5.159057,5.346543,5.521107,5.714145,5.885108,6.083306),
               '5' = c(0.460496,0.748785,1.001390,1.239675,1.469008,1.684729,
                       1.895557,2.124153,2.322674,2.526781,2.744438,2.949004,
                       3.147604,3.336262,3.544633,3.740248,3.949054,4.136169,
                       4.327286,4.532917,4.718904,4.908332,5.101896,5.303462,
                       5.495721,5.688849,5.866095,6.068351,6.242770,6.444772),
               '1' = c(0.740138,1.072101,1.352099,1.626695,1.866702,2.125950,
                       2.342252,2.589244,2.809778,3.033944,3.268031,3.491102,
                       3.708033,3.903995,4.116829,4.317087,4.554650,4.734714,
                       4.974172,5.156282,5.369309,5.576596,5.759427,5.973941,
                       6.203718,6.393582,6.572949,6.771058,6.977607,7.186491))
  if (x < dist[[as.character(a)]][d]) {return(FALSE)}
  else {return(TRUE)}
}
