# change point detection of landsat data including 150 locations and urbanization mapping
# Tue Oct  2 21:27:40 2018 ------------------------------
#library(fdapace)

#################################################################################
#SWIR is helpful in distinguishing urban and bare land
# Band 1 - Blue
# Band 2 - Green
# Band 3 - Red
# Band 4 - Near Infrared
# Band 5 - Middle-infrared ??
# Band 7 - Short-wave Infrared
# NDVI = (NIR − Red)/(NIR + Red) = (4-3)/(4+3)
# MNDWI = (Green−MIR)/(Green+MIR) = (2-5)/(2+5)
#################################################################################

#############################################################################
##########################change-point detection#############################
#############################################################################
# test statistics T_N, regular one
TNforward <- function(x, lambda, xi) {
  d <- length(lambda)
  N <- nrow(xi)
  t <- 0
  for (i in 1:d) {
    if (floor(x*N)==0) {
      xi.test <- 0 - x*sum(xi[,i])
    } else {
      xi.test <- sum(xi[1:(floor(x*N)),i]) - x*sum(xi[,i])
    }
    t <- t + (1/lambda[i])*(xi.test^2)
  }
  return(t/N)
}

#test statistics S_Nd, following K_d distribution
SNd <- function(lambda, xi) {
  xi <- as.matrix(xi) # make sure it's a matrix
  N <- nrow(xi)
  tn <- unlist(sapply(1:N, function(k) {TNforward(k/N, lambda, xi)}))
  return(mean(tn))
}
# estimate change-point using continuous grid
thetaN <- function(lambda, xi, precision = 10) {
  xx <- seq(0, 1, length.out = 1000) 
  xi <- as.matrix(xi) # make sure it's a matrix
  N <- nrow(xi)
  tn <- round(sapply(xx, function(x) {TNforward(x, lambda, xi)}), precision)
  max.idx <- which(tn == max(tn))
  P1 <- max(floor(min(xx[max.idx])*N),1)
  P2 <- ceiling(max(xx[max.idx])*N)
  TN <- max(tn)
  return(list(P1 = P1, P2 = P2, TN = TN))
  #xx[xx.idx[,max.idx]]
}
# estimate change-point using discrete grid
thetaNdc <- function(lambda, xi, output = "single") {
  xi <- as.matrix(xi) # make sure it's a matrix
  N <- nrow(xi)
  xx <- seq(1, N)/N
  tn <- sapply(xx, function(x) {TNforward(x, lambda, xi)})
  max.idx <- which.max(tn) # the first maximum
  TN <- max(tn)
  if (output == "single") {
    return(list(P = max.idx, TN = TN))
  } else {
    if (max.idx == 1 || max.idx == N) {
      max.idx2 <- max.idx + ifelse(max.idx == 1, 1, -1)
    } else {
      cands <- c(tn[max.idx-1], tn[max.idx+1])
      max.idx2 <- max.idx + c(-1,1)[which.max(cands)]
    }
    return(list(P1 = min(max.idx, max.idx2), P2 = max(max.idx, max.idx2), TN = TN))
  }
}

###############################################################################
#################### change-point estimation ##################################
###############################################################################
## preprocessing the xi estimation results
Prepxi <- function(xi.res, id, year.range = NULL, FVE = 85, beta = NULL, K = NULL) {
  xi.resi <- xi.res[[id]]
  yeari <- xi.resi$years
  if (is.null(year.range)) {
    year.range <- range(yeari)
  }
  yr.idx <- which(yeari <=  year.range[2] & yeari >= year.range[1]) # select the range with change-point 
  lambda <- xi.resi$lam
  xiest <- xi.resi$xi[yr.idx,]
  # select number of principles
  if (is.null(K)) { # determine the K by FVE
    K <- min(which(xi.resi$cumFVE >= FVE)) 
  } else {
    K <- min(K, length(lambda))
  }
  xiest2 <- as.matrix(xiest) # make sure it is a matrix
  if (!is.null(beta)) {
    xiest2 <- sapply(1:length(lambda), function(x) {1/(1+exp(-beta*(xiest[,x])))})
  }
  res <- list(xi = xiest2, xiold = xiest, lam = lambda, K = K, years = yeari[yr.idx])
  return(res)
}
## function wraps up the procedure of estimating change-point
## xi.res is the list of FPCA results, is a list consist of xi, lambda, cumulative FVE and years
## pcs is the principle components considered, only applicable for separate testing 
# Mfchange <- function(xi.res, id, year.range = NULL, FVE = 85, beta = NULL, K = NULL,
#                     method = "sum") {
#   prep <- Prepxi(xi.res, id, year.range, FVE, beta, K)
#   xiest <- prep$xi; lambda <- prep$lam
#   K <- prep$K; yeari <- prep$years
#   # calculate statistics, either weighted sum or ensembles
#   if (method == "max") {
#     res <- NULL
#     for (k in 1:K) {
#       k.est <- thetaNdc(lambda = lambda[k], xi = xiest[,k])
#       res <- rbind(res, c(yeari[k.est$P], k.est$TN))
#     }
#     dom.pc <- which.max(res[,2])
#     return(c(res[dom.pc,], K))
#   } else if (method == "sum") {
#     k.est <- thetaNdc(lambda[1:K], xiest[,1:K])
#     return(c(yeari[k.est$P], k.est$TN, K))
#   } else {
#     stop("Please set a valid method for change-point estimation!")
#   }
# }
Mfchange <- function(xi.res, id, year.range = NULL, FVE = 85, beta = NULL, K = NULL,
                     method = "sum") {
  prep <- Prepxi(xi.res, id, year.range, FVE, beta, K)
  xiest <- prep$xi; lambda <- prep$lam
  K <- prep$K; yeari <- prep$years
  # calculate statistics, either weighted sum or ensembles
  if (method == "max") {
    res <- NULL
    for (k in 1:K) {
      k.est <- thetaN(lambda = lambda[k], xi = xiest[,k])
      res <- rbind(res, c(yeari[k.est$P2], k.est$TN))
    }
    dom.pc <- which.max(res[,2])
    return(c(res[dom.pc,], K))
  } else if (method == "sum") {
    k.est <- thetaN(lambda[1:K], xiest[,1:K])
    return(c(yeari[k.est$P2], k.est$TN, K))
  } else {
    stop("Please set a valid method for change-point estimation!")
  }
}

## change-point estimation by univariate FPCA: ensemble three indices or sum three functions
## the first input is a list of different component results, eg. xi.ndvi, xi.mndwi, xi.b7
Ufchange <- function(comps.list, id, year.range = NULL, FVE = 85, beta = NULL, K = NULL, 
                     method = "max") {
  num.comps <- length(comps.list)
  p.list <- rep(list(NULL), num.comps)
  # check if K is given 
  if (!is.null(K)) {
    # determine the K for each index
    Ks <- rep(K%/%num.comps, num.comps)
    Ks <- Ks + sample(c(rep(1, K%%num.comps), rep(0, num.comps - K%%num.comps)), replace = FALSE)
    if (sum(Ks) != K) {
      stop("distribution of K is not correct!")
    }
    for (i in 1:num.comps) {
      p.list[[i]] <- Prepxi(comps.list[[i]], id, year.range, FVE, beta, K = Ks[i]) 
    }
  } else {
    for (i in 1:num.comps) {
      p.list[[i]] <- Prepxi(comps.list[[i]], id, year.range, FVE, beta, K = NULL) # determine K by FVE
    }
  }
  N <- nrow(p.list[[1]]$xi); yeari <- p.list[[1]]$years # same N and years for all indices
  totalK <- sum(sapply(p.list, '[[', 'K'))
  if (method == "max") {
    lambda <- unlist(lapply(p.list, function(x) {x$lam[1:(x$K)]}))
    xiest <- do.call('cbind', lapply(p.list, function(x) {x$xi[, 1:(x$K), drop=FALSE]}))
    res <- NULL
    for (k in totalK) {
      k.est <- thetaNdc(lambda = lambda[k], xi = xiest[,k])
      res <- rbind(res, c(yeari[k.est$P], k.est$TN))
    }
    dom.pc <- which.max(res[,2])
    return(c(res[dom.pc,], totalK))
  } else if (method == "sum") {
    ## get the TN function
    xx <- seq(1, N)/N
    tn <- rep(0, N)
    for (i in 1:num.comps) {
      pi <- p.list[[i]]
      Ki <- pi$K
      tn <- tn + sapply(xx, function(x) {TNforward(x, pi$lam[1:Ki], pi$xi[,1:Ki,drop=FALSE])})
    }
    max.idx <- which.max(tn) # the first maximum
    TN <- max(tn)
    return(c(yeari[max.idx], TN, totalK))
  } else {
    stop("Please set a valid method for change-point estimation!")
  }
}

## function for estimating change-point by regression method
# dat is annual time series
regthetaN <- function(dat, id, year.range = NULL, bands) {
  if (is.null(year.range)) {
    year.range <- range(dat$year[dat$pointID == id])
  }
  dati <- subset(dat, pointID == id & year <= year.range[2] & year >= year.range[1])
  values <- c(); year1s <- c(); year2s <- c()
  for (b in bands) {
    lmfit <- lm(as.formula(paste0(b, "~year")), dati)
    P1 <- which.min(lmfit$residuals)
    P2 <- which.max(lmfit$residuals)
    values <- c(values, abs(dati[[b]][P1] - dati[[b]][P2]))
    year1s <- c(year1s, dati$year[P1])
    year2s <- c(year2s, dati$year[P2])
  }
  type <- which.max(values)
  return(c(pointID = id, year1 = min(year1s[type], year2s[type]), year2 =  max(year1s[type], year2s[type]), value = max(values)))
}


#############################################################################
######################## change-point testing ###############################
#############################################################################
## BSdist is a list of K: 
## each is a ecdf of kth component if max, a ecdf of first k components if sum
Mftest <- function(xi.res, id, year.range = NULL, FVE = 85, beta = NULL, K = NULL,
                   method = "sum", BSdist  = list(), Kd = FALSE, alpha = 0.05) {
  prep <- Prepxi(xi.res, id, year.range, FVE, beta, K)
  xiest <- prep$xi; lambda <- prep$lam
  K <- min(prep$K, length(BSdist))
  # calculate statistics, either weighted sum or ensembles
  if (method == "max") { 
    # BSdist is a list of ecdf, each corresponds to kth component if bootstrap, otherwise cumulative asymptotic distribution
    snds <- sapply(seq_len(K), function(k) {SNd(lambda[k], xiest[,k])})
    snd <- max(snds)
    if (Kd) {
      pval <- 1 - (BSdist[[1]](snd))^K
    } else {
      p_values <- rep(NA, K)
      for (k in 1:K) {
        p_values[k] <- BSdist[[k]](snd)
      }
      pval <- 1-prod(p_values)
    }
    return(c(pval<alpha, K))
  } else if (method == "sum") {
    # BSdist is a list of ecdf, each corresponds to first K components no matter bootstrap or not
    snd <- SNd(lambda[1:K], xiest[,1:K])
    pval <- 1-BSdist[[K]](snd)
    return(c(pval<alpha, K))
  } else {
    stop("Please set a valid method for change-point estimation!")
  }
}

## BS.list is a list of bands or indices, 
## each corresponds to a band/index: list of each component's ecdf if max, 
## bootstrap results, matrix if sum
Uftest <- function(comps.list, id, year.range = NULL, FVE = 85, beta = NULL, K = NULL, 
                   method = "max", BSdist.list = list(), alpha = 0.05) {
  if (length(comps.list) != length(BSdist.list)) {
    stop("The list of pc scores and list of dist. should match!")
  }
  num.comps <- length(comps.list)
  p.list <- rep(list(NULL), num.comps)
  # check if K is given 
  if (!is.null(K)) {
    # determine the K for each index
    Ks <- rep(K%/%num.comps, num.comps)
    Ks <- Ks + sample(c(rep(1, K%%num.comps), rep(0, num.comps - K%%num.comps)), replace = FALSE)
    if (sum(Ks) != K) {
      stop("distribution of K is not correct!")
    }
    for (i in 1:num.comps) {
      p.list[[i]] <- Prepxi(comps.list[[i]], id, year.range, FVE, beta, K = Ks[i]) 
    }
  } else {
    for (i in 1:num.comps) {
      p.list[[i]] <- Prepxi(comps.list[[i]], id, year.range, FVE, beta, K = NULL) # determine K by FVE
    }
  }
  N <- nrow(p.list[[1]]$xi) # same N and years for all indices
  Ks <- sapply(p.list, '[[', 'K')
  totalK <- sum(Ks)
  non.idx <- which(Ks != 0)
  if (method == "max") {
    # BSdist.list is a list of all bands list, each ecdf corresponds to the kth component 
    pvals <- c()
    for (i in non.idx) {
      BSdisti <- BSdist.list[[i]][1:Ks[i]]
      pi <- p.list[[i]]
      sndsi <- sapply(seq(Ks[i]), function(k) {SNd(pi$lam[k], pi$xi[,k])})
      sndi <- max(sndsi)
      p_valuesi <- sapply(seq(Ks[i]), function(k) {BSdisti[[k]](sndi)})
      pvals <- append(pvals, 1-prod(p_valuesi))
    }
    bronf <- alpha/length(non.idx)
    return(c(any(pvals<bronf), totalK))
  } else if (method == "sum") {
    # for the sum case, BSdist.list is the simulation results from each band
    # in case some sampling result in missing values
    min.row <- min(sapply(non.idx, function(i) {nrow(BSdist.list[[i]])}))
    BSres <- rowSums(sapply(non.idx, function(i) {BSdist.list[[i]][1:min.row,Ks[i]]}))
    snd <- sum(sapply(p.list[non.idx], function(x) {SNd(x$lam[1:(x$K)], x$xi[, 1:(x$K)])}))
    pval <- 1-ecdf(BSres)(snd)
    return(c(pval<alpha, totalK))
  } else {
    stop("Please set a valid method for change-point estimation!")
  }
}

## another version of Uftest function that uses simulated Kd distribution 
## for sum ensemble, it is not valid because of dependence among indices
UftestX <- function(comps.list, id, year.range = NULL, FVE = 85, beta = NULL, K = NULL, 
                   Kdist, alpha = 0.05) {
  num.comps <- length(comps.list)
  p.list <- rep(list(NULL), num.comps)
  # check if K is given 
  if (!is.null(K)) {
    # determine the K for each index
    Ks <- rep(K%/%num.comps, num.comps)
    Ks <- Ks + sample(c(rep(1, K%%num.comps), rep(0, num.comps - K%%num.comps)), replace = FALSE)
    if (sum(Ks) != K) {
      stop("distribution of K is not correct!")
    }
    for (i in 1:num.comps) {
      p.list[[i]] <- Prepxi(comps.list[[i]], id, year.range, FVE, beta, K = Ks[i]) 
    }
  } else {
    for (i in 1:num.comps) {
      p.list[[i]] <- Prepxi(comps.list[[i]], id, year.range, FVE, beta, K = NULL) # determine K by FVE
    }
  }
  N <- nrow(p.list[[1]]$xi) # same N and years for all indices
  Ks <- sapply(p.list, '[[', 'K')
  totalK <- sum(Ks)
  non.idx <- which(Ks != 0)
  pvals <- c()
  for (i in non.idx) {
    pi <- p.list[[i]]
    sndsi <- sapply(seq(Ks[i]), function(k) {SNd(pi$lam[k], pi$xi[,k])})
    sndi <- max(sndsi)
    pvals <- append(pvals, 1 - (Kdist[[1]](sndi))^(Ks[i]))
  }
  bronf <- alpha/length(non.idx)
  return(c(any(pvals<bronf), totalK))
}

## function for testing by regression 
# ts: time series of yearly average, data frame of all locations(ids)
# three indices: if one has significant slope, then reject hypothesis
# regTest <- function(ts, id) {
#   tsi <- subset(ts, pointID == id & year <= 2016 & year >= 2001)
#   pvs <- rep(0, 3)
#   bands <- c("NDVI", "MNDWI", "B7")
#   for (i in 1:3) {
#     lmfit <- lm(as.formula(paste0(bands[i], "~year")), data = tsi)
#     pvs[i] <- summary(lmfit)$coefficients[2,4]
#   }
#   return(min(pvs))
# }


## bootstrap sampling using FPCA
## Ly is centered functional data series
fBTSampling <- function(Ly, Lt, K, bs.n = 100, resample.n = 1000) {
  n <- length(Ly)
  b.idx <- sample(1:n, size = bs.n, replace = TRUE)
  Lyb <- Ly[b.idx]; Ltb <- Lt[b.idx]
  pcab <- FPCA(Lyb, Ltb, optns = list(dataType = "Sparse", kernel = "epan", methodBwMu = 'GCV', FVEthreshold = 1,
                                      maxK = K, methodMuCovEst = "smooth", methodXi = "CE"))
  xi <- pcab$xiEst; lambda <- pcab$lambda
  bwMu <- pcab$bwMu; bwCov <- pcab$bwCov
  # calculate Snd for k = 1,..., 24
  res1 <- rep(NA, K)
  res1[seq_len(pcab$selectK)] <- sapply(seq_len(pcab$selectK), function(k) {SNd(lambda[1:k], xi[,1:k])})
  # resample 1000 times
  res2 <- foreach(i = 1:(resample.n-1), .combine = 'rbind') %dopar% {
    b.idx <- sample(1:n, size = bs.n, replace = TRUE)
    Lyb <- Ly[b.idx]
    Ltb <- Lt[b.idx]
    pcab <- FPCA(Lyb, Ltb, optns = list(dataType = "Sparse", userBwMu = bwMu, userBwCov = bwCov, FVEthreshold = 1,
                                        kernel = "epan", maxK = K, methodMuCovEst = "smooth", methodXi = "CE"))
    xi <- pcab$xiEst
    lambda <- pcab$lambda
    # calculate Snd for k = 1,..., 10
    sndK <- rep(NA, K)
    sndK[seq_len(pcab$selectK)] <- sapply(seq_len(pcab$selectK), function(k) {SNd(lambda[1:k], xi[,1:k])})
    sndK
  }
  res <- rbind(res1, res2)
  return(res)
}

## bootstrap sampling by MFPCA
## Ly is centered functional time series
mfBTSampling <- function(Ly, Lt, K, bs.n = 100, resample.n = 1000) {
  n <- length(Ly)
  b.idx <- sample(1:n, size = bs.n, replace = TRUE)
  Lyb <- Ly[b.idx]; Ltb <- Lt[b.idx]
  pcab <- RFPCA(Lyb, Ltb, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean",
                                       kernel = "epan", maxK = K, methodMuCovEst = "smooth", methodXi = "CE"))
  xi <- pcab$xi; lambda <- pcab$lam
  bwMu <- pcab$userBwMu; bwCov <- pcab$userBwCov
  # calculate Snd for k = 1,..., 24
  res1 <- rep(NA, K)
  res1[seq_len(pcab$K)] <- sapply(seq_len(pcab$K), function(k) {SNd(lambda[1:k], xi[,1:k])})
  # resample 1000 times
  res2 <- foreach(i = 1:(resample.n-1), .combine = 'rbind') %dopar% {
    b.idx <- sample(1:n, size = bs.n, replace = TRUE)
    Lyb <- Ly[b.idx]
    Ltb <- Lt[b.idx]
    pcab <- RFPCA(Lyb, Ltb, optns = list(dataType = "Sparse", userBwMu = bwMu, userBwCov = bwCov, mfdName = "Euclidean",
                                         kernel = "epan", maxK = K, methodMuCovEst = "smooth", methodXi = "CE"))
    xi <- pcab$xi; lambda <- pcab$lam
    # calculate Snd for k = 1,..., 10
    sndK <- rep(NA, K)
    sndK[seq_len(pcab$K)] <- sapply(seq_len(pcab$K), function(k) {SNd(lambda[1:k], xi[,1:k])})
    sndK
  }
  res <- rbind(res1, res2)
  return(res)
}
