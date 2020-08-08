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
TNforward <- function(x, xicov, xi) {
  K <- nrow(xicov)
  n <- nrow(xi)
  if (floor(x*n)==0) {
    cusum <- rep(0, K) - x*colSums(xi)
  } else {
    cusum <- colSums(xi[1:floor(x*n), ,drop=FALSE]) - x * colSums(xi)
  }
  res <- tryCatch((1/n) * t(cusum) %*% solve(xicov) %*% cusum, error = function(e) {NA})
  return(res)
}

#test statistics S_Nd, following K_d distribution
SNd <- function(xicov, xi) {
  xi <- as.matrix(xi) # make sure it is a matrix
  xicov <- as.matrix(xicov)
  n <- nrow(xi)
  tn <- unlist(sapply(1:n, function(k) {TNforward(k/n, xicov, xi)}))
  return(mean(tn))
}

# estimate change-point using continuous grid
thetaN <- function(xicov, xi, precision = 10) {
  xx <- seq(0, 1, length.out = 1000) 
  xi <- as.matrix(xi) # make sure it's a matrix
  xicov <- as.matrix(xicov)
  n <- nrow(xi)
  tn <- round(sapply(xx, function(x) {TNforward(x, xicov, xi)}), precision)
  max.idx <- which(tn == max(tn))
  P1 <- max(floor(min(xx[max.idx])*n),1)
  P2 <- ceiling(max(xx[max.idx])*n)
  TN <- max(tn)
  return(list(P1 = P1, P2 = P2, TN = TN))
  #xx[xx.idx[,max.idx]]
}
# estimate change-point using discrete grid
thetaNdc <- function(xicov, xi, output = "single") {
  xi <- as.matrix(xi) # make sure it's a matrix
  xicov <- as.matrix(xicov)
  n <- nrow(xi)
  if (output == "single") { 
    xx <- seq(1, n)/n
    tn <- sapply(xx, function(x) {TNforward(x, xicov, xi)})
    max.idx <- which.max(tn) # the first maximum
    TN <- max(tn)
    return(list(P = max.idx, TN = TN))
  } else { # used in change-point estimation for real data
    xx <- seq(0, 1, length.out = 1000)
    tn <- round(sapply(xx, function(x) {TNforward(x, xicov, xi)}), 10)
    max.idx <- which(tn == max(tn))
    P2 <- ceiling(max(xx[max.idx])*n)
    TN <- max(tn)
    return(list(P = P2, TN = TN))
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
  res <- list(xi = xiest2, xiold = xiest, lam = lambda, xicov = xi.resi$xicov, K = K, years = yeari[yr.idx])
  return(res)
}

Mfchange <- function(xi.res, id, year.range = NULL, FVE = 85, beta = NULL, K = NULL,
                     method = "sum", est.output = "single") {
  prep <- Prepxi(xi.res, id, year.range, FVE, beta, K)
  xiest <- prep$xi; xicov <- prep$xicov
  K <- prep$K; yeari <- prep$years
  # calculate statistics, either weighted sum or ensembles
  if (method == "max") {
    res <- NULL
    for (k in 1:K) {
      k.est <- thetaNdc(xicov = xicov[k,k], xi = xiest[,k], output = est.output)
      res <- rbind(res, c(yeari[k.est$P], k.est$TN))
    }
    dom.pc <- which.max(res[,2])
    return(c(res[dom.pc,], K))
  } else if (method == "sum") {
    k.est <- thetaNdc(xicov = xicov[1:K, 1:K], xi = xiest[,1:K], output = est.output)
    return(c(yeari[k.est$P], k.est$TN, K))
  } else {
    stop("Please set a valid method for change-point estimation!")
  }
}

## change-point estimation by univariate FPCA: ensemble three indices or sum three functions
## the first input is a list of different component results, eg. xi.ndvi, xi.mndwi, xi.b7
Ufchange <- function(comps.list, id, year.range = NULL, FVE = 85, beta = NULL, K = NULL, 
                     method = "max", est.output = "single") {
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
    lambda <- unlist(lapply(p.list, function(x) {diag(x$xicov)[1:(x$K)]}))
    xiest <- do.call('cbind', lapply(p.list, function(x) {x$xi[, 1:(x$K), drop=FALSE]}))
    res <- NULL
    for (k in 1:totalK) {
      k.est <- thetaNdc(xicov = lambda[k], xi = xiest[,k], output = est.output)
      res <- rbind(res, c(yeari[k.est$P], k.est$TN))
    }
    dom.pc <- which.max(res[,2])
    return(c(res[dom.pc,], totalK))
  } else if (method == "sum") {
    ## get the TN function
    if (est.output == "single") {
      xx <- seq(1, N)/N
      tn <- rep(0, N)
    } else {
      xx <- seq(0, 1, length.out = 1000)
      tn <- rep(0, 1000)
    }
    for (i in 1:num.comps) {
      pi <- p.list[[i]]
      Ki <- pi$K
      tn <- tn + sapply(xx, function(x) {TNforward(x, pi$xicov[1:Ki, 1:Ki, drop=FALSE], pi$xi[,1:Ki,drop=FALSE])})
    }
    if (est.output == "single") {
      max.idx <- which.max(tn) # the first maximum
      TN <- max(tn)
      return(c(yeari[max.idx], TN, totalK))
    } else {
      max.idx <- which(tn == max(tn))
      P2 <- ceiling(max(xx[max.idx])*n)
      TN <- max(tn)
      return(c(yeari[P2], TN, totalK))
    }
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
  xiest <- prep$xi; xicov <- prep$xicov
  K <- min(prep$K, length(BSdist))
  # calculate statistics, either weighted sum or ensembles
  if (method == "max") { 
    # BSdist is a list of ecdf, each corresponds to kth component if bootstrap, otherwise cumulative asymptotic distribution
    snds <- sapply(seq_len(K), function(k) {SNd(xicov[k,k], xiest[,k])})
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
    #return(c(pval<alpha, K))
    return(c(pval, K))
  } else if (method == "sum") {
    # BSdist is a list of ecdf, each corresponds to first K components no matter bootstrap or not
    snd <- SNd(xicov[1:K, 1:K], xiest[,1:K])
    pval <- 1-BSdist[[K]](snd)
    #return(c(pval<alpha, K))
    return(c(pval, K))
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
  if (method == "max") {
    # BSdist.list is a list of all bands list, each ecdf corresponds to the kth component 
    Ks <- sapply(1:num.comps, function(i) {min(p.list[[i]]$K, length(BSdist.list[[i]]))})		
    totalK <- sum(Ks)		
    non.idx <- which(Ks != 0)
    pvals <- c()
    for (i in non.idx) {
      BSdisti <- BSdist.list[[i]][1:Ks[i]]
      pi <- p.list[[i]]
      sndsi <- sapply(seq(Ks[i]), function(k) {SNd(pi$xicov[k,k], pi$xi[,k])})
      sndi <- max(sndsi)
      p_valuesi <- sapply(seq(Ks[i]), function(k) {BSdisti[[k]](sndi)})
      pvals <- append(pvals, 1-prod(p_valuesi))
    }
    bronf <- alpha/length(non.idx)
    #return(c(any(pvals<bronf), totalK))
    return(c(min(pvals), totalK))
  } else if (method == "sum") {
    # for the sum case, BSdist.list is the simulation results from each band
    # in case some sampling result in missing values
    Ks <- sapply(1:num.comps, function(i) {min(p.list[[i]]$K, ncol(BSdist.list[[i]]))})		
    totalK <- sum(Ks)		
    non.idx <- which(Ks != 0)
    min.row <- min(sapply(non.idx, function(i) {nrow(BSdist.list[[i]])}))
    BSres <- rowSums(sapply(non.idx, function(i) {BSdist.list[[i]][1:min.row, Ks[i]]}))		 
    snd <- sum(sapply(non.idx, function(i) {
      pi <- p.list[[i]]		
      SNd(pi$xicov[1:Ks[i], 1:Ks[i]], pi$xi[, 1:Ks[i]])		
    }))
    pval <- 1-ecdf(BSres)(snd)
    #return(c(pval<alpha, totalK))
    return(c(pval, totalK))
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
    sndsi <- sapply(seq(Ks[i]), function(k) {SNd(pi$xicov[k,k], pi$xi[,k])})
    sndi <- max(sndsi)
    # if (sndi > sndi.max) {
    #   sndi.max <- sndi 
    # }
    pvals <- append(pvals, 1 - (Kdist[[1]](sndi))^(Ks[i]))
  }
  bronf <- alpha/length(non.idx)
  #return(c(any(pvals<bronf), totalK))
  return(c(min(pvals), totalK))
}

## function for testing by regression 
# ts: time series of yearly average, data frame of all locations(ids)
# three indices: if one has significant slope, then reject hypothesis
regTest <- function(dat, id, year.range = NULL, bands, alpha = 0.05) {
  if (is.null(year.range)) {
    year.range <- range(dat$year[dat$pointID == id])
  }
  dati <- subset(dat, pointID == id & year <= year.range[2] & year >= year.range[1])
  pvs <- NULL
  for (b in bands) {
    lmfit <- lm(as.formula(paste0(b, "~year")), dati)
    pvs <- c(pvs, summary(lmfit)$coefficients[2,4])
  }
  bronf <- alpha/length(bands)
  return(c(pointID = id, det = any(pvs < bronf)))
}


## bootstrap sampling using FPCA
## Ly is centered functional data series
fBTSampling <- function(Ly, Lt, K, bs.n = 100, resample.n = 1000) {
  n <- length(Ly)
  b.idx <- sample(1:n, size = bs.n, replace = TRUE)
  Lyb <- Ly[b.idx]; Ltb <- Lt[b.idx]
  pcab <- NA
  while (!is.list(pcab)) {
    pcab <- tryCatch(FPCA(Lyb, Ltb, optns = list(dataType = "Sparse", kernel = "epan", methodBwMu = "GCV", methodBwCov = "GCV",
                                      maxK = K, methodMuCovEst = "smooth", methodXi = "CE", FVEthreshold = 1)),
                     error = function(e) {NA})
  }
  xi <- pcab$xiEst
  bwMu <- pcab$bwMu; bwCov <- pcab$bwCov
  xicov <- lapply(Ltb, function(tt) {
    Phi <- ConvertSupport(fromGrid = pcab$workGrid, toGrid = tt, phi = pcab$phi)
    if (length(tt) == 1) {
      Phi <- t(Phi)
    }
    Cov <- ConvertSupport(fromGrid = pcab$workGrid, toGrid = as.numeric(tt), Cov = as.matrix(pcab$fittedCov))
    Lam <- diag(pcab$lambda, nrow = pcab$selectK)
    Lam %*% t(Phi) %*% solve(Cov+diag(pcab$sigma2, nrow = length(tt))) %*% Phi %*% Lam
  })
  xicov <- (1/length(Lt))*Reduce('+', xicov)
  # calculate Snd for k = 1,..., 24
  cat("initial fBTS", '\n')
  res1 <- list(S = rep(NA, K), A = rep(NA, K), lams = rep(NA, K))
  res1$S[seq_len(pcab$selectK)] <- sapply(seq_len(pcab$selectK), function(k) {SNd(xicov[1:k, 1:k], xi[,1:k])})
  res1$A[seq_len(pcab$selectK)] <- sapply(seq_len(pcab$selectK), function(k) {SNd(xicov[k, k], xi[,k])})
  res1$lams[seq_len(pcab$selectK)] <- pcab$lambda
  # resample 1000 times
  cat("real fBTS", '\n')
  res2 <- foreach(i = 1:(resample.n-1)) %dopar% {
    b.idx <- sample(1:n, size = bs.n, replace = TRUE)
    Lyb <- Ly[b.idx]
    Ltb <- Lt[b.idx]
    pcab <- tryCatch(FPCA(Lyb, Ltb, optns = list(dataType = "Sparse", kernel = "epan", userBwMu = bwMu, userBwCov = bwCov,
                          maxK = K, methodMuCovEst = "smooth", methodXi = "CE", FVEthreshold = 1)), 
                     error = function(e) {NA})
    if (!is.list(pcab)) {return(NULL)}
    xicov <- lapply(Ltb, function(tt) {
      Phi <- ConvertSupport(fromGrid = pcab$workGrid, toGrid = tt, phi = pcab$phi)
      if (length(tt) == 1) {
        Phi <- t(Phi)
      }
      Cov <- ConvertSupport(fromGrid = pcab$workGrid, toGrid = as.numeric(tt), Cov = as.matrix(pcab$fittedCov))
      Lam <- diag(pcab$lambda, nrow = pcab$selectK)
      Lam %*% t(Phi) %*% solve(Cov+diag(pcab$sigma2, nrow = length(tt))) %*% Phi %*% Lam
    })
    xicov <- (1/length(Lt))*Reduce('+', xicov)
    xi <- pcab$xiEst
    # calculate Snd for k = 1,..., 10
    res2i <- list(S = rep(NA, K), A = rep(NA, K), lams = rep(NA, K))
    res2i$S[seq_len(pcab$selectK)] <- sapply(seq_len(pcab$selectK), function(k) {SNd(xicov[1:k, 1:k], xi[,1:k])})
    res2i$A[seq_len(pcab$selectK)] <- sapply(seq_len(pcab$selectK), function(k) {SNd(xicov[k, k], xi[,k])})
    res2i$lams[seq_len(pcab$selectK)] <- pcab$lambda
    return(res2i)
  }
  cat("all fBTS finished", '\n')
  #saveRDS(list(res1, res2), file = "./debug_fBTS.rds")
  res <- list()
  res$S <- rbind(res1$S, do.call('rbind', lapply(res2, function(x) {x$S})))
  res$A <- rbind(res1$A, do.call('rbind', lapply(res2, function(x) {x$A})))
  res$lams <- rbind(res1$lams, do.call('rbind', lapply(res2, function(x) {x$lams})))
  return(res)
}

## bootstrap sampling by MFPCA
## Ly is centered functional time series
mfBTSampling <- function(Ly, Lt, K, bs.n = 100, resample.n = 1000) {
  n <- length(Ly); p <- nrow(Ly[[1]])
  b.idx <- sample(1:n, size = bs.n, replace = TRUE)
  Lyb <- Ly[b.idx]; Ltb <- Lt[b.idx]
  pcab <- NA
  while (!is.list(pcab)) {
    pcab <- tryCatch(RFPCA(Lyb, Ltb, optns = list(dataType = "Sparse", userBwMu = "GCV", userBwCov = "GCV", mfdName = "Euclidean",
                           kernel = "epan", maxK = K, methodMuCovEst = "smooth", methodXi = "CE")),
                     error = function(e) {NA})
  }
  xi <- pcab$xi; lambda <- pcab$lam
  bwMu <- pcab$userBwMu; bwCov <- pcab$userBwCov
  xicov <- lapply(Ltb, function(tt) {
    m <- length(tt)
    Phi <- apply(pcab$phiObsTrunc[match(tt, pcab$obsGrid),,,drop=FALSE], 3, function(x) {c(t(x))})
    Lam <- diag(pcab$lam, nrow = pcab$K)
    Cov0 <- pcab$covObs[match(tt, pcab$obsGrid), match(tt, pcab$obsGrid),,,drop=FALSE]
    Cov <- matrix(NA, nrow = m*p, ncol = m*p)
    for (i in 1:m) {
      indi <- (p*(i-1)+1):(p*i)
      for (j in 1:m) {
        indj <- (p*(j-1)+1):(p*j)
        Cov[indi, indj] <- Cov0[i,j,,] 
      }
    }
    Lam %*% t(Phi) %*% solve(Cov+diag(pcab$sigma2, nrow = m*p)) %*% Phi %*% Lam
  })
  xicov <- (1/length(Ltb))*Reduce('+', xicov)
  # calculate Snd for k = 1,..., 24
  res1 <- list(S = rep(NA, K), A = rep(NA, K), lams = rep(NA, K))
  res1$S[seq_len(pcab$K)] <- sapply(seq_len(pcab$K), function(k) {SNd(xicov[1:k, 1:k], xi[,1:k])})
  res1$A[seq_len(pcab$K)] <- sapply(seq_len(pcab$K), function(k) {SNd(xicov[k, k], xi[,k])})
  res1$lams[seq_len(pcab$K)] <- pcab$lam
  # resample 1000 times
  res2 <- foreach(i = 1:(resample.n-1)) %dopar% {
    b.idx <- sample(1:n, size = bs.n, replace = TRUE)
    Lyb <- Ly[b.idx]
    Ltb <- Lt[b.idx]
    pcab <- tryCatch(RFPCA(Lyb, Ltb, optns = list(dataType = "Sparse", userBwMu = bwMu, userBwCov = bwCov, mfdName = "Euclidean",
                           kernel = "epan", maxK = K, methodMuCovEst = "smooth", methodXi = "CE")),
                     error = function(e) {NA})
    if (!is.list(pcab)) {
      return(NULL)
    }
    xi <- pcab$xi; lambda <- pcab$lam
    xicov <- lapply(Ltb, function(tt) {
      m <- length(tt)
      Phi <- apply(pcab$phiObsTrunc[match(tt, pcab$obsGrid),,,drop=FALSE], 3, function(x) {c(t(x))})
      Lam <- diag(pcab$lam, nrow = pcab$K)
      Cov0 <- pcab$covObs[match(tt, pcab$obsGrid), match(tt, pcab$obsGrid),,,drop=FALSE]
      Cov <- matrix(NA, nrow = m*p, ncol = m*p)
      for (i in 1:m) {
        indi <- (p*(i-1)+1):(p*i)
        for (j in 1:m) {
          indj <- (p*(j-1)+1):(p*j)
          Cov[indi, indj] <- Cov0[i,j,,]
        }
      }
      Lam %*% t(Phi) %*% solve(Cov+diag(pcab$sigma2, nrow = m*p)) %*% Phi %*% Lam
    })
    xicov <- (1/length(Ltb))*Reduce('+', xicov)
    # calculate Snd for k = 1,..., 10
    res2i <- list(S = rep(NA, K), A = rep(NA, K), lams = rep(NA, K))
    res2i$S[seq_len(pcab$K)] <- sapply(seq_len(pcab$K), function(k) {SNd(xicov[1:k, 1:k], xi[,1:k])})
    res2i$A[seq_len(pcab$K)] <- sapply(seq_len(pcab$K), function(k) {SNd(xicov[k, k], xi[,k])})
    res2i$lams[seq_len(pcab$K)] <- pcab$lam
    return(res2i)
  }
  res <- list()
  res$S <- rbind(res1$S, do.call('rbind', lapply(res2, function(x) {x$S})))
  res$A <- rbind(res1$A, do.call('rbind', lapply(res2, function(x) {x$A})))
  res$lams <- rbind(res1$lams, do.call('rbind', lapply(res2, function(x) {x$lams})))
  return(res)
}
