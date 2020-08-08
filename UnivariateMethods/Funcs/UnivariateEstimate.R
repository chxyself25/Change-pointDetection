# function for doing change-point estimation (estimating start and end of the urbanization process)
# 3 methods avaliable: MIE, MCE, SCE

## change-point estimation by univariate FPCA: ensemble three indices or sum three functions
## the first input is a list of different component results, eg. ndvi, mndwi, b7
Ufchange <- function(comps.list, year.range = NULL, FVE = 85, K = NULL, 
                     method = "SCE") {
  p.list <- PrepFPCA(comps.list, year.range, FVE, K)
  N <- nrow(p.list[[1]]$xi); yeari <- p.list[[1]]$years # same N and years for all indices
  Ks <- sapply(p.list, '[[', 'K')
  totalK <- sum(Ks) # the totalK might be less than K given
  non.idx <- which(Ks != 0)
  if (method == "IME") {
    lambda <- unlist(lapply(p.list[non.idx], function(x) {diag(x$xicov)[1:(x$K)]}))
    xiest <- do.call('cbind', lapply(p.list[non.idx], function(x) {x$xi[, 1:(x$K), drop=FALSE]}))
    res <- NULL
    for (k in totalK) {
      k.est <- thetaN(xicov = lambda[k], xi = xiest[,k])
      res <- rbind(res, c(yeari[k.est$P1], yeari[k.est$P2], k.est$TN))
    }
    dom.pc <- which.max(res[,3])
    return(c(res[dom.pc,], totalK)) 
  } else if (method == "CSE") {
    xx <- seq(0, 1, length.out = 1000)
    tn <- rep(0, 1000)
    for (pi in p.list[non.idx]) {
      Ki <- pi$K
      tn <- tn + sapply(xx, function(x) {TNforward(x, pi$xicov[1:Ki, 1:Ki, drop=FALSE], pi$xi[,1:Ki,drop=FALSE])})
    }
    max.idx <- max(which(tn == max(tn)))
    P1 <- max(floor(min(xx[max.idx])*N),1)
    P2 <- ceiling(max(xx[max.idx])*N)
    TN <- max(tn)
    return(c(yeari[P1], yeari[P2], TN, totalK))
  } else if (method == "ISE"){
    res <- NULL
    for (pi in p.list[non.idx]) {
      Ki <- pi$K
      i.est <- thetaN(xicov = pi$xicov[1:Ki, 1:Ki], xi = pi$xi[, 1:Ki])
      res <- rbind(res, c(yeari[i.est$P1], yeari[i.est$P2], i.est$TN))
    }
    dom.idx <- which.max(res[,3])
    return(c(res[dom.idx,], totalK))
  }else {
    stop("Please set a valid method for change-point estimation!")
  }
}


## test passed: test sample: pixel 1
# xi.res <- readRDS("./UnivariateMethods/DSM_test_fpca.rds")
# Kd_dist <- readRDS("./DataResults/ForTesting/Kd_ecdf_25.rds")
# Ufchange(xi.res, year.range = NULL, FVE = 85, K = NULL, method = 'MCE')
