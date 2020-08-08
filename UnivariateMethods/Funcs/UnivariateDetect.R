# function for doing urban detection (hypothesis testing)
# 3 methods avaliable: MCE, MIE, SCE


# change-point detection function using asymptotic distribution
# SCE is not valid here, because asymptotic distribution cannot be obtained
# input is list of FPCA results, correponding to different bands
UftestA <- function(comps.list, year.range = NULL, FVE = 85, K = NULL,
                    Kdist, alpha = 0.05, method = 'IME') {
  p.list <- PrepFPCA(comps.list, year.range, FVE, K, Kd.dim = length(Kdist))
  N <- nrow(p.list[[1]]$xi) # same N and years for all indices
  Ks <- sapply(p.list, '[[', 'K')
  totalK <- sum(Ks)
  non.idx <- which(Ks != 0)
  pvals <- c()
  if (method == 'IME') {
    for (i in non.idx) {
      pi <- p.list[[i]]
      sndsi <- sapply(seq(Ks[i]), function(k) {SNd(pi$xicov[k,k], pi$xi[,k])})
      sndi <- max(sndsi)
      pvals <- append(pvals, 1 - (Kdist[[1]](sndi))^(Ks[i]))
    }
  } else if (method == 'ISE') {
    for (i in non.idx) {
      pi <- p.list[[i]]
      Ki <- Ks[i]
      sndi <- SNd(pi$xicov[1:Ki, 1:Ki], pi$xi[, 1:Ki])
      pvals <- append(pvals, 1 - Kdist[[Ki]](sndi))
    }
  } else {
    stop("The input detection method is invalid!")
  }
  bronf <- alpha/length(non.idx)
  return(c(any(pvals<bronf), totalK))
  #return(min(pvals), totalK)
}

## test passed: test sample, pixel 1
# xi.res <- readRDS("./UnivariateMethods/DSM_test_fpca.rds")
# Kd_dist <- readRDS("./DataResults/ForTesting/Kd_ecdf_25.rds")
# UftestA(xi.res, year.range = NULL, FVE = 85, K = NULL, Kd_dist, alpha = 0.05, method = 'MIE')
