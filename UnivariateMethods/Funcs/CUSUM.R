# functions for calculating Tn and Sn statistics
# test statistics T_N, regular one
# TNforward <- function(x, lambda, xi) {
#   d <- length(lambda)
#   N <- nrow(xi)
#   t <- 0
#   for (i in 1:d) {
#     if (floor(x*N)==0) {
#       xi.test <- 0 - x*sum(xi[,i])
#     } else {
#       xi.test <- sum(xi[1:(floor(x*N)),i]) - x*sum(xi[,i])
#     }
#     t <- t + (1/lambda[i])*(xi.test^2)
#   }
#   return(t/N)
# }
TNforward <- function(x, xicov, xi) {
  K <- nrow(xicov)
  n <- nrow(xi)
  if (floor(x*n)==0) {
    cusum <- rep(0, K) - x*colSums(xi)
  } else {
    cusum <- colSums(xi[1:floor(x*n), ,drop=FALSE]) - x * colSums(xi)
  }
  res <- (1/n) * t(cusum) %*% solve(xicov) %*% cusum
  return(res)
}

# estimate change-point using discrete grid
thetaN <- function(xicov, xi, output = "multiple") {
  xi <- as.matrix(xi) # make sure it's a matrix
  xicov <- as.matrix(xicov)
  n <- nrow(xi)
  if (output == "single") {
    xx <- seq(1, n)/n
    tn <- sapply(xx, function(x) {TNforward(x, xicov, xi)})
    max.idx <- which.max(tn)
    TN <- max(tn)
    return(list(P = max.idx, TN = TN))
  } else {
    xx <- seq(0, 1, length.out = 1000) 
    tn <- sapply(xx, function(x) {TNforward(x, xicov, xi)})
    max.idx <- max(which(tn == max(tn)))
    P1 <- max(floor(min(xx[max.idx])*n),1)
    P2 <- ceiling(max(xx[max.idx])*n)
    TN <- max(tn)
    return(list(P1 = P1, P2 = P2, TN = TN))
  }
}

#test statistics S_Nd, following K_d distribution, xi should be guaranteed to be a matrix
SNd <- function(xicov, xi) {
  xi <- as.matrix(xi) # make sure it is a matrix
  xicov <- as.matrix(xicov)
  n <- nrow(xi)
  tn <- unlist(sapply(1:n, function(k) {TNforward(k/n, xicov, xi)}))
  return(mean(tn))
}

#thetaN(xc[[1]]$xicov[2,2], xc[[1]]$xi[,1])
#SNd(xc[[1]]$xicov[2,2], xc[[1]]$xi[,1])
