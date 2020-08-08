# simulate theoretical Kd distribution
# run on server
library(doParallel)
library(fdapace, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6")
library(Sim.DiffProc, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6")
registerDoParallel(cores = 50)

###################### simulate Kd distribution #########################################
tt <- seq(0, 1, length.out = 1001)
res <- foreach(i = 1:100000, .combine = "rbind") %dopar% {
  BBs <- matrix(NA, ncol = 25, nrow = 1001)
  for (j in 1:25) {
    BBs[,j] <- BB(N=1000, M=1, t0=0, T=1)
  }
  BBs2 <- apply(BBs, 1, function(x) {cumsum(x^2)})
  apply(BBs2, 1, function(x) {trapzRcpp(tt, x)})
}

cat("simulation finished!", "\n")
saveRDS(res, file = "./Kd_sim_100K_25.rds")
cat("get the empirical cdf...", '\n')
res2 <- apply(res, 2, ecdf)
saveRDS(res2, file = "./Kd_ecdf_25.rds")

