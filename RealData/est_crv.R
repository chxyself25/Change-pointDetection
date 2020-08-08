source("../Funcs/urban_detect_correct.R")
source("../Funcs/eval_func.R")
source("./est_data_load.R")
library(doParallel)
library(dplyr, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6")
registerDoParallel(cores = 50)
FVEs <- c(26:99, 99.5, 99.9)
 
metric_func <- function(det, year) {
  mean(abs(det-year))
}
methods <- c('mMAX', 'mSUM', 'uMAX', 'uSUM')
dat <- readRDS("./Estimate_P2_All_FVEs.rds")
res <- NULL
for (m in methods) {
  cat(m, 3, "\n")
  datm <- subset(dat, method == m & Dimension == 3)
  resm <- CVfchange(ids, splits, datm, "fve", FVEs, n = 500)
  res <- rbind(res, data.frame(MAE = resm, method = m, Dimension = 3))
  cat(m, 5, "\n")
  datm <- subset(dat, method == m & Dimension == 5)
  resm <- CVfchange(ids, splits, datm, "fve", FVEs, n = 500)
  res <- rbind(res, data.frame(MAE = resm, method = m, Dimension = 5))
}
rownames(res) <- NULL
saveRDS(res, file = "./CRV_P2_All_FVEs_MAE.rds")

metric_func <- function(det, year) {
  sqrt(mean((det-year)^2))
}
res <- NULL
for (m in methods) {
  cat(m, 3, "\n")
  datm <- subset(dat, method == m & Dimension == 3)
  resm <- CVfchange(ids, splits, datm, "fve", FVEs, n = 500)
  res <- rbind(res, data.frame(RMSE = resm, method = m, Dimension = 3))
  cat(m, 5, "\n")
  datm <- subset(dat, method == m & Dimension == 5)
  resm <- CVfchange(ids, splits, datm, "fve", FVEs, n = 500)
  res <- rbind(res, data.frame(RMSE = resm, method = m, Dimension = 5))
}
rownames(res) <- NULL
saveRDS(res, file = "./CRV_P2_All_FVEs_RMSE.rds")

metric_func <- function(det, year) {
  mean(abs(det-year) > 1)
}
res <- NULL
for (m in methods) {
  cat(m, 3, "\n")
  datm <- subset(dat, method == m & Dimension == 3)
  resm <- CVfchange(ids, splits, datm, "fve", FVEs, n = 500)
  res <- rbind(res, data.frame(Acc = resm, method = m, Dimension = 3))
  cat(m, 5, "\n")
  datm <- subset(dat, method == m & Dimension == 5)
  resm <- CVfchange(ids, splits, datm, "fve", FVEs, n = 500)
  res <- rbind(res, data.frame(Acc = resm, method = m, Dimension = 5))
}
rownames(res) <- NULL
saveRDS(res, file = "./CRV_P2_All_FVEs_Acc.rds")




