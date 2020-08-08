source("../Funcs/urban_detect_correct.R")
source("../Funcs/eval_func.R")
source("./est_data_load.R")
library(doParallel)
registerDoParallel(cores = 50)
#FVEs <- c(46:99, 99.5, 99.9)
FVEs <- c(26:45)

## estimation using functional methods
if (FALSE) {
res <- readRDS("./Estimate_P2_All_FVEs.rds")
for (f in FVEs) {
  cat("doing ", f, "\n")
  for (b in c(3,5)) {
  cat("mMAX", "\n")
  resf <- fchange_multi(get(paste0("mxi.", b)), ids = ids, year.range = c(2001, 2016), FVE = f, beta = NULL, K = NULL,
                      method = "max", type = "mv", est.output = "multiple")
  resf <- cbind(resf, data.frame(fve = f, year = ref$year, method = "mMAX", Dimension = b))
  res <- rbind(res, resf)
  cat("mSUM", "\n")
  resf <- fchange_multi(get(paste0("mxi.", b)), ids = ids, year.range = c(2001, 2016), FVE = f, beta = NULL, K = NULL,
                      method = "sum", type = "mv", est.output = "multiple")
  resf <- cbind(resf, data.frame(fve = f, year = ref$year, method = "mSUM", Dimension = b))
  res <- rbind(res, resf)
  cat("uMAX", "\n")
  resf <- fchange_multi(get(paste0("xi.u", b)), ids = ids, year.range = c(2001, 2016), FVE = f, beta = NULL, K = NULL,
                      method = "max", type = "uv", est.output = "multiple")
  resf <- cbind(resf, data.frame(fve = f, year = ref$year, method = "uMAX", Dimension = b))
  res <- rbind(res, resf)
  cat("uSUM", "\n")
  resf <- fchange_multi(get(paste0("xi.u", b)), ids = ids, year.range = c(2001, 2016), FVE = f, beta = NULL, K = NULL,
                      method = "sum", type = "uv", est.output = "multiple")
  resf <- cbind(resf, data.frame(fve = f, year = ref$year, method = "uSUM", Dimension = b))
  res <- rbind(res, resf)
}
}
cat(unique(res$fve), "\n")
saveRDS(res, file = "./Estimate_P2_All_FVEs.rds")
}

## estimation using regression methods
dat <- readRDS("../DataResults/CleanData/Landsat_All_3Indices_5Bands_Filtered_avg.rds")
bands <- c("NDVI", "MNDWI", "B7")
res3 <- foreach(id = ids, .combine = "rbind") %dopar% {
  resi <- regthetaN(dat, id, year.range = c(2001, 2016), bands = bands)
  data.frame(pointID = id, det1 = resi[2], det2 = resi[3])
}
res3$year <- ref$year; res3$Dimension = 3
bands <- c("B2", "B3", "B4", "B5", "B7")
res5 <- foreach(id = ids, .combine = "rbind") %dopar% {
  resi <- regthetaN(dat, id, year.range = c(2001, 2016), bands = bands)
  data.frame(pointID = id, det1 = resi[2], det2 = resi[3])
}
res5$year <- ref$year; res5$Dimension = 5
saveRDS(rbind(res3, res5), file = "./Estimate_Reg_All.rds")
