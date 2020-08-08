source("../Funcs/urban_detect_correct.R")
source("../Funcs/eval_func.R")
source("./test_data_load.R")
library(doParallel)
registerDoParallel(cores = 50)
Kdist <- readRDS("../Kd_simulation/Kd_ecdf_25.rds")
FVEs <- c(26:99, 99.5, 99.9)
#FVEs <- c(26:45)

## testing using asymptotic distrbution
if (FALSE) {
res <- readRDS("./Detect_A_All_FVEs.rds")
for (f in FVEs) {
  cat("doing ", f, "\n")
  for (b in c(3,5)) {
  resf <- ftest_multi(get(paste0("mxi.", b)), ids = ids, year.range = c(2001, 2016), FVE = f, beta = NULL, K = NULL,
                      method = "max", type = "mv", BSres = Kdist, Kd = TRUE)
  resf <- cbind(resf, data.frame(fve = f, year = ref$year, method = "mMAX", Dimensions = b, dist = "Asymptotic"))
  res <- rbind(res, resf)
  
  resf <- ftest_multi(get(paste0("mxi.", b)), ids = ids, year.range = c(2001, 2016), FVE = f, beta = NULL, K = NULL,
                      method = "sum", type = "mv", BSres = Kdist, Kd = TRUE)
  resf <- cbind(resf, data.frame(fve = f, year = ref$year, method = "mSUM", Dimensions = b, dist = "Asymptotic"))
  res <- rbind(res, resf)
  
  resf <- ftest_multi(get(paste0("xi.u", b)), ids = ids, year.range = c(2001, 2016), FVE = f, beta = NULL, K = NULL,
                      method = "max", type = "uv", BSres = Kdist, Kd = TRUE)
  resf <- cbind(resf, data.frame(fve = f, year = ref$year, method = "uMAX", Dimensions = b, dist = "Asymptotic"))
  res <- rbind(res, resf)
}
}
cat(unique(res$fve), "\n")
saveRDS(res, file = "./Detect_A_All_FVEs.rds")
}

## testing using bootstrap distribution
if (TRUE) {
#res <- readRDS("./Detect_BTS_All_FVEs.rds")
res <- NULL
for (f in FVEs) {
  cat("doing ", f, "\n")
  for (b in c(3,5)) {
  resf <- ftest_multi(get(paste0("mxi.", b)), ids = ids, year.range = c(2001, 2016), FVE = f, beta = NULL, K = NULL,
                      method = "max", type = "mv", BSres = get(paste0("mBTS", b)), Kd = FALSE)
  resf <- cbind(resf, data.frame(fve = f, year = ref$year, method = "mMAX", Dimensions = b, dist = "Bootstrap"))
  res <- rbind(res, resf)
  
  resf <- ftest_multi(get(paste0("mxi.", b)), ids = ids, year.range = c(2001, 2016), FVE = f, beta = NULL, K = NULL,
                      method = "sum", type = "mv", BSres = get(paste0("mBTS", b)), Kd = FALSE)
  resf <- cbind(resf, data.frame(fve = f, year = ref$year, method = "mSUM", Dimensions = b, dist = "Bootstrap"))
  res <- rbind(res, resf)

  resf <- ftest_multi(get(paste0("xi.u", b)), ids = ids, year.range = c(2001, 2016), FVE = f, beta = NULL, K = NULL,
                      method = "max", type = "uv", BSres = get(paste0("uBTS", b)), Kd = FALSE)
  resf <- cbind(resf, data.frame(fve = f, year = ref$year, method = "uMAX", Dimensions = b, dist = "Bootstrap"))
  res <- rbind(res, resf)
  
  resf <- ftest_multi(get(paste0("xi.u", b)), ids = ids, year.range = c(2001, 2016), FVE = f, beta = NULL, K = NULL,
                      method = "sum", type = "uv", BSres = get(paste0("uBTS", b)), Kd = FALSE)
  resf <- cbind(resf, data.frame(fve = f, year = ref$year, method = "uSUM", Dimensions = b, dist = "Bootstrap"))
  res <- rbind(res, resf)
}
}
cat(unique(res$fve), "\n")
saveRDS(res, file = "./Detect_BTS_All_FVEs_check.rds")
}

## testing using regression method
if (FALSE) {
dat <- readRDS("../DataResults/CleanData/Landsat_All_3Indices_5Bands_Filtered_avg.rds")
bands <- c("NDVI", "MNDWI", "B7")
res3 <- foreach(id = ids, .combine = "rbind") %dopar% {
  resi <- regTest(dat, id, year.range = c(2001, 2016), bands = bands)
  data.frame(pointID = id, det = resi[2])
}
res3$year <- ref$year; res3$Dimension = 3
bands <- c("B2", "B3", "B4", "B5", "B7")
res5 <- foreach(id = ids, .combine = "rbind") %dopar% {
  resi <- regTest(dat, id, year.range = c(2001, 2016), bands = bands)
  data.frame(pointID = id, det = resi[2])
}
res5$year <- ref$year; res5$Dimension = 5
saveRDS(rbind(res3, res5), file = "./Detect_Reg_All.rds")
}


