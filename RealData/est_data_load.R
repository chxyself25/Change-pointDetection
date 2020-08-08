## prepare script for testing of change-point
## prepare the data for estimation cross validation
ref1 <- read.csv(file = "../DataResults/LandsatRaw_ForDrZhu_merged_reference.csv")
ref1$pointID <- ref1$pointID + 120
ref2 <- read.csv(file = "../DataResults/LandsatRaw_ForDrZhu_merged0_reference.csv")
ref2$pointID <- ref2$pointID + 263
# remove ambiguous cases and unchange sites
rm.idx <- c(16,17,94,143) + 120
rm.idx <- c(rm.idx, 34+263)
# test if removing these id could improve results
ref <- subset(rbind(ref1, ref2), !is.na(year) & !(pointID %in% rm.idx) & year > 2000)
ids <- ref$pointID
splits <- c(rep(37,7), rep(38,3))

## univariate FPCA results: 3 indices or 5 bands
xi.ndvi <- readRDS("../DataResults/FPCAResults/xilam_est_All_uFPCA_fsd_ndvi.rds")
xi.mndwi <- readRDS("../DataResults/FPCAResults/xilam_est_All_uFPCA_fsd_mndwi.rds")
xi.b7 <- readRDS("../DataResults/FPCAResults/xilam_est_All_uFPCA_fsd_b7.rds")
xi.u3 <- list(ndvi = xi.ndvi, mndwi = xi.mndwi, b7 = xi.b7)
xi.b2 <- readRDS("../DataResults/FPCAResults/xilam_est_All_uFPCA_fsd_b2.rds")
xi.b3 <- readRDS("../DataResults/FPCAResults/xilam_est_All_uFPCA_fsd_b3.rds")
xi.b4 <- readRDS("../DataResults/FPCAResults/xilam_est_All_uFPCA_fsd_b4.rds")
xi.b5 <- readRDS("../DataResults/FPCAResults/xilam_est_All_uFPCA_fsd_b5.rds")
xi.u5 <- list(b2 = xi.b2, b3 = xi.b3, b4 = xi.b4, b5 = xi.b5, b7 = xi.b7)

# multivariate FPCA results
mxi.3 <- readRDS("../DataResults/FPCAResults/xilam_est_All_mFPCA_3Indices.rds")
mxi.5 <- readRDS("../DataResults/FPCAResults/xilam_est_All_mFPCA_5Bands.rds")


