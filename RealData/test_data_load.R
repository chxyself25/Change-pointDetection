## prepare script for testing of change-point
## prepare locations ids for change-point detection
ref1 <- read.csv(file = "../DataResults/LandsatRaw_ForDrZhu_merged_reference.csv")
ref <- data.frame(pointID = 1:503, year = c(rep(0, 120), as.numeric(!is.na(ref1$year)), rep(1, 240)))
## remove the point not valid, ambigious
rm.idx <- c(16, 17, 94, 143) + 120
rm.idx <- c(rm.idx, c(34+263, 235+263)) # not valid for FPCA: not enough curves and change-year before 2001
ref <- subset(ref, !pointID %in% rm.idx)
ids <- ref$pointID
splits <- c(rep(50, 7), rep(49, 3))
zero.ids <- ids[ref$year == 0] # no change-point sites id

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

# bootstrap distribution
mBTS3 <- readRDS("../DataResults/ForTesting/bootstrap1000_All_mv_3Indices.rds")
mBTS5 <- readRDS("../DataResults/ForTesting/bootstrap1000_All_mv_5Bands.rds")
uBTSb2 <- readRDS("../DataResults/ForTesting/bootstrap1000_All_uv_fsd_b2_check.rds")
uBTSb3 <- readRDS("../DataResults/ForTesting/bootstrap1000_All_uv_fsd_b3.rds")
uBTSb4 <- readRDS("../DataResults/ForTesting/bootstrap1000_All_uv_fsd_b4.rds")
uBTSb5 <- readRDS("../DataResults/ForTesting/bootstrap1000_All_uv_fsd_b5.rds")
uBTSb7 <- readRDS("../DataResults/ForTesting/bootstrap1000_All_uv_fsd_b7.rds")
uBTSndvi <- readRDS("../DataResults/ForTesting/bootstrap1000_All_uv_fsd_ndvi_check.rds")
uBTSmndwi <- readRDS("../DataResults/ForTesting/bootstrap1000_All_uv_fsd_mndwi.rds")
uBTS3 <- list(uBTSndvi, uBTSmndwi, uBTSb7)
uBTS5 <- list(uBTSb2, uBTSb3, uBTSb4, uBTSb5, uBTSb7)

# xc1 <- readRDS("./DataResults/ForTesting/bootstrap1000_All_uv_fsd_mndwi_1-25.rds")
# xc2 <- readRDS("./DataResults/ForTesting/bootstrap1000_All_uv_fsd_mndwi_26-503.rds")
# xc2$S[1:25] <- xc1$S; xc2$A[1:25] <- xc1$A
# which(sapply(xc2$S, is.null))
# which(sapply(xc2$A, is.null))
# saveRDS(xc2, file = "./DataResults/ForTesting/bootstrap1000_All_uv_fsd_mndwi.rds")


