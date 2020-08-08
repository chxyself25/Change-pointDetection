# main function: conduct change-point detection and estimation by univariate methods
# save intermediate results if necessary
# require: fdapace, signal, dplyr

# function for change-point detection and estimation for one pixel data
# dat is the preprocessed data from one pixel
Ufchangepoint <- function(dat, bands = c('ndvi', 'mndwi', 'b7'), Kdist = list(), detect.method = 'ISE', est.method = 'CSE', 
                          year.range = NULL, FVE = c(85, 85), K = c(NULL, NULL), alpha = 0.05) {
  # Univariate FPCA
  pca <- GetFPCA(dat, bands)
  # change-point detection
  detect.res <- UftestA(pca, year.range, FVE = FVE[1], K = K[1], Kdist, alpha, method = detect.method)
  if (detect.res[1]) {
    # change-point estimation
    est.res <- Ufchange(pca, year.range, FVE = FVE[2], K = K[2], method = est.method)
  } else {
    est.res <- NULL
  }
  
  return(list(detect = detect.res, est = est.res))
}


# wrap up function for change-point detection and estimation for one dataset (multiple pixels)
UfchangepointAll <- function(raw_data, var.names, bands = c('ndvi', 'mndwi', 'b7'), Kdist = list(), detect.method = 'ISE', est.method = 'CSE', 
                             year.range = NULL, FVE = c(85, 85), K = c(NULL, NULL), alpha = 0.05,
                             save.results = TRUE, use.results = TRUE, dir = getwd()) {
  # preprocess raw data
  cat("Organize the raw dataset..", '\n')
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  data.file <- paste0(dir, "/raw_data_preprocess.rds")
  if (file.exists(data.file) & use.results) {
    data <- readRDS(data.file)
  } else {
    data <- Preprocess(raw_data, var.names)
    if (save.results) {
      saveRDS(data, data.file)
    }
  }
  
  ids <- unique(data$pointID)
  # univariate FPCA
  cat("Doing univariate FPCA...", '\n')
  pca.file <- paste0(dir, "/univariate_fpca_", paste(bands, collapse = "_"), ".rds")
  if (file.exists(pca.file) & use.results) {
    pca <- readRDS(pca.file)
  } else {
    pca = foreach(i = ids) %dopar% {
      datai <- subset(data, pointID == i)
      #cat("Doing FPCA for pixel ", i, "\n")
      pcai <- GetFPCA(datai, bands)
      pcai
    }
    #pca <- list()
    #for (i in 1:length(ids)) {
    #  cat("Doing FPCA for ", i, "\n")
    #  datai <- subset(data, pointID == ids[i])
    #  pca[[i]] <- GetFPCA(datai, bands)
    #}
    if (save.results) {
      saveRDS(pca, pca.file)
    }
  }
  rm(data)

  nna.idx <- which(sapply(pca, function(x) {length(x) == 3}))
  # change-point detection and estimation
  cat("Change-point detection and estimation for each pixel...", '\n')
  res <- foreach(i = nna.idx, .combine = 'rbind') %dopar% {
    pcai <- pca[[i]]
    # change-point detection
    if (!is.null(detect.method)) {
      detect.res <- UftestA(pcai, year.range, FVE = FVE[1], K = K[1], Kdist, alpha, method = detect.method)
    }  else {detect.res = rep(NA, 2)}
    if (!is.null(est.method)) {
     # change-point estimation
     est.res <- Ufchange(pcai, year.range, FVE = FVE[2], K = K[2], method = est.method)
    } else {
     est.res = rep(NA, 4)
    }
    #est.res <- Ufchange(pcai, year.range, FVE = FVE[2], K = K[2], method = est.method)
    data.frame(pointID = ids[i], detect = detect.res[1], detectK = detect.res[2], 
               P1 = est.res[1], P2 = est.res[2], estTN = est.res[3], estK = est.res[4])
    #data.frame(pointID = ids[i], P1 = est.res[1], P2 = est.res[2], estTN = est.res[3], estK = est.res[4])
  }
  rm(pca)
  
  return(res)
}

## test passed: test sample, pixel 1 and 100
# test_dat <- readRDS("./UnivariateMethods/DSM_test.rds")
# Kd_dist <- readRDS("./DataResults/ForTesting/Kd_ecdf_25.rds")
# data = Preprocess(test_dat, var.names = c('name', 'doy', 'year', 'B1', 'B2', 'B3', 'B4', 'B5', 'B7'))
# datai <- subset(data, pointID == 100)
# xc <- Ufchangepoint(datai, bands = c('ndvi', 'mndwi', 'b7'), Kdist = Kd_dist, detect.method = 'MIE', 
#                     est.method = 'SCE', year.range = NULL, FVE = c(99, 85))



