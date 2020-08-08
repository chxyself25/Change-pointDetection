# functions for preprocessing the data
# calculate indices, apply SG-filter (do not need normalize here, since multivariate method is not used)

#################################################################################
#SWIR is helpful in distinguishing urban and bare land
# Band 1 - Blue
# Band 2 - Green
# Band 3 - Red
# Band 4 - Near Infrared
# Band 5 - Middle-infrared ??
# Band 7 - Short-wave Infrared
# NDVI = (NIR − Red)/(NIR + Red) = (4-3)/(4+3)
# MNDWI = (Green−MIR)/(Green+MIR) = (2-5)/(2+5)
#################################################################################

# var.names is the variable name corresponding to pointID, doy, year, bands
Preprocess <- function(raw_data, var.names = c('pointID', 'doy', 'year', 'B1', 'B2', 'B3', 'B4', 'B5', 'B7')) {
  # standardize variable names
  correct.nm <- c('pointID', 'doy', 'year', 'B1', 'B2', 'B3', 'B4', 'B5', 'B7')
  raw_data <- raw_data[var.names]
  names(raw_data)[var.names != correct.nm] <- correct.nm[var.names != correct.nm]
  # average duplicates
  data <- raw_data %>% 
    group_by(pointID, year, doy) %>%
    summarize(B1 = mean(B1), B2 = mean(B2), B3 = mean(B3), B4 = mean(B4), B5 = mean(B5), B7 = mean(B7)) %>% as.data.frame
  # calculate index
  data$NDVI <- with(data, (B4-B3)/(B4+B3))
  data$MNDWI <- with(data, (B2-B5)/(B2+B5))
  data <- na.omit(data)
  # outlier removal
  for (b in c(paste('B', c(1:5, 7), sep = ""), 'NDVI', 'MNDWI')) {
    data[[tolower(b)]] <- sgolayfilt(data[[b]])
  }
  return(data)
}


## test passed: test samples, DSM1 (cost 44s) 
# library(signal)
# library(dplyr)
# test_dat <- readRDS("./UnivariateMethods/DSM_test.rds")
# test_dat <- readRDS("./UnivariateMethods/DSM1.rds")
# start <- Sys.time()
# xc = Preprocess(test_dat, var.names = c('name', 'doy', 'year', 'B1', 'B2', 'B3', 'B4', 'B5', 'B7'))
# cat('cost: ', Sys.time() - start, '\n')

