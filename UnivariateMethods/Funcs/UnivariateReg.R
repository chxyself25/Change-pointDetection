## regression method implementation on univariate functional data


regthetaN <- function(dat, id, year.range = NULL, bands) {
  if (is.null(year.range)) {
    year.range <- range(dat$year[dat$pointID == id])
  }
  dati <- subset(dat, pointID == id & year <= year.range[2] & year >= year.range[1])
  values <- c(); year1s <- c(); year2s <- c()
  for (b in bands) {
    lmfit <- lm(as.formula(paste0(b, "~year")), dati)
    P1 <- which.min(lmfit$residuals)
    P2 <- which.max(lmfit$residuals)
    values <- c(values, abs(dati[[b]][P1] - dati[[b]][P2]))
    year1s <- c(year1s, dati$year[P1])
    year2s <- c(year2s, dati$year[P2])
  }
  type <- which.max(values)
  return(c(pointID = id, year1 = min(year1s[type], year2s[type]), year2 =  max(year1s[type], year2s[type]), value = max(values)))
}

URegEstAll <- function(raw_data, var.names, bands = c('ndvi', 'mndwi', 'b7'), year.range = NULL, save.results = TRUE, use.results = TRUE, dir = getwd()) {
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
  
  # aggregate data into year level
  data2 <- split(data, f = list(as.factor(data$pointID), as.factor(data$year)))
  data3 <- lapply(data2, function(x) {
    if (nrow(x) > 0) { 
    data.frame(pointID = unique(x$pointID), year = unique(x$year), b2 = mean(x$b2), b3 = mean(x$b3), b4 = mean(x$b4), b5 = mean(x$b5), b7 = mean(x$b7), ndvi = mean(x$ndvi), mndwi = mean(x$mndwi))} else {NULL}
})
  data3 <- do.call("rbind", data3)
  rownames(data3) <- NULL
  res <- foreach(id = ids, .combine = 'rbind') %dopar% {
    resi <- regthetaN(data3, id, year.range, bands)
    resi
  }
  return(res)
}
