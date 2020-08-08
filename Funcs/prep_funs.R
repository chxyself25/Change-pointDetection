## some preprocessing data functions
## normalization, interpolation
############################################################################
#####################pre-process data functions#############################
############################################################################
#doy <- min(dat.25$doy):max(dat.25$doy)
#doy <- 60:330
# function for preprocessing: dat is a data frame with year, doy, 6 bands, interpolating data
prep_data <- function(dat, band, xx, threshold = 10) {
  year <- sort(unique(dat$year))
  # first interpolate on 1 to 3*366 days
  mat <- matrix(NA, ncol = length(xx), nrow = length(year))
  for (yr in year) {
    yr.dat <- subset(dat, year == yr)
    mat[which(year == yr), ] <- IntpFill(fromGrid = yr.dat$doy, toGrid = xx, y = yr.dat[[band]])
    # doyeval <- lapply(yr.dat$doy, function(x) {(max(x-10,1)):(min(x+10, max(xx)))})
    # doyeval <- sort(unique(unlist(doyeval)))
    # mat[which(year == yr), which(xx %in% doyeval)] <- fitted(locpol(as.formula(paste0(band, "~doy")), data = yr.dat, deg = 1, bw = 15, kernel = gaussK,
    #                                   xeval = doyeval))
  }
  # remove the doys which has less or equal to 10 after interpolation
  nna.col <- apply(mat, 2, function(x) {length(which(!is.na(x)))})
  doy <- xx[which(nna.col > threshold)]
  mat <- mat[, which(nna.col > threshold)]
  # delete rows having no nonmissing values
  nna.row <- apply(mat, 1, function(x) !all(is.na(x)))
  mat <- mat[nna.row, ]
  year <- year[nna.row]
  return(list(mat = mat, doy = doy, year = year))
}

# Normalizing dataset, after prep_data()
fd_normalize <- function(mat, doy, M) {
  limits <- c(min(doy), max(doy))
  h <- diff(limits) / M
  xx <- c(limits[1], 
          seq(limits[1] + h / 2, limits[2] - h / 2, length.out=M - 1), 
          limits[2])
  N <- length(xx);
  # midpoint <- seq(limits[1], limits[2], length.out=M)
  # normalized data
  res <- matrix(NA, ncol = ncol(mat), nrow = nrow(mat))
  for (i in 2:(N-1)) {
    ids <- ((doy >= xx[i-1]) & (doy < xx[i]))
    res[,ids] <- (mat[,ids]-mean(mat[,ids], na.rm = TRUE))/sd(mat[,ids], na.rm = TRUE)
  }
  # for the last bin, include the left and right end point
  ids =  ((doy >= xx[i]) &(doy <= xx[i+1]))
  res[,ids] <- (mat[,ids]-mean(mat[,ids], na.rm = TRUE))/sd(mat[,ids], na.rm = TRUE)
  return(res)
}

# functional normalize, with respect to sparse dataset: sample mean and standard deviation (same with dense)
fd_normalize_sp <- function(dat, band, M) {
  ids <- unique(dat$pointID)
  sd.band <- NULL
  for (id in ids) {
    dati <- subset(dat, pointID == id)
    udoys <- sort(unique(dati$doy))
    uyears <- unique(dati$year)
    band.mat <- matrix(NA, ncol = length(udoys), nrow = 31)
    for (yr in uyears) {
      doy.yr <- subset(dati, year == yr)
      band.mat[which(yr == uyears), which(udoys %in% doy.yr$doy)] <- doy.yr[[band]]
    }
    sd.band.mat <- c(t(fd_normalize(band.mat, udoys, M = M)))
    sd.band <- c(sd.band, na.omit(sd.band.mat))
  }
  return(sd.band)
}

## sparse functional normalize: smoothing method to get mean function and variance
fd_normalize_sparse <- function(dat, band) {
  ids <- unique(dat$pointID)
  res <- foreach(id = ids, .combine = 'c') %dopar% {
    dati <- subset(dat, pointID == id)
    uyears <- unique(dati$year)
    Ly <- lapply(split(dati, f = as.factor(dati$year)), function(x) {x[[band]]})
    Lt <- lapply(split(dati, f = as.factor(dati$year)), function(x) {x$doy})
    band.pca <- FPCA(Ly, Lt, optns = list(dataType = "Sparse"))
    norm.band <- c()
    for (yr in uyears) {
      doy.yr <- subset(dati, year == yr)
      sd.band <- IntpFill(fromGrid = band.pca$workGrid, toGrid = doy.yr$doy, 
                          y = diag(band.pca$fittedCov) + band.pca$sigma2)
      mu.band <- IntpFill(fromGrid = band.pca$workGrid, toGrid = doy.yr$doy, 
                          y = band.pca$mu)
      norm.band <- c(norm.band, (doy.yr[[band]] - mu.band)/sqrt(sd.band))
    }
    norm.band
  }
  return(res)
}

#' Convert support of a mu/phi/cov etc. to and from obsGrid and workGrid
#' 
#' Convert the support of a given function 1-D or 2-D function from 'fromGrd' to 'toGrid'.
#' Both grids need to be sorted. This is a interpolation/convenience function.


IntpFill <- function(fromGrid, toGrid, y=NULL, Cov=NULL) {
  # sort the grids in case they are ordered
  y <- y[order(fromGrid)]
  fromGrid <- sort(fromGrid)
  # In case the range of toGrid is larger than fromGrid due to numeric error
  buff <- .Machine$double.eps * max(abs(fromGrid)) * 3
  if (abs(toGrid[1] - fromGrid[1]) < buff)
    toGrid[1] <- fromGrid[1]
  if (abs(toGrid[length(toGrid)] - fromGrid[length(fromGrid)]) < buff)
    toGrid[length(toGrid)] <- fromGrid[length(fromGrid)]
  # if ( ( fromGrid[1] - buff  >  toGrid[1]) || 
  #      ( fromGrid[length(fromGrid)] + buff < toGrid[length(toGrid)]) ) {
  #   stop("Insufficient size of 'fromGrid'.")}
  
  if (is.null(y)) {
    gd <- expand.grid(X=toGrid, Y=toGrid)
    res <- matrix(interp2lin(fromGrid, fromGrid, Cov, gd$X, gd$Y), nrow=length(toGrid))
  } else {
    togrid <- toGrid
    if (fromGrid[1] >  toGrid[1]) {
      togrid[toGrid < fromGrid[1]] <- NA
    }
    if (fromGrid[length(fromGrid)] < toGrid[length(toGrid)]) {
      togrid[toGrid > fromGrid[length(fromGrid)]] <- NA 
    }
    res <- rep(NA, length(togrid))
    if (length(fromGrid) == 1) {
      res[which.min(abs(togrid-fromGrid))] <- y
    }else {
      res[!is.na(togrid)] <- approx(fromGrid, y, togrid[!is.na(togrid)], method='linear')$y 
    }
    return(res)
  }
}
