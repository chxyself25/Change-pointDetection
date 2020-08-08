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

# densing data by kernel smoothing, bandwidth small
# x and y should be corresponding and sorted
LocFill <- function(x, toGrid, y = NULL) {
  locfit <- locpol(y~x, data = data.frame(x=x, y=y), deg = 1, bw = 1, kernel = EpaK, xeval = toGrid)
}
# r_rot1 <- locpol(yy~xx, data = data.frame(xx=xx, yy=yy), deg = 1, bw = hcv, kernel = gaussK, xeval = lat.eval)


# testing function
