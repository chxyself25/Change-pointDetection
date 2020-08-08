# function for doing cross validation and evaluation 

## function that estimate change-point for multiple sites: ensemble method with transformation
fchange_multi <- function(xi.res, ids, year.range = NULL, 
                          FVE = 85, beta = NULL, K = NULL, method = "sum", type = "mv", est.output = "single") {
  res <- NULL
  if (type == "mv") { # only one index as input
    res <- foreach(id = ids, .combine = 'rbind') %dopar% {
      resi <- Mfchange(xi.res, id, year.range, FVE, beta, K, method, est.output)
      # res <- rbind(res, data.frame(pointID = id, det = resi[1], value = resi[2], K = resi[3]))
      data.frame(pointID = id, det = resi[1], value = resi[2], K = resi[3])
    }
  } else if (type == "uv") {
    res <- foreach(id = ids, .combine = 'rbind') %dopar% {
      resi <- Ufchange(xi.res, id, year.range, FVE, beta, K, method, est.output)
      # res <- rbind(res, data.frame(pointID = id, det = resi[1], value = resi[2], K = resi[3]))
      data.frame(pointID = id, det = resi[1], value = resi[2], K = resi[3])
    } 
  } else {
    stop("Please enter 'mv' or 'uv' for the type!")
  }
  return(res)
}

## function that tests change-point for multiple sites: 
## for Kd = TRUE, the BSres is the simulated asymptotic distribution (list of ecdf functions)
## for Kd = FALSE, BSres is the bootstrap sampling distribution (matrix of sampling results)
ftest_multi <- function(xi.res, ids, year.range = NULL, FVE = 85, beta = NULL, K = NULL, 
                        method = "sum", type = "uv", BSres= list(), Kd = FALSE, alpha = 0.05) {
  if (type == "mv") {
    if (!Kd) {
      res <- foreach(id = ids, .combine = 'rbind') %dopar% {
        # process the bootstrap results and convert it to ecdf
        if (method == "max") {
          #BSresi <- t(apply(cbind(0, BSres[[id]]), 1, diff))
          BSresi <- BSres$A[[id]]
        } else {
          BSresi <- BSres$S[[id]]
        }
        nna.idx <- apply(BSresi, 2, function(x) {!all(is.na(x))})
        BSdisti <- apply(BSresi[, 1:max(which(nna.idx))], 2, ecdf)
        resi <- Mftest(xi.res, id, year.range, FVE, beta, K, method, BSdisti, Kd, alpha)
        data.frame(pointID = id, det = resi[1], K = resi[2])
      }
    } else {
      res <- foreach(id = ids, .combine = 'rbind') %dopar% {
        resi <- Mftest(xi.res, id, year.range, FVE, beta, K, method, BSres, Kd, alpha)
        data.frame(pointID = id, det = resi[1], K = resi[2])
      }
    }
  } else if (type == "uv") {
    if (!Kd) {
      res <- foreach(id = ids, .combine = 'rbind') %dopar% {
        # process the bootstrap results and convert it to ecdf
        if (method == "max") {
          BSdisti <- lapply(BSres, function(x) {
            #temp <- t(apply(cbind(0, x[[id]]), 1, diff))
            temp <- x$A[[id]]
            nna.idx <- apply(temp, 2, function(y) {!all(is.na(y))})
            apply(temp[, 1:max(which(nna.idx))], 2, ecdf)
          })
        } else {
          BSdisti <- lapply(BSres, function(x) {
            temp <- x$S[[id]]
            nna.idx <- apply(temp, 2, function(y) {!all(is.na(y))})
            temp[, 1:max(which(nna.idx))]
          })
        }
        resi <- Uftest(xi.res, id, year.range, FVE, beta, K, method, BSdisti, alpha)
        data.frame(pointID = id, det = resi[1], K = resi[2])
      }
    } else { # this option only applies on max ensemble
      if (method == "sum") {
        stop("only max ensemble can use asymptotic distribution!") 
      }
      res <- foreach(id = ids, .combine = 'rbind') %dopar% {
        resi <- UftestX(xi.res, id, year.range, FVE, beta, K, BSres, alpha)
        data.frame(pointID = id, det = resi[1], K = resi[2])
      }
    }
  } else {
    stop("Please enter 'mv' or 'uv' for the type!")
  }
  return(res)
}

## function for doing cross validation on tuning FVE, using the already fitted results
## there should be a metric_func defined beforehand, it can be accuracy within one year tolerance, mse
CVfchange <- function(ids, splits, dat, tune.vname, tune.vars, n = 200, zero.ids = NULL) {
  filter_arg <- paste0(tune.vname, " == tune.vars[j]")
  M <- length(splits)
  ## cross validation with 10 folds
  cv.res <- foreach(cv = 1:n, .combine = "rbind") %dopar% {
    # split data into 10 sets randomly
    idstemp <- ids
    sets <- list()
    if (!is.null(zero.ids)) {
      for (i in 1:(M-1)) {
        sets[[i]] <- sample(idstemp, size = splits[i])
        while(mean(sets[[i]] %in% zero.ids)==1 || mean(sets[[i]] %in% zero.ids)==0) {
          sets[[i]] <- sample(idstemp, size = splits[i])
        }
        idstemp <- setdiff(idstemp, sets[[i]])
      }
    } else {
      for (i in 1:(M-1)) {
        sets[[i]] <- sample(idstemp, size = splits[i])
        idstemp <- setdiff(idstemp, sets[[i]])
      }
    }
    sets[[M]] <- idstemp
    res <- NULL
    for (i in 1:M) {
      idstr <- setdiff(ids, sets[[i]])
      idstt <- sets[[i]]
      # tune FVE on train set
      error <- c()
      for (j in 1:length(tune.vars)) {
        argj <- paste0("pointID %in% idstr & ", filter_arg)
        resf <- filter_(dat, argj)
        error <- c(error, metric_func(resf$det, resf$year))
        # mse <- c(mse, mean((resf$det - resf$year)^2))
      }
      j <- which.min(error)
      # test set
      ttf <- filter_(dat, paste0("pointID %in% idstt & ", filter_arg))
      res <- rbind(res, ttf)
    }
    # mean((res$det - res$year)^2) # , mean(res[[tune.vname]]))
    metric_func(res$det, res$year)
    # type I and type II error
    # c(sum(res$det < 0.05 & res$year == 0)/sum(res$year == 0), sum(res$det > 0.05 & res$year == 1)/sum(res$year == 1))
  }
  return(cv.res)
}

## function for doing cross validation on tuning FVE, using the already fitted results
## estf1 is the main method, estf2 is the method to be compared with (use determined K from esf1)
## ids has to be sorted
CVfchange2 <- function(ids, splits, tune1, tune2, eval1, eval2, tune.vname, eval.vname, tune.vars, 
                       n = 200, int.req = FALSE, zero.ids = NULL) {
  filter_arg <- paste0(tune.vname, " == tune.vars[j]")
  M <- length(splits)
  ## cross validation with 10 folds
  cv.res <- foreach(cv = 1:n, .combine = 'rbind') %dopar% {
    # split data into 10 sets randomly
    idstemp <- ids
    sets <- list()
    if (!is.null(zero.ids)) {
      for (i in 1:(M-1)) {
        sets[[i]] <- sample(idstemp, size = splits[i])
        while(mean(sets[[i]] %in% zero.ids)==1 || mean(sets[[i]] %in% zero.ids)==0) {
          sets[[i]] <- sample(idstemp, size = splits[i])
        }
        idstemp <- setdiff(idstemp, sets[[i]])
      }
    } else {
      for (i in 1:(M-1)) {
        sets[[i]] <- sample(idstemp, size = splits[i])
        idstemp <- setdiff(idstemp, sets[[i]])
      }
    }
    sets[[M]] <- idstemp
    res1 <- NULL; res2 <- NULL
    for (i in 1:M) {
      idstr <- setdiff(ids, sets[[i]])
      idstt <- sets[[i]]
      # tune FVE on train set
      err1 <- c(); err2 <- c()
      for (j in 1:length(tune.vars)) {
        argj <- paste0("pointID %in% idstr & ", filter_arg)
        resf <- filter_(tune1, argj)
        err1 <- c(err1, metric_func(resf$det, resf$year))
        resf <- filter_(tune2, argj)
        err2 <- c(err2, metric_func(resf$det, resf$year))
      }
      # test on method 1
      j <- which.min(err1)
      ttf1 <- filter_(tune1, paste0("pointID %in% idstt & ", filter_arg))
      # test on method 2
      j <- which.min(err2)
      ttf2 <- filter_(tune2, paste0("pointID %in% idstt & ", filter_arg))
      res1 <- rbind(res1, ttf1)
      res2 <- rbind(res2, ttf2)
    }
    # take the average of parameters selected
    fix.para <- 0.5 * (res1[[eval.vname]][order(res1$pointID)] + res2[[eval.vname]][order(res2$pointID)])
    if (int.req) {
      fix.para <- ceiling(fix.para)
    }
    # comK <- 0.5 * (res1$K[order(res1$pointID)] + res2$K[order(res2$pointID)])
    # comK <- max(res1$K[order(res1$pointID)], res2$K[order(res2$pointID)])
    # method 1 performance
    pmc1 <- left_join(data.frame(pointID = ids, fix = fix.para), eval1, by = c('pointID', 'fix' = eval.vname))
    # method 2 performance
    pmc2 <- left_join(data.frame(pointID = ids, fix = fix.para), eval2, by = c('pointID', 'fix' = eval.vname))
    c(metric_func(pmc1$det, pmc1$year), metric_func(pmc2$det, pmc2$year))
    # list(pmc1, pmc2)
  }
  return(cv.res)
}