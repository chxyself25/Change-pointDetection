sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
library(doParallel)
library(signal, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(dplyr, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(fdapace, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
library(Rcpp, lib.loc = "/home/xchang/R/x86_64-pc-linux-gnu-library/3.6/")
sourceDir("./Funcs/", trace = FALSE)
sourceDir("../fdapace/R/", trace = FALSE)
registerDoParallel(cores = 25)

Kd_dist <- readRDS("../Kd_simulation/Kd_ecdf_25.rds")
names.used = c('pID', 'doy', 'year', 'B1', 'B2', 'B3', 'B4', 'B5', 'B7')
FVEs <- c(55:99, 99.9)
methods <- c('IME', 'ISE', 'CSE')
#bands <- c("b2", "b3", "b4", "b5", "b7")
bands <- c("ndvi", "mndwi", "b7")

## read data which matches with reference data
if (!file.exists("Ankeny_validate_data.rds")) {
ref <- read.csv("Ankeney_ValSmp_interpreted.csv")
ids <- ref$pointID
dat <- NULL
for (seg in 1:14) {
  file.name <- paste0("./Ankeny", seg, ".rds")
  akn.seg <- readRDS(file.name)
  dat <- rbind(dat, subset(akn.seg, pID %in% ids))
}
saveRDS(dat, file = "Ankeny_validate_data.rds")
} else {dat <- readRDS("Ankeny_validate_data.rds"); ids <- unique(dat$pointID)}

cat("doing estimation by regression method", "\n")
regres <- URegEstAll(dat, var.names = names.used, bands, year.range = NULL, save.results = TRUE, use.results = FALSE, dir = "./Ankeny_validation/")
saveRDS(regres, file = "./Ankeny_validation/ankeny_validate_reg3.rds")

cat("doing estimation by IME method", "\n")
res.seg <- UfchangepointAll(dat, var.names = names.used, bands, Kdist =  Kd_dist, detect.method = NULL, est.method = 'IME',
                 year.range = NULL, FVE = c(55, 55), K = c(NULL, NULL), alpha = 0.05, save.results = TRUE, use.results = TRUE, dir = "./Ankeny_validation/")
res <- NULL
for (m in methods) {
  for (f in FVEs) {
    res.seg <- UfchangepointAll(dat, var.names = names.used, bands = c('ndvi', 'mndwi', 'b7'), Kdist =  Kd_dist, detect.method = NULL, est.method = m,
                 year.range = NULL, FVE = c(f, f), K = c(NULL, NULL), alpha = 0.05, save.results = FALSE, use.results = TRUE, dir = "./Ankeny_validation/")
    res.seg$method = m; res.seg$fve = f
    res <- rbind(res, res.seg)
  }
}
saveRDS(res, file = "./Ankeny_validation/ankeny_validate_fchange3.rds")


