################ some help functions ############

# replace method names with the one defined in paper
change_method_names <- function(methods, original = NULL, new = NULL) {
  if (is.null(original)) {
    original = c('MMCE', 'MSCE', "UMCE", "USCE")
  }
  if (is.null(new)) {
    new <- c("mMAX", "mSUM", "uMAX", "uSUM")
  }
  if (length(original) != length(new)) {
    stop("the input arguments do not match!")
  }
  methods <- as.character(methods)
  for (i in 1:length(original)) {
    methods <- gsub(original[i], new[i], methods)
    #methods <- replace(methods, which(methods == original[i]), rep(new[i], sum(methods == original[i])))
  }
  return(methods)
}