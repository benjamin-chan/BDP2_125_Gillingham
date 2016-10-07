simulateCorrData <- function (R, n, prefix = "x", seed = NULL) {
  require(magrittr)
  require(data.table)
  U <- t(chol(R))
  ncol <- dim(U)[1]
  set.seed(seed)
  Z <- matrix(rnorm(ncol * n, 0, 1), nrow = ncol, ncol = n)
  X <- t(U %*% Z)
  D <-
    data.table(X) %>%
    setnames(names(.), gsub("^V", prefix, names(.))) %>%
    round(2)
  list(nominalCorr = R,
       empericalCorr = cor(X),
       data.table = D)
}
