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


importDataToList <- function (f) {
  require(readxl)
  require(magrittr)
  require(dplyr)
  col_types <- c(rep("text", 6), "numeric")
  D <- read_excel(f, col_types = col_types)
  names(D) <- D %>% names() %>% tolower() %>% gsub("\\s", "_", .)
  D <-
    D %>%
    filter(!is.na(id)) %>%
    mutate(activity = factor(activity, levels = c("rest", "rest/fasted", "ex"))) %>%   # Reorder factor
    mutate(logValue = log10(value))  # log transform
  L <- list(file = f,
            file.size = file.size(f),
            file.mtime = file.mtime(f),
            excel_sheets = excel_sheets(f),
            nrow = nrow(D),
            ncol = ncol(D),
            names = names(D),
            head = head(D),
            data = D)
  L
}


tableChr <- function (D) {
  table(D$id)
  isChr <- sapply(D, class) == "character"
  for (i in 1:sum(isChr)) {
    j <- which(isChr)[i]
    L <- list(variable = names(j),
              table = table(D[, j]))
    show(L)
  }
}
