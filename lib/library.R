importDataToList <- function (f) {
  require(readxl)
  require(magrittr)
  require(dplyr)
  D <- read_excel(f, sheet = 1, na = c("", "n/a"))
  names(D) <- D %>% names() %>% tolower() %>% gsub("\\s", "_", .)
  D <-
    D %>%
    filter(!is.na(id)) %>%
    mutate(genotype = case_when(grepl("wt", condition) ~ "WT",
                                grepl("ko", condition) ~ "KO",
                                TRUE ~ NA_character_),
           condition = case_when(grepl("rest", condition) ~ "Rest",
                                 grepl("c7", condition) ~ "C7",
                                 grepl("c8", condition) ~ "C8",
                                 grepl("ex", condition) ~ "Exercise",
                                 TRUE ~ NA_character_)) %>%
    mutate(genotype = factor(genotype, levels = c("WT", "KO")),
           condition = factor(condition,
                              levels = c("Rest", "C7", "C8", "Exercise"),
                              labels = c("Rest", "C7", "C8", "Exercise"))) %>%
    mutate(condition = factor(condition)) %>%
    mutate(metabolite = gsub("\\W", "", metabolite)) %>%
    mutate(metabolite = factor(metabolite)) %>%
    select(-c(aim, group)) %>%
    filter(!is.na(z_value))
  L <- list(file = f,
            file.size = file.size(f),
            file.mtime = file.mtime(f),
            excel_sheets = excel_sheets(f),
            dim = dim(D),
            names = names(D),
            head = head(D),
            data = D)
  L
}


summarizeOutcome <- function (D) {
  x1 <- summary(D$value)
  x2 <- summary(D$logValue)
  x3 <- summary(D$zValue)
  x4 <- summary(D$zLogValue)
  L <- data.frame(rbind(x1, x2, x3, x4))
  rownames(L) <- c("raw", "log-transform", "normalized raw", "normalized log-transform")
  L
}


estimateModel <- function (data, fixed, random) {
  require(magrittr)
  require(nlme)
  ctrl <- lmeControl(opt = "optim",
                     maxIter = 500, msMaxIter = 500,
                     tolerance = 1e-6, niterEM = 25, msMaxEval = 200, msTol = 1e-7)
  cs <-
    corCompSymm(form = random, fixed = FALSE) %>%
    Initialize(data = data)
  Dim(cs)
  lme(fixed,
      data = data,
      random = random,
      correlation = cs,
      control = ctrl)
}


adjustPvalue <- function (beta) {
  require(magrittr)
  require(dplyr)
  beta %>%
    mutate(p.adjustBH = p.adjust(p.value, method = "BH"),
           sig = p.adjustBH < 0.05)
}


runClusters <- function (df, metabolites, fixed, random, xvar, contrastValue, ctrl) {
  require(magrittr)
  require(dplyr)
  require(doParallel)
  require(data.table)
  genotypes <- df %>% pull(genotype) %>% levels()
  lookup <-
    expand.grid(metabolites, genotypes, stringsAsFactors = FALSE) %>%
    data.frame %>%
    rename(metabolite = Var1, genotype = Var2)
  n <- nrow(lookup)
  cl <- makeCluster(3)
  registerDoParallel(cl)
  L <- foreach (i = 1:n) %dopar% {
    require(magrittr)
    require(dplyr)
    require(nlme)
    dfi <-
      df %>%
      mutate(metabolite = relevel(metabolite, lookup$metabolite[i])) %>%
      mutate(genotype = relevel(genotype, lookup$genotype[i]))
    cs <- corCompSymm(form = random, fixed = FALSE) %>% Initialize(data = dfi)
    M <- dfi %>% lme(fixed, data = ., random = random, correlation = cs, control = ctrl)
    M %>%
      anova(Terms = xvar) %>%
      data.frame(contrast = contrastValue,
                 metabolite = lookup$metabolite[i],
                 genotype = lookup$genotype[i],
                 beta = M %>% fixef %>% .[names(.) == paste0(xvar, contrastValue)],
                 se = M %>% summary %>% .$tTable %>% data.frame %>%
                        select(matches("Std.Error")) %>%
                        .[row.names(.) == paste0(xvar, contrastValue), ],
                 .) %>%
      mutate(lowerCL = beta + qnorm(0.025) * se,
             upperCL = beta + qnorm(0.975) * se)
  }
  stopCluster(cl)
  rbindlist(L) %>%
    mutate(p.adjustBH = p.adjust(p.value, method = "BH"),
           sig = p.adjustBH < 0.05)
}


testContrast <- function (nlmeObj, contrast) {
  require(multcomp)
  g <- glht(nlmeObj, linfct = c(contrast)) %>% summary()
  x <- g[["test"]][["coefficients"]]
  s <- g[["test"]][["sigma"]]
  p <- g[["test"]][["pvalues"]][[1]]
  result <- data.frame(contrast = contrast,
                       coefficient = x,
                       sigma = s,
                       pvalue = p,
                       stringsAsFactors = FALSE)
  rownames(result) <- NULL
  result
}
