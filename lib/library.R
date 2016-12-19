importDataToList <- function (f) {
  require(readxl)
  require(magrittr)
  require(dplyr)
  D <- read_excel(f, sheet = 1) %>% .[, 2:22]
  names(D) <- D %>% names() %>% tolower() %>% gsub("\\s", "_", .)
  vars <- c("id",
            "genotype",
            "activity",
            "chow",
            "metabolite_type",
            "metabolite",
            "value",
            "logValue",
            "zValue",
            "zLogValue",
            "important")
  D <-
    D[, !is.na(names(D))] %>%
    select(which(!is.na(names(.)))) %>%
    rename(logValue = log) %>%
    rename(zValue = z_of_normalize_values) %>%
    rename(zLogValue = z_of_log_values) %>%
    mutate(important = as.logical(`important_metabolite_(1_for_important)`)) %>%
    select(one_of(vars)) %>%
    filter(!is.na(id)) %>%
    filter(important == TRUE) %>%
    mutate(genotype = factor(genotype,
                             levels = c("+/+", "-/-"),
                             labels = c("WT", "KO"))) %>%   # Reorder factor
    mutate(activity = factor(activity,
                             levels = c("rest", "rest/fasted", "ex"),
                             labels = c("Rest", "Rest", "Exercise"))) %>%   # Reorder factor
    mutate(activity = droplevels(activity)) %>%
    mutate(chow = factor(chow,
                         levels = c("reg", "White (C7)", "yellow (C8)"),
                         labels = c("Regular", "White (C7)", "Yellow (C8)"))) %>%   # Reorder factor
    mutate(chow = droplevels(chow)) %>%
    mutate(metabolite_type = factor(metabolite_type,
                                    levels = c("acylcarnitines", "amino acids", "Amino Acids", "organic acids"),
                                    labels = c("Acylcarnitines", "Amino acids", "Amino acids", "Organic acids"))) %>%   # Reorder factor
    mutate(metabolite_type = droplevels(metabolite_type)) %>% 
    mutate(metabolite = factor(metabolite))
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
  L <- data.frame(rbind(x1, x2))
  rownames(L) <- c("nominal", "log-transform")
  L
}


runClusters <- function (df, metabolites, fixed, xvar, contrastValue) {
  require(magrittr)
  require(dplyr)
  require(doParallel)
  require(data.table)
  genotypes <- c("WT", "KO")
  lookup <-
    expand.grid(metabolites, genotypes, stringsAsFactors = FALSE) %>%
    data.frame %>%
    rename(metabolite = Var1, genotype = Var2)
  n <- nrow(lookup)
  cl <- makeCluster(10)
  registerDoParallel(cl)
  L <- foreach (i = 1:n) %dopar% {
    require(magrittr)
    require(dplyr)
    require(nlme)
    dfi <-
      df %>%
      mutate(metabolite = relevel(metabolite, lookup$metabolite[i])) %>%
      mutate(genotype = relevel(genotype, lookup$genotype[i]))
    random <- formula(~ 1 | id)
    ctrl <- lmeControl(opt = "optim",
                       maxIter = 500, msMaxIter = 500)
    cs <- corSymm(form = random, fixed = FALSE) %>% Initialize(data = dfi)
    M <- dfi %>% lme(fixed, data = ., random = random, correlation = cs, control = ctrl)
    M %>%
      anova(Terms = xvar) %>%
      data.frame(contrast = contrastValue,
                 metabolite = lookup$metabolite[i],
                 genotype = lookup$genotype[i],
                 beta = M %>% fixef %>% .[names(.) == paste0(xvar, contrastValue)],
                 .)
  }
  stopCluster(cl)
  rbindlist(L)
}
