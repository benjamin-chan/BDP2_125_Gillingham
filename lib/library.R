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


contrast <- function (fixed, df, xvar, contrastValue, refMetabolite, refGenotype) {
  require(magrittr)
  require(dplyr)
  require(nlme)
  df <-
    df %>%
    mutate(metabolite = relevel(metabolite, refMetabolite)) %>%
    mutate(genotype = relevel(genotype, refGenotype))
  cs <- corSymm(form = random, fixed = FALSE) %>% Initialize(data = df)
  M <- df %>% lme(fixed, data = ., random = random, correlation = cs, control = ctrl)
  M %>%
    anova(Terms = xvar) %>%
    data.frame(contrast = contrastValue,
               metabolite = refMetabolite,
               genotype = refGenotype,
               beta = M %>% fixef %>% .[names(.) == paste0(xvar, contrastValue)],
               .)
}


contrastGenotype <- function(refGenotype, df, xvar, contrastValue) {
  if (xvar == "activity") {
    rbind(contrast(fixed, df, xvar, contrastValue, "3-HYDROXYBUTYRIC", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "arginine", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "CITRIC", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "FUMARIC", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "glutamine", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "isoleucine", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "LACTIC", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "LCAC total", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "leucine", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "MALIC", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "MCAC Total", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "METHYLSUCCINIC", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "PYRUVIC_P2P", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "SUCCINIC-2", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "valine", refGenotype))
  } else if (xvar == "chow") {
    rbind(contrast(fixed, df, xvar, contrastValue, "3-HYDROXYBUTYRIC", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "arginine", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "CITRIC", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "FUMARIC", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "glutamine", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "isoleucine", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "LACTIC", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "LC even AC total", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "LC odd AC total", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "leucine", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "MALIC", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "MCAC total", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "METHYLSUCCINIC", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "SUCCINIC-2", refGenotype),
          contrast(fixed, df, xvar, contrastValue, "valine", refGenotype))
  }
}
