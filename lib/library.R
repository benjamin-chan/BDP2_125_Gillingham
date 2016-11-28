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
    mutate(genotype = factor(genotype,
                             levels = c("-/-", "+/+"),
                             labels = c("Neg/Neg", "Pos/Pos"))) %>%   # Reorder factor
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
    mutate(logValue = log10(value))  # log transform
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
