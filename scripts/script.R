library(checkpoint)
checkpoint("2019-04-01", use.knitr = TRUE)

files <- c("header.yaml",
           "preamble.Rmd",
           "read.Rmd",
           "model.Rmd")
f <- file("master.Rmd", open = "w")
for (i in 1:length(files)) {
    x <- readLines(sprintf("scripts/%s", files[i]))
    writeLines(x, f)
    if (i < length(files)) {writeLines("\n---\n", f)}
}
close(f)
library(knitr)
library(rmarkdown)
opts_chunk$set(fig.path = "figures/")
knit("master.Rmd", output = "docs/BDP2_125_Gillingham.md")
pandoc("docs/BDP2_125_Gillingham.md", format = "html")
file.remove("master.Rmd")
