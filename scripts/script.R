Sys.time0 <- Sys.time()

sink("script.log")

files <- c("simulation.Rmd")

f <- file("master.Rmd", open = "w")
for (i in 1:length(files)) {
    x <- readLines(files[i])
    writeLines(x, f)
    if (i < length(files)) { writeLines("\n---\n", f) }
    }
close(f)

library(knitr)
library(rmarkdown)
knit("master.Rmd", output = "index.md")
file.remove("master.Rmd")

sink()

capture.output(list(completionDateTime = Sys.time(),
                    executionTime = Sys.time() - Sys.time0,
                    sessionInfo = sessionInfo()),
               file = "script.log",
               append = TRUE)
