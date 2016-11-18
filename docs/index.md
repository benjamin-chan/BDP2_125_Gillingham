---
title: "Metabolomics of very long-chain aclCoA dehydrogenase knockout mice"
date: "2016-11-18 15:27:17"
author: Benjamin Chan (chanb@ohsu.edu)
output:
  html_document:
    toc: true
    theme: simplex
---

---

# Preamble

Set working directory.


```r
setwd("~/Projects/GillinghamMetabolomics/scripts")
```

Load libraries.


```r
library(readxl)
library(magrittr)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(nlme)
```

```
## 
## Attaching package: 'nlme'
```

```
## The following object is masked from 'package:dplyr':
## 
##     collapse
```

```r
library(ggplot2)
```

Reproducibility steps.


```r
sessionInfo()
```

```
## R version 3.3.1 (2016-06-21)
## Platform: x86_64-redhat-linux-gnu (64-bit)
## Running under: CentOS release 6.8 (Final)
## 
## locale:
##  [1] LC_CTYPE=en_US.iso885915       LC_NUMERIC=C                  
##  [3] LC_TIME=en_US.iso885915        LC_COLLATE=en_US.iso885915    
##  [5] LC_MONETARY=en_US.iso885915    LC_MESSAGES=en_US.iso885915   
##  [7] LC_PAPER=en_US.iso885915       LC_NAME=C                     
##  [9] LC_ADDRESS=C                   LC_TELEPHONE=C                
## [11] LC_MEASUREMENT=en_US.iso885915 LC_IDENTIFICATION=C           
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
## [1] ggplot2_2.1.0     nlme_3.1-128      dplyr_0.5.0       magrittr_1.5     
## [5] readxl_0.1.1      rmarkdown_1.0     knitr_1.14        checkpoint_0.3.16
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.7      lattice_0.20-34  digest_0.6.10    assertthat_0.1  
##  [5] plyr_1.8.4       grid_3.3.1       R6_2.1.3         gtable_0.2.0    
##  [9] DBI_0.5-1        formatR_1.4      scales_0.4.0     evaluate_0.9    
## [13] stringi_1.1.1    tools_3.3.1      stringr_1.1.0    munsell_0.4.3   
## [17] colorspace_1.2-6 htmltools_0.3.5  methods_3.3.1    tibble_1.2
```

```r
set.seed(as.integer(as.Date("2016-11-18")))
```

Source user-defined functions.


```r
sapply(list.files("../lib", full.names = TRUE), source)
```

```
##         ../lib/library.R
## value   ?               
## visible FALSE
```

---

# Read data

Import the data.
Data files are locally stored.


```r
D1 <- importData("../data/raw/Sample stats formatting - Gillingham.xlsx")
```

Check the data.


```r
# tableChr(D1)
summary(D1$value)
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##  0.00000  0.01702  0.04656  0.57780  0.14720 22.46000
```

---

# Model

Parameters.

* `genotype` is wild type `WT` or knockout `KO`
* `activity` is `rest` or `exhaustive`
* `n` is the number of mice per (`genotype`, `activity`)
* `m` is the number of metabolites
* `j` is number of metabolite types
* `k` is number of metabolites per type
* `rho1` is the within-metabolite correlation
* `rho2` is the between-metabolite correlation

Define the fixed effects part of the model.


```r
# fixed <- formula(value ~
#                    genotype +
#                    activity +
#                    metabolite_type +
#                    genotype * metabolite_type +
#                    activity * metabolite_type)
fixed <- formula(value ~ activity)
```

Define the random effects part of the model.
Use `corCompSymm`, *compound symmetry structure corresponding to a constant correlation*.


```r
random <- formula(~ 1 | id)
cs <-
  corCompSymm(form = random, fixed = FALSE) %>%
  Initialize(data = D1)
Dim(cs)
```

```
## $N
## [1] 237
## 
## $M
## [1] 3
## 
## $maxLen
## [1] 79
## 
## $sumLenSq
## [1] 18723
## 
## $len
## groups
## 1010 1012 1134 
##   79   79   79 
## 
## $start
## [1]   0  79 158
```

```r
# summary(cs)
```

Estimate model.
Default to `correlation = NULL`, corresponding to no within-group correlations.


```r
M <- lme(fixed, data = D1, random = random, correlation = NULL)
summary(M)
```

```
## Linear mixed-effects model fit by REML
##  Data: D1 
##        AIC      BIC    logLik
##   1067.939 1081.777 -529.9693
## 
## Random effects:
##  Formula: ~1 | id
##          (Intercept) Residual
## StdDev: 0.0001363578 2.261847
## 
## Fixed effects: list(fixed) 
##                   Value Std.Error  DF   t-value p-value
## (Intercept)   0.6010433 0.1799429 234  3.340189  0.0010
## activityrest -0.0697032 0.3116703   1 -0.223644  0.8599
##  Correlation: 
##              (Intr)
## activityrest -0.577
## 
## Standardized Within-Group Residuals:
##        Min         Q1        Med         Q3        Max 
## -0.2657312 -0.2535524 -0.2346657 -0.1821647  9.6641987 
## 
## Number of Observations: 237
## Number of Groups: 3
```

Estimate model.
Specify the correlation structure.


```r
ctrl <- lmeControl(maxIter = 500)
M <- lme(fixed, data = D1, random = random, correlation = cs, control = ctrl)
summary(M)
```

```
## Linear mixed-effects model fit by REML
##  Data: D1 
##        AIC      BIC    logLik
##   1069.403 1086.701 -529.7014
## 
## Random effects:
##  Formula: ~1 | id
##         (Intercept) Residual
## StdDev:  0.06114966 2.254205
## 
## Correlation Structure: Compound symmetry
##  Formula: ~1 | id 
##  Parameter estimate(s):
##          Rho 
## -0.009858475 
## Fixed effects: list(fixed) 
##                   Value  Std.Error  DF t-value p-value
## (Intercept)   0.6010433 0.09643696 234  6.2325  0.0000
## activityrest -0.0697032 0.16703371   1 -0.4173  0.7483
##  Correlation: 
##              (Intr)
## activityrest -0.577
## 
## Standardized Within-Group Residuals:
##        Min         Q1        Med         Q3        Max 
## -0.2752325 -0.2567077 -0.2345799 -0.1801539  9.6883610 
## 
## Number of Observations: 237
## Number of Groups: 3
```
