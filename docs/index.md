---
title: "Metabolomics of very long-chain aclCoA dehydrogenase knockout mice"
date: "2016-11-23 14:22:29"
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
library(knitr)
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

Import Aim 1: Exercise data


```r
L1 <- importDataToList("../data/raw/Ultragenyx Aim 1 combined.xlsx")
```

```
## Warning in `levels<-`(`*tmp*`, value = if (nl == nL) as.character(labels)
## else paste0(labels, : duplicated levels in factors are deprecated

## Warning in `levels<-`(`*tmp*`, value = if (nl == nL) as.character(labels)
## else paste0(labels, : duplicated levels in factors are deprecated
```

```r
L1[["file"]]
```

```
## [1] "../data/raw/Ultragenyx Aim 1 combined.xlsx"
```

```r
L1[["dim"]]
```

```
## [1] 6854    8
```

```r
L1[["names"]]
```

```
## [1] "id"              "genotype"        "activity"        "chow"           
## [5] "metabolite_type" "metabolite"      "value"           "logValue"
```

```r
L1[["head"]]
```

```
## # A tibble: 6 × 8
##       id genotype activity    chow metabolite_type
##    <chr>   <fctr>   <fctr>  <fctr>           <chr>
## 1 1208B   Pos/Pos Exercise Regular  acylcarnitines
## 2 1208B   Pos/Pos Exercise Regular  acylcarnitines
## 3 1208B   Pos/Pos Exercise Regular  acylcarnitines
## 4 1208B   Pos/Pos Exercise Regular  acylcarnitines
## 5 1208B   Pos/Pos Exercise Regular  acylcarnitines
## 6 1208B   Pos/Pos Exercise Regular  acylcarnitines
## # ... with 3 more variables: metabolite <chr>, value <dbl>, logValue <dbl>
```

```r
D1 <- L1[["data"]]
```


```r
L2 <- importDataToList("../data/raw/Ultragenyx Aim 2 combined.xlsx")
```

```
## Warning in `levels<-`(`*tmp*`, value = if (nl == nL) as.character(labels)
## else paste0(labels, : duplicated levels in factors are deprecated

## Warning in `levels<-`(`*tmp*`, value = if (nl == nL) as.character(labels)
## else paste0(labels, : duplicated levels in factors are deprecated
```

```r
L2[["file"]]
```

```
## [1] "../data/raw/Ultragenyx Aim 2 combined.xlsx"
```

```r
L2[["dim"]]
```

```
## [1] 7436    8
```

```r
L2[["names"]]
```

```
## [1] "id"              "genotype"        "activity"        "chow"           
## [5] "metabolite_type" "metabolite"      "value"           "logValue"
```

```r
L2[["head"]]
```

```
## # A tibble: 6 × 8
##      id genotype activity       chow metabolite_type
##   <chr>   <fctr>   <fctr>     <fctr>           <chr>
## 1  1120  Neg/Neg Exercise White (C7)  acylcarnitines
## 2  1120  Neg/Neg Exercise White (C7)  acylcarnitines
## 3  1120  Neg/Neg Exercise White (C7)  acylcarnitines
## 4  1120  Neg/Neg Exercise White (C7)  acylcarnitines
## 5  1120  Neg/Neg Exercise White (C7)  acylcarnitines
## 6  1120  Neg/Neg Exercise White (C7)  acylcarnitines
## # ... with 3 more variables: metabolite <chr>, value <dbl>, logValue <dbl>
```

```r
D2 <- L2[["data"]]
```


## Check data

Check the `value` and `logValue`.


```r
kable(summarizeOutcome(D1))
```



|              |       Min.| X1st.Qu.|  Median|    Mean| X3rd.Qu.|      Max.|
|:-------------|----------:|--------:|-------:|-------:|--------:|---------:|
|nominal       |  0.0000001|  0.03156|  0.1242| 87.5100|   1.3430| 81090.000|
|log-transform | -6.9550000| -1.50100| -0.9060| -0.7525|   0.1281|     4.909|

```r
ggplot(D1) +
  aes(x = value) +
  geom_density(alpha = 1/2, fill = "blue") +
  scale_x_log10("log10 scale") +
  facet_grid(genotype ~ activity) +
  theme_bw()
```

![plot of chunk densitiesAim1](../figures/densitiesAim1-1.png)


```r
kable(summarizeOutcome(D2))
```



|              |       Min.| X1st.Qu.|  Median|     Mean| X3rd.Qu.|      Max.|
|:-------------|----------:|--------:|-------:|--------:|--------:|---------:|
|nominal       |  0.0000001|   0.0309|  0.1291| 276.5000|   1.2660| 4.949e+05|
|log-transform | -6.9550000|  -1.5100| -0.8891|  -0.7612|   0.1025| 5.695e+00|

```r
ggplot(D2) +
  aes(x = value) +
  geom_density(alpha = 1/2, fill = "blue") +
  scale_x_log10("log10 scale") +
  facet_grid(genotype ~ chow) +
  theme_bw()
```

![plot of chunk densitiesAim2](../figures/densitiesAim2-1.png)


Check fixed effects factors.


```r
tableFixed(D1)
```

```
## Error in tableFixed(D1): could not find function "show"
```

```r
tableFixed(D2)
```

```
## Error in tableFixed(D2): could not find function "show"
```

---

# Model


## Aim 1: Exercise

**Rest versus exhaustive exercise**

Estimate model.
Default to `correlation = NULL`, corresponding to no within-group correlations.


```r
fixed <- formula(logValue ~
                   genotype +
                   activity +
                   metabolite_type +
                   genotype * activity +
                   genotype * metabolite_type)
random <- formula(~ 1 | id)
M <- lme(fixed, data = D1, random = random, correlation = NULL)
summary(M)
```

```
## Linear mixed-effects model fit by REML
##  Data: D1 
##        AIC      BIC    logLik
##   23849.06 23917.37 -11914.53
## 
## Random effects:
##  Formula: ~1 | id
##          (Intercept) Residual
## StdDev: 5.140006e-05 1.373658
## 
## Fixed effects: list(fixed) 
##                                                   Value  Std.Error   DF
## (Intercept)                                  -1.4713162 0.07331735 6801
## genotypePos/Pos                               0.0742170 0.10414511   45
## activityExercise                             -0.0062902 0.06510023 6801
## metabolite_typeamino acids                    1.5020471 0.06213088 6801
## metabolite_typeorganic acids                  1.1050970 0.06885615 6801
## genotypePos/Pos:activityExercise              0.0134917 0.09258194 6801
## genotypePos/Pos:metabolite_typeamino acids   -0.0853005 0.08798927 6801
## genotypePos/Pos:metabolite_typeorganic acids  0.0066969 0.09762673 6801
##                                                 t-value p-value
## (Intercept)                                  -20.067775  0.0000
## genotypePos/Pos                                0.712631  0.4798
## activityExercise                              -0.096623  0.9230
## metabolite_typeamino acids                    24.175532  0.0000
## metabolite_typeorganic acids                  16.049358  0.0000
## genotypePos/Pos:activityExercise               0.145727  0.8841
## genotypePos/Pos:metabolite_typeamino acids    -0.969443  0.3324
## genotypePos/Pos:metabolite_typeorganic acids   0.068597  0.9453
##  Correlation: 
##                                              (Intr) gntP/P actvtE
## genotypePos/Pos                              -0.704              
## activityExercise                             -0.888  0.625       
## metabolite_typeamino acids                   -0.693  0.488  0.499
## metabolite_typeorganic acids                 -0.625  0.440  0.450
## genotypePos/Pos:activityExercise              0.624 -0.889 -0.703
## genotypePos/Pos:metabolite_typeamino acids    0.489 -0.694 -0.352
## genotypePos/Pos:metabolite_typeorganic acids  0.441 -0.614 -0.318
##                                              mtblt_typma mtblt_typra
## genotypePos/Pos                                                     
## activityExercise                                                    
## metabolite_typeamino acids                                          
## metabolite_typeorganic acids                  0.490                 
## genotypePos/Pos:activityExercise             -0.351      -0.317     
## genotypePos/Pos:metabolite_typeamino acids   -0.706      -0.346     
## genotypePos/Pos:metabolite_typeorganic acids -0.346      -0.705     
##                                              gP/P:E gntypPs/Ps:mtblt_typma
## genotypePos/Pos                                                           
## activityExercise                                                          
## metabolite_typeamino acids                                                
## metabolite_typeorganic acids                                              
## genotypePos/Pos:activityExercise                                          
## genotypePos/Pos:metabolite_typeamino acids    0.501                       
## genotypePos/Pos:metabolite_typeorganic acids  0.439  0.485                
## 
## Standardized Within-Group Residuals:
##        Min         Q1        Med         Q3        Max 
## -5.0852593 -0.5139823  0.0642195  0.4939858  3.8402323 
## 
## Number of Observations: 6854
## Number of Groups: 47
```

Estimate model.
Specify the correlation structure using `cs`.
Use `corSymm`, *general correlation matrix, with no additional structure*.


```r
cs <-
  corSymm(form = random, fixed = TRUE) %>%
  Initialize(data = D1)
Dim(cs)
```

```
## $N
## [1] 6854
## 
## $M
## [1] 47
## 
## $maxLen
## [1] 164
## 
## $sumLenSq
## [1] 1052486
## 
## $len
## groups
## 1208B   1208B 1207B   1207B 1206B   1206B 1205B   1205B 1204B   1204B 
##     79     85     79     85     79     85     79     85     79     85 
##   1212   1211   1210   1209   1202   1201   1172   1171   1170   1164 
##    164    164    164    164    130    164    164    164    164    164 
##   1163   1151   1150   1135   1134   1095   1094   1092   1091   1090 
##    164    164    164    164    164    164    164    164    164    164 
##   1077   1076   1073   1066   1060   1046   1034   1033   1030   1029 
##    164    164    164    164    164    164    164    164    164    164 
##   1028   1019   1018   1017   1014   1012   1010 
##    164    164    164    164    164    164    164 
## 
## $start
##  [1]    0   79  164  243  328  407  492  571  656  735  820  984 1148 1312
## [15] 1476 1606 1770 1934 2098 2262 2426 2590 2754 2918 3082 3246 3410 3574
## [29] 3738 3902 4066 4230 4394 4558 4722 4886 5050 5214 5378 5542 5706 5870
## [43] 6034 6198 6362 6526 6690
```

```r
M <- lme(fixed, data = D1, random = random, correlation = cs)
summary(M)
```

```
## Linear mixed-effects model fit by REML
##  Data: D1 
##        AIC      BIC    logLik
##   23849.06 23917.37 -11914.53
## 
## Random effects:
##  Formula: ~1 | id
##          (Intercept) Residual
## StdDev: 5.140488e-05 1.373658
## 
## Correlation Structure: General
##  Formula: ~1 | id 
##  Parameter estimate(s):
##  Correlation: 
##     1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
## 2   0                                                                   
## 3   0 0                                                                 
## 4   0 0 0                                                               
## 5   0 0 0 0                                                             
## 6   0 0 0 0 0                                                           
## 7   0 0 0 0 0 0                                                         
## 8   0 0 0 0 0 0 0                                                       
## 9   0 0 0 0 0 0 0 0                                                     
## 10  0 0 0 0 0 0 0 0 0                                                   
## 11  0 0 0 0 0 0 0 0 0 0                                                 
## 12  0 0 0 0 0 0 0 0 0 0  0                                              
## 13  0 0 0 0 0 0 0 0 0 0  0  0                                           
## 14  0 0 0 0 0 0 0 0 0 0  0  0  0                                        
## 15  0 0 0 0 0 0 0 0 0 0  0  0  0  0                                     
## 16  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0                                  
## 17  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0                               
## 18  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0                            
## 19  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0                         
## 20  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0                      
## 21  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0                   
## 22  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0                
## 23  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0             
## 24  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0          
## 25  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0       
## 26  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
## 27  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 28  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 29  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 30  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 31  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 32  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 33  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 34  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 35  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 36  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 37  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 38  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 39  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 40  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 41  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 42  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 43  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 44  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 45  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 46  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 47  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 48  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 49  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 50  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 51  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 52  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 53  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 54  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 55  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 56  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 57  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 58  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 59  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 60  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 61  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 62  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 63  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 64  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 65  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 66  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 67  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 68  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 69  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 70  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 71  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 72  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 73  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 74  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 75  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 76  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 77  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 78  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 79  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 80  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 81  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 82  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 83  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 84  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 85  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 86  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 87  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 88  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 89  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 90  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 91  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 92  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 93  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 94  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 95  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 96  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 97  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 98  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 99  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 100 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 101 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 102 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 103 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 104 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 105 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 106 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 107 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 108 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 109 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 110 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 111 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 112 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 113 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 114 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 115 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 116 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 117 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 118 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 119 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 120 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 121 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 122 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 123 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 124 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 125 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 126 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 127 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 128 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 129 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 130 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 131 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 132 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 133 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 134 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 135 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 136 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 137 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 138 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 139 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 140 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 141 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 142 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 143 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 144 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 145 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 146 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 147 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 148 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 149 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 150 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 151 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 152 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 153 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 154 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 155 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 156 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 157 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 158 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 159 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 160 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 161 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 162 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 163 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 164 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
##     27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49
## 2                                                                       
## 3                                                                       
## 4                                                                       
## 5                                                                       
## 6                                                                       
## 7                                                                       
## 8                                                                       
## 9                                                                       
## 10                                                                      
## 11                                                                      
## 12                                                                      
## 13                                                                      
## 14                                                                      
## 15                                                                      
## 16                                                                      
## 17                                                                      
## 18                                                                      
## 19                                                                      
## 20                                                                      
## 21                                                                      
## 22                                                                      
## 23                                                                      
## 24                                                                      
## 25                                                                      
## 26                                                                      
## 27                                                                      
## 28  0                                                                   
## 29  0  0                                                                
## 30  0  0  0                                                             
## 31  0  0  0  0                                                          
## 32  0  0  0  0  0                                                       
## 33  0  0  0  0  0  0                                                    
## 34  0  0  0  0  0  0  0                                                 
## 35  0  0  0  0  0  0  0  0                                              
## 36  0  0  0  0  0  0  0  0  0                                           
## 37  0  0  0  0  0  0  0  0  0  0                                        
## 38  0  0  0  0  0  0  0  0  0  0  0                                     
## 39  0  0  0  0  0  0  0  0  0  0  0  0                                  
## 40  0  0  0  0  0  0  0  0  0  0  0  0  0                               
## 41  0  0  0  0  0  0  0  0  0  0  0  0  0  0                            
## 42  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                         
## 43  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                      
## 44  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                   
## 45  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                
## 46  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0             
## 47  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0          
## 48  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0       
## 49  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
## 50  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 51  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 52  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 53  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 54  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 55  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 56  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 57  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 58  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 59  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 60  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 61  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 63  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 64  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 65  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 66  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 67  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 68  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 69  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 70  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 71  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 72  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 73  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 74  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 75  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 76  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 77  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 78  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 79  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 80  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 81  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 82  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 83  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 84  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 85  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 86  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 87  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 88  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 89  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 90  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 91  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 92  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 93  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 94  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 95  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 96  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 97  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 98  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 99  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 100 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 101 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 102 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 103 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 104 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 105 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 106 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 107 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 108 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 109 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 110 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 111 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 112 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 113 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 114 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 115 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 116 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 117 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 118 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 119 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 120 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 121 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 122 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 123 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 124 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 125 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 126 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 127 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 128 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 129 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 130 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 131 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 132 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 133 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 134 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 135 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 136 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 137 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 138 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 139 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 140 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 141 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 142 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 143 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 144 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 145 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 146 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 147 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 148 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 149 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 150 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 151 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 152 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 153 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 154 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 155 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 156 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 157 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 158 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 159 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 160 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 161 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 162 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 163 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 164 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
##     50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72
## 2                                                                       
## 3                                                                       
## 4                                                                       
## 5                                                                       
## 6                                                                       
## 7                                                                       
## 8                                                                       
## 9                                                                       
## 10                                                                      
## 11                                                                      
## 12                                                                      
## 13                                                                      
## 14                                                                      
## 15                                                                      
## 16                                                                      
## 17                                                                      
## 18                                                                      
## 19                                                                      
## 20                                                                      
## 21                                                                      
## 22                                                                      
## 23                                                                      
## 24                                                                      
## 25                                                                      
## 26                                                                      
## 27                                                                      
## 28                                                                      
## 29                                                                      
## 30                                                                      
## 31                                                                      
## 32                                                                      
## 33                                                                      
## 34                                                                      
## 35                                                                      
## 36                                                                      
## 37                                                                      
## 38                                                                      
## 39                                                                      
## 40                                                                      
## 41                                                                      
## 42                                                                      
## 43                                                                      
## 44                                                                      
## 45                                                                      
## 46                                                                      
## 47                                                                      
## 48                                                                      
## 49                                                                      
## 50                                                                      
## 51  0                                                                   
## 52  0  0                                                                
## 53  0  0  0                                                             
## 54  0  0  0  0                                                          
## 55  0  0  0  0  0                                                       
## 56  0  0  0  0  0  0                                                    
## 57  0  0  0  0  0  0  0                                                 
## 58  0  0  0  0  0  0  0  0                                              
## 59  0  0  0  0  0  0  0  0  0                                           
## 60  0  0  0  0  0  0  0  0  0  0                                        
## 61  0  0  0  0  0  0  0  0  0  0  0                                     
## 62  0  0  0  0  0  0  0  0  0  0  0  0                                  
## 63  0  0  0  0  0  0  0  0  0  0  0  0  0                               
## 64  0  0  0  0  0  0  0  0  0  0  0  0  0  0                            
## 65  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                         
## 66  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                      
## 67  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                   
## 68  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                
## 69  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0             
## 70  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0          
## 71  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0       
## 72  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
## 73  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 74  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 75  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 76  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 77  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 78  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 79  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 80  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 81  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 82  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 83  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 84  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 85  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 86  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 87  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 88  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 89  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 90  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 91  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 92  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 93  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 94  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 95  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 96  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 97  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 98  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 99  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 100 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 101 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 102 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 103 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 104 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 105 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 106 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 107 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 108 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 109 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 110 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 111 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 112 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 113 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 114 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 115 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 116 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 117 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 118 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 119 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 120 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 121 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 122 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 123 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 124 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 125 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 126 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 127 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 128 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 129 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 130 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 131 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 132 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 133 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 134 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 135 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 136 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 137 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 138 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 139 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 140 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 141 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 142 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 143 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 144 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 145 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 146 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 147 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 148 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 149 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 150 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 151 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 152 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 153 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 154 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 155 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 156 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 157 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 158 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 159 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 160 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 161 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 162 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 163 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 164 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
##     73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95
## 2                                                                       
## 3                                                                       
## 4                                                                       
## 5                                                                       
## 6                                                                       
## 7                                                                       
## 8                                                                       
## 9                                                                       
## 10                                                                      
## 11                                                                      
## 12                                                                      
## 13                                                                      
## 14                                                                      
## 15                                                                      
## 16                                                                      
## 17                                                                      
## 18                                                                      
## 19                                                                      
## 20                                                                      
## 21                                                                      
## 22                                                                      
## 23                                                                      
## 24                                                                      
## 25                                                                      
## 26                                                                      
## 27                                                                      
## 28                                                                      
## 29                                                                      
## 30                                                                      
## 31                                                                      
## 32                                                                      
## 33                                                                      
## 34                                                                      
## 35                                                                      
## 36                                                                      
## 37                                                                      
## 38                                                                      
## 39                                                                      
## 40                                                                      
## 41                                                                      
## 42                                                                      
## 43                                                                      
## 44                                                                      
## 45                                                                      
## 46                                                                      
## 47                                                                      
## 48                                                                      
## 49                                                                      
## 50                                                                      
## 51                                                                      
## 52                                                                      
## 53                                                                      
## 54                                                                      
## 55                                                                      
## 56                                                                      
## 57                                                                      
## 58                                                                      
## 59                                                                      
## 60                                                                      
## 61                                                                      
## 62                                                                      
## 63                                                                      
## 64                                                                      
## 65                                                                      
## 66                                                                      
## 67                                                                      
## 68                                                                      
## 69                                                                      
## 70                                                                      
## 71                                                                      
## 72                                                                      
## 73                                                                      
## 74  0                                                                   
## 75  0  0                                                                
## 76  0  0  0                                                             
## 77  0  0  0  0                                                          
## 78  0  0  0  0  0                                                       
## 79  0  0  0  0  0  0                                                    
## 80  0  0  0  0  0  0  0                                                 
## 81  0  0  0  0  0  0  0  0                                              
## 82  0  0  0  0  0  0  0  0  0                                           
## 83  0  0  0  0  0  0  0  0  0  0                                        
## 84  0  0  0  0  0  0  0  0  0  0  0                                     
## 85  0  0  0  0  0  0  0  0  0  0  0  0                                  
## 86  0  0  0  0  0  0  0  0  0  0  0  0  0                               
## 87  0  0  0  0  0  0  0  0  0  0  0  0  0  0                            
## 88  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                         
## 89  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                      
## 90  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                   
## 91  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                
## 92  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0             
## 93  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0          
## 94  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0       
## 95  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
## 96  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 97  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 98  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 99  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 100 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 101 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 102 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 103 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 104 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 105 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 106 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 107 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 108 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 109 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 110 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 111 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 112 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 113 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 114 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 115 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 116 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 117 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 118 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 119 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 120 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 121 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 122 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 123 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 124 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 125 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 126 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 127 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 128 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 129 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 130 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 131 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 132 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 133 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 134 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 135 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 136 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 137 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 138 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 139 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 140 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 141 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 142 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 143 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 144 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 145 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 146 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 147 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 148 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 149 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 150 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 151 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 152 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 153 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 154 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 155 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 156 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 157 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 158 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 159 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 160 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 161 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 162 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 163 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 164 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
##     96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113
## 2                                                                      
## 3                                                                      
## 4                                                                      
## 5                                                                      
## 6                                                                      
## 7                                                                      
## 8                                                                      
## 9                                                                      
## 10                                                                     
## 11                                                                     
## 12                                                                     
## 13                                                                     
## 14                                                                     
## 15                                                                     
## 16                                                                     
## 17                                                                     
## 18                                                                     
## 19                                                                     
## 20                                                                     
## 21                                                                     
## 22                                                                     
## 23                                                                     
## 24                                                                     
## 25                                                                     
## 26                                                                     
## 27                                                                     
## 28                                                                     
## 29                                                                     
## 30                                                                     
## 31                                                                     
## 32                                                                     
## 33                                                                     
## 34                                                                     
## 35                                                                     
## 36                                                                     
## 37                                                                     
## 38                                                                     
## 39                                                                     
## 40                                                                     
## 41                                                                     
## 42                                                                     
## 43                                                                     
## 44                                                                     
## 45                                                                     
## 46                                                                     
## 47                                                                     
## 48                                                                     
## 49                                                                     
## 50                                                                     
## 51                                                                     
## 52                                                                     
## 53                                                                     
## 54                                                                     
## 55                                                                     
## 56                                                                     
## 57                                                                     
## 58                                                                     
## 59                                                                     
## 60                                                                     
## 61                                                                     
## 62                                                                     
## 63                                                                     
## 64                                                                     
## 65                                                                     
## 66                                                                     
## 67                                                                     
## 68                                                                     
## 69                                                                     
## 70                                                                     
## 71                                                                     
## 72                                                                     
## 73                                                                     
## 74                                                                     
## 75                                                                     
## 76                                                                     
## 77                                                                     
## 78                                                                     
## 79                                                                     
## 80                                                                     
## 81                                                                     
## 82                                                                     
## 83                                                                     
## 84                                                                     
## 85                                                                     
## 86                                                                     
## 87                                                                     
## 88                                                                     
## 89                                                                     
## 90                                                                     
## 91                                                                     
## 92                                                                     
## 93                                                                     
## 94                                                                     
## 95                                                                     
## 96                                                                     
## 97  0                                                                  
## 98  0  0                                                               
## 99  0  0  0                                                            
## 100 0  0  0  0                                                         
## 101 0  0  0  0  0                                                      
## 102 0  0  0  0  0   0                                                  
## 103 0  0  0  0  0   0   0                                              
## 104 0  0  0  0  0   0   0   0                                          
## 105 0  0  0  0  0   0   0   0   0                                      
## 106 0  0  0  0  0   0   0   0   0   0                                  
## 107 0  0  0  0  0   0   0   0   0   0   0                              
## 108 0  0  0  0  0   0   0   0   0   0   0   0                          
## 109 0  0  0  0  0   0   0   0   0   0   0   0   0                      
## 110 0  0  0  0  0   0   0   0   0   0   0   0   0   0                  
## 111 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0              
## 112 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0          
## 113 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0      
## 114 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 115 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 116 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 117 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 118 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 119 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 120 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 121 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 122 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 123 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 124 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 125 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 126 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 127 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 128 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 129 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 130 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 131 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 132 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 133 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 134 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 135 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 136 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 137 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 138 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 139 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 140 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 141 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 142 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 143 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 144 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 145 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 146 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 147 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 148 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 149 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 150 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 151 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 152 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 153 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 154 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 155 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 156 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 157 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 158 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 159 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 160 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 161 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 162 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 163 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 164 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
##     114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130
## 2                                                                      
## 3                                                                      
## 4                                                                      
## 5                                                                      
## 6                                                                      
## 7                                                                      
## 8                                                                      
## 9                                                                      
## 10                                                                     
## 11                                                                     
## 12                                                                     
## 13                                                                     
## 14                                                                     
## 15                                                                     
## 16                                                                     
## 17                                                                     
## 18                                                                     
## 19                                                                     
## 20                                                                     
## 21                                                                     
## 22                                                                     
## 23                                                                     
## 24                                                                     
## 25                                                                     
## 26                                                                     
## 27                                                                     
## 28                                                                     
## 29                                                                     
## 30                                                                     
## 31                                                                     
## 32                                                                     
## 33                                                                     
## 34                                                                     
## 35                                                                     
## 36                                                                     
## 37                                                                     
## 38                                                                     
## 39                                                                     
## 40                                                                     
## 41                                                                     
## 42                                                                     
## 43                                                                     
## 44                                                                     
## 45                                                                     
## 46                                                                     
## 47                                                                     
## 48                                                                     
## 49                                                                     
## 50                                                                     
## 51                                                                     
## 52                                                                     
## 53                                                                     
## 54                                                                     
## 55                                                                     
## 56                                                                     
## 57                                                                     
## 58                                                                     
## 59                                                                     
## 60                                                                     
## 61                                                                     
## 62                                                                     
## 63                                                                     
## 64                                                                     
## 65                                                                     
## 66                                                                     
## 67                                                                     
## 68                                                                     
## 69                                                                     
## 70                                                                     
## 71                                                                     
## 72                                                                     
## 73                                                                     
## 74                                                                     
## 75                                                                     
## 76                                                                     
## 77                                                                     
## 78                                                                     
## 79                                                                     
## 80                                                                     
## 81                                                                     
## 82                                                                     
## 83                                                                     
## 84                                                                     
## 85                                                                     
## 86                                                                     
## 87                                                                     
## 88                                                                     
## 89                                                                     
## 90                                                                     
## 91                                                                     
## 92                                                                     
## 93                                                                     
## 94                                                                     
## 95                                                                     
## 96                                                                     
## 97                                                                     
## 98                                                                     
## 99                                                                     
## 100                                                                    
## 101                                                                    
## 102                                                                    
## 103                                                                    
## 104                                                                    
## 105                                                                    
## 106                                                                    
## 107                                                                    
## 108                                                                    
## 109                                                                    
## 110                                                                    
## 111                                                                    
## 112                                                                    
## 113                                                                    
## 114                                                                    
## 115 0                                                                  
## 116 0   0                                                              
## 117 0   0   0                                                          
## 118 0   0   0   0                                                      
## 119 0   0   0   0   0                                                  
## 120 0   0   0   0   0   0                                              
## 121 0   0   0   0   0   0   0                                          
## 122 0   0   0   0   0   0   0   0                                      
## 123 0   0   0   0   0   0   0   0   0                                  
## 124 0   0   0   0   0   0   0   0   0   0                              
## 125 0   0   0   0   0   0   0   0   0   0   0                          
## 126 0   0   0   0   0   0   0   0   0   0   0   0                      
## 127 0   0   0   0   0   0   0   0   0   0   0   0   0                  
## 128 0   0   0   0   0   0   0   0   0   0   0   0   0   0              
## 129 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0          
## 130 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0      
## 131 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 132 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 133 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 134 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 135 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 136 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 137 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 138 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 139 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 140 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 141 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 142 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 143 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 144 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 145 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 146 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 147 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 148 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 149 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 150 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 151 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 152 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 153 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 154 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 155 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 156 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 157 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 158 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 159 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 160 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 161 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 162 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 163 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 164 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
##     131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147
## 2                                                                      
## 3                                                                      
## 4                                                                      
## 5                                                                      
## 6                                                                      
## 7                                                                      
## 8                                                                      
## 9                                                                      
## 10                                                                     
## 11                                                                     
## 12                                                                     
## 13                                                                     
## 14                                                                     
## 15                                                                     
## 16                                                                     
## 17                                                                     
## 18                                                                     
## 19                                                                     
## 20                                                                     
## 21                                                                     
## 22                                                                     
## 23                                                                     
## 24                                                                     
## 25                                                                     
## 26                                                                     
## 27                                                                     
## 28                                                                     
## 29                                                                     
## 30                                                                     
## 31                                                                     
## 32                                                                     
## 33                                                                     
## 34                                                                     
## 35                                                                     
## 36                                                                     
## 37                                                                     
## 38                                                                     
## 39                                                                     
## 40                                                                     
## 41                                                                     
## 42                                                                     
## 43                                                                     
## 44                                                                     
## 45                                                                     
## 46                                                                     
## 47                                                                     
## 48                                                                     
## 49                                                                     
## 50                                                                     
## 51                                                                     
## 52                                                                     
## 53                                                                     
## 54                                                                     
## 55                                                                     
## 56                                                                     
## 57                                                                     
## 58                                                                     
## 59                                                                     
## 60                                                                     
## 61                                                                     
## 62                                                                     
## 63                                                                     
## 64                                                                     
## 65                                                                     
## 66                                                                     
## 67                                                                     
## 68                                                                     
## 69                                                                     
## 70                                                                     
## 71                                                                     
## 72                                                                     
## 73                                                                     
## 74                                                                     
## 75                                                                     
## 76                                                                     
## 77                                                                     
## 78                                                                     
## 79                                                                     
## 80                                                                     
## 81                                                                     
## 82                                                                     
## 83                                                                     
## 84                                                                     
## 85                                                                     
## 86                                                                     
## 87                                                                     
## 88                                                                     
## 89                                                                     
## 90                                                                     
## 91                                                                     
## 92                                                                     
## 93                                                                     
## 94                                                                     
## 95                                                                     
## 96                                                                     
## 97                                                                     
## 98                                                                     
## 99                                                                     
## 100                                                                    
## 101                                                                    
## 102                                                                    
## 103                                                                    
## 104                                                                    
## 105                                                                    
## 106                                                                    
## 107                                                                    
## 108                                                                    
## 109                                                                    
## 110                                                                    
## 111                                                                    
## 112                                                                    
## 113                                                                    
## 114                                                                    
## 115                                                                    
## 116                                                                    
## 117                                                                    
## 118                                                                    
## 119                                                                    
## 120                                                                    
## 121                                                                    
## 122                                                                    
## 123                                                                    
## 124                                                                    
## 125                                                                    
## 126                                                                    
## 127                                                                    
## 128                                                                    
## 129                                                                    
## 130                                                                    
## 131                                                                    
## 132 0                                                                  
## 133 0   0                                                              
## 134 0   0   0                                                          
## 135 0   0   0   0                                                      
## 136 0   0   0   0   0                                                  
## 137 0   0   0   0   0   0                                              
## 138 0   0   0   0   0   0   0                                          
## 139 0   0   0   0   0   0   0   0                                      
## 140 0   0   0   0   0   0   0   0   0                                  
## 141 0   0   0   0   0   0   0   0   0   0                              
## 142 0   0   0   0   0   0   0   0   0   0   0                          
## 143 0   0   0   0   0   0   0   0   0   0   0   0                      
## 144 0   0   0   0   0   0   0   0   0   0   0   0   0                  
## 145 0   0   0   0   0   0   0   0   0   0   0   0   0   0              
## 146 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0          
## 147 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0      
## 148 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 149 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 150 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 151 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 152 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 153 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 154 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 155 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 156 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 157 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 158 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 159 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 160 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 161 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 162 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 163 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 164 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
##     148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163
## 2                                                                  
## 3                                                                  
## 4                                                                  
## 5                                                                  
## 6                                                                  
## 7                                                                  
## 8                                                                  
## 9                                                                  
## 10                                                                 
## 11                                                                 
## 12                                                                 
## 13                                                                 
## 14                                                                 
## 15                                                                 
## 16                                                                 
## 17                                                                 
## 18                                                                 
## 19                                                                 
## 20                                                                 
## 21                                                                 
## 22                                                                 
## 23                                                                 
## 24                                                                 
## 25                                                                 
## 26                                                                 
## 27                                                                 
## 28                                                                 
## 29                                                                 
## 30                                                                 
## 31                                                                 
## 32                                                                 
## 33                                                                 
## 34                                                                 
## 35                                                                 
## 36                                                                 
## 37                                                                 
## 38                                                                 
## 39                                                                 
## 40                                                                 
## 41                                                                 
## 42                                                                 
## 43                                                                 
## 44                                                                 
## 45                                                                 
## 46                                                                 
## 47                                                                 
## 48                                                                 
## 49                                                                 
## 50                                                                 
## 51                                                                 
## 52                                                                 
## 53                                                                 
## 54                                                                 
## 55                                                                 
## 56                                                                 
## 57                                                                 
## 58                                                                 
## 59                                                                 
## 60                                                                 
## 61                                                                 
## 62                                                                 
## 63                                                                 
## 64                                                                 
## 65                                                                 
## 66                                                                 
## 67                                                                 
## 68                                                                 
## 69                                                                 
## 70                                                                 
## 71                                                                 
## 72                                                                 
## 73                                                                 
## 74                                                                 
## 75                                                                 
## 76                                                                 
## 77                                                                 
## 78                                                                 
## 79                                                                 
## 80                                                                 
## 81                                                                 
## 82                                                                 
## 83                                                                 
## 84                                                                 
## 85                                                                 
## 86                                                                 
## 87                                                                 
## 88                                                                 
## 89                                                                 
## 90                                                                 
## 91                                                                 
## 92                                                                 
## 93                                                                 
## 94                                                                 
## 95                                                                 
## 96                                                                 
## 97                                                                 
## 98                                                                 
## 99                                                                 
## 100                                                                
## 101                                                                
## 102                                                                
## 103                                                                
## 104                                                                
## 105                                                                
## 106                                                                
## 107                                                                
## 108                                                                
## 109                                                                
## 110                                                                
## 111                                                                
## 112                                                                
## 113                                                                
## 114                                                                
## 115                                                                
## 116                                                                
## 117                                                                
## 118                                                                
## 119                                                                
## 120                                                                
## 121                                                                
## 122                                                                
## 123                                                                
## 124                                                                
## 125                                                                
## 126                                                                
## 127                                                                
## 128                                                                
## 129                                                                
## 130                                                                
## 131                                                                
## 132                                                                
## 133                                                                
## 134                                                                
## 135                                                                
## 136                                                                
## 137                                                                
## 138                                                                
## 139                                                                
## 140                                                                
## 141                                                                
## 142                                                                
## 143                                                                
## 144                                                                
## 145                                                                
## 146                                                                
## 147                                                                
## 148                                                                
## 149 0                                                              
## 150 0   0                                                          
## 151 0   0   0                                                      
## 152 0   0   0   0                                                  
## 153 0   0   0   0   0                                              
## 154 0   0   0   0   0   0                                          
## 155 0   0   0   0   0   0   0                                      
## 156 0   0   0   0   0   0   0   0                                  
## 157 0   0   0   0   0   0   0   0   0                              
## 158 0   0   0   0   0   0   0   0   0   0                          
## 159 0   0   0   0   0   0   0   0   0   0   0                      
## 160 0   0   0   0   0   0   0   0   0   0   0   0                  
## 161 0   0   0   0   0   0   0   0   0   0   0   0   0              
## 162 0   0   0   0   0   0   0   0   0   0   0   0   0   0          
## 163 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0      
## 164 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## Fixed effects: list(fixed) 
##                                                   Value  Std.Error   DF
## (Intercept)                                  -1.4713162 0.07331735 6801
## genotypePos/Pos                               0.0742170 0.10414511   45
## activityExercise                             -0.0062902 0.06510023 6801
## metabolite_typeamino acids                    1.5020471 0.06213088 6801
## metabolite_typeorganic acids                  1.1050970 0.06885615 6801
## genotypePos/Pos:activityExercise              0.0134917 0.09258194 6801
## genotypePos/Pos:metabolite_typeamino acids   -0.0853005 0.08798927 6801
## genotypePos/Pos:metabolite_typeorganic acids  0.0066969 0.09762673 6801
##                                                 t-value p-value
## (Intercept)                                  -20.067775  0.0000
## genotypePos/Pos                                0.712631  0.4798
## activityExercise                              -0.096623  0.9230
## metabolite_typeamino acids                    24.175532  0.0000
## metabolite_typeorganic acids                  16.049358  0.0000
## genotypePos/Pos:activityExercise               0.145727  0.8841
## genotypePos/Pos:metabolite_typeamino acids    -0.969443  0.3324
## genotypePos/Pos:metabolite_typeorganic acids   0.068597  0.9453
##  Correlation: 
##                                              (Intr) gntP/P actvtE
## genotypePos/Pos                              -0.704              
## activityExercise                             -0.888  0.625       
## metabolite_typeamino acids                   -0.693  0.488  0.499
## metabolite_typeorganic acids                 -0.625  0.440  0.450
## genotypePos/Pos:activityExercise              0.624 -0.889 -0.703
## genotypePos/Pos:metabolite_typeamino acids    0.489 -0.694 -0.352
## genotypePos/Pos:metabolite_typeorganic acids  0.441 -0.614 -0.318
##                                              mtblt_typma mtblt_typra
## genotypePos/Pos                                                     
## activityExercise                                                    
## metabolite_typeamino acids                                          
## metabolite_typeorganic acids                  0.490                 
## genotypePos/Pos:activityExercise             -0.351      -0.317     
## genotypePos/Pos:metabolite_typeamino acids   -0.706      -0.346     
## genotypePos/Pos:metabolite_typeorganic acids -0.346      -0.705     
##                                              gP/P:E gntypPs/Ps:mtblt_typma
## genotypePos/Pos                                                           
## activityExercise                                                          
## metabolite_typeamino acids                                                
## metabolite_typeorganic acids                                              
## genotypePos/Pos:activityExercise                                          
## genotypePos/Pos:metabolite_typeamino acids    0.501                       
## genotypePos/Pos:metabolite_typeorganic acids  0.439  0.485                
## 
## Standardized Within-Group Residuals:
##        Min         Q1        Med         Q3        Max 
## -5.0852593 -0.5139823  0.0642195  0.4939858  3.8402323 
## 
## Number of Observations: 6854
## Number of Groups: 47
```


## Aim 2: Chow

**MCT chow versus experimental chow with triheptanoin**

Estimate model.
Default to `correlation = NULL`, corresponding to no within-group correlations.


```r
fixed <- formula(logValue ~
                   genotype +
                   chow +
                   metabolite_type +
                   genotype * chow +
                   genotype * metabolite_type)
random <- formula(~ 1 | id)
M <- lme(fixed, data = D2, random = random, correlation = NULL)
summary(M)
```

```
## Linear mixed-effects model fit by REML
##  Data: D2 
##        AIC      BIC    logLik
##   26233.49 26302.62 -13106.75
## 
## Random effects:
##  Formula: ~1 | id
##          (Intercept) Residual
## StdDev: 9.091544e-05 1.407474
## 
## Fixed effects: list(fixed) 
##                                                   Value  Std.Error   DF
## (Intercept)                                  -1.3072088 0.04005958 7388
## genotypePos/Pos                              -0.1382761 0.05665281   40
## chowYellow (C8)                              -0.1224580 0.04616532   40
## metabolite_typeAmino Acids                    1.4396406 0.05359902 7388
## metabolite_typeorganic acids                  1.0233612 0.06037107 7388
## genotypePos/Pos:chowYellow (C8)               0.0580403 0.06528762   40
## genotypePos/Pos:metabolite_typeAmino Acids    0.0753990 0.07580046 7388
## genotypePos/Pos:metabolite_typeorganic acids  0.1235755 0.08537759 7388
##                                                t-value p-value
## (Intercept)                                  -32.63161  0.0000
## genotypePos/Pos                               -2.44076  0.0192
## chowYellow (C8)                               -2.65260  0.0114
## metabolite_typeAmino Acids                    26.85946  0.0000
## metabolite_typeorganic acids                  16.95118  0.0000
## genotypePos/Pos:chowYellow (C8)                0.88899  0.3793
## genotypePos/Pos:metabolite_typeAmino Acids     0.99470  0.3199
## genotypePos/Pos:metabolite_typeorganic acids   1.44740  0.1478
##  Correlation: 
##                                              (Intr) gntP/P cY(C8) mtb_AA
## genotypePos/Pos                              -0.707                     
## chowYellow (C8)                              -0.576  0.407              
## metabolite_typeAmino Acids                   -0.499  0.353  0.000       
## metabolite_typeorganic acids                 -0.443  0.313  0.000  0.331
## genotypePos/Pos:chowYellow (C8)               0.407 -0.576 -0.707  0.000
## genotypePos/Pos:metabolite_typeAmino Acids    0.353 -0.499  0.000 -0.707
## genotypePos/Pos:metabolite_typeorganic acids  0.313 -0.443  0.000 -0.234
##                                              mtbl_a gP/P:( gP/P:A
## genotypePos/Pos                                                  
## chowYellow (C8)                                                  
## metabolite_typeAmino Acids                                       
## metabolite_typeorganic acids                                     
## genotypePos/Pos:chowYellow (C8)               0.000              
## genotypePos/Pos:metabolite_typeAmino Acids   -0.234  0.000       
## genotypePos/Pos:metabolite_typeorganic acids -0.707  0.000  0.331
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -5.03533870 -0.47485373  0.05099419  0.50934516  4.29841962 
## 
## Number of Observations: 7436
## Number of Groups: 44
```

Estimate model.
Specify the correlation structure using `cs`.
Use `corSymm`, *general correlation matrix, with no additional structure*.


```r
cs <-
  corSymm(form = random, fixed = TRUE) %>%
  Initialize(data = D2)
Dim(cs)
```

```
## $N
## [1] 7436
## 
## $M
## [1] 44
## 
## $maxLen
## [1] 169
## 
## $sumLenSq
## [1] 1256684
## 
## $len
## groups
##       1120       1126       1127       1128       1142       1143 
##        169        169        169        169        169        169 
##       1144       1145       1146       1158       1159       1176 
##        169        169        169        169        169        169 
##       1177       1179       1180       1190       1191       1193 
##        169        169        169        169        169        169 
##       1194       1199       1200 1192A/1198       1107       1113 
##        169        169        169        169        169        169 
##       1114       1115       1117       1118       1119       1204 
##        169        169        169        169        169        169 
##       1205       1206       1208       1101       1102       1103 
##        169        169        169        169        169        169 
##       1184       1185       1186       1195       1196       1197 
##        169        169        169        169        169        169 
##       1203 1192/1198B 
##        169        169 
## 
## $start
##  [1]    0   84  168  252  336  420  504  588  672  756  840  924 1008 1092
## [15] 1176 1260 1344 1428 1512 1596 1680 1764 1848 1932 2016 2100 2184 2268
## [29] 2352 2436 2520 2604 2688 2772 2856 2940 3024 3108 3192 3276 3360 3444
## [43] 3528 3612
```

```r
M <- lme(fixed, data = D2, random = random, correlation = cs)
summary(M)
```

```
## Linear mixed-effects model fit by REML
##  Data: D2 
##        AIC      BIC    logLik
##   26233.49 26302.62 -13106.75
## 
## Random effects:
##  Formula: ~1 | id
##          (Intercept) Residual
## StdDev: 9.091859e-05 1.407474
## 
## Correlation Structure: General
##  Formula: ~1 | id 
##  Parameter estimate(s):
##  Correlation: 
##     1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
## 2   0                                                                   
## 3   0 0                                                                 
## 4   0 0 0                                                               
## 5   0 0 0 0                                                             
## 6   0 0 0 0 0                                                           
## 7   0 0 0 0 0 0                                                         
## 8   0 0 0 0 0 0 0                                                       
## 9   0 0 0 0 0 0 0 0                                                     
## 10  0 0 0 0 0 0 0 0 0                                                   
## 11  0 0 0 0 0 0 0 0 0 0                                                 
## 12  0 0 0 0 0 0 0 0 0 0  0                                              
## 13  0 0 0 0 0 0 0 0 0 0  0  0                                           
## 14  0 0 0 0 0 0 0 0 0 0  0  0  0                                        
## 15  0 0 0 0 0 0 0 0 0 0  0  0  0  0                                     
## 16  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0                                  
## 17  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0                               
## 18  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0                            
## 19  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0                         
## 20  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0                      
## 21  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0                   
## 22  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0                
## 23  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0             
## 24  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0          
## 25  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0       
## 26  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
## 27  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 28  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 29  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 30  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 31  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 32  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 33  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 34  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 35  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 36  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 37  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 38  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 39  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 40  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 41  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 42  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 43  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 44  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 45  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 46  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 47  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 48  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 49  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 50  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 51  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 52  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 53  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 54  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 55  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 56  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 57  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 58  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 59  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 60  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 61  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 62  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 63  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 64  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 65  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 66  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 67  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 68  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 69  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 70  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 71  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 72  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 73  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 74  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 75  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 76  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 77  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 78  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 79  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 80  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 81  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 82  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 83  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 84  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 85  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 86  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 87  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 88  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 89  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 90  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 91  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 92  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 93  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 94  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 95  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 96  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 97  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 98  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 99  0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 100 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 101 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 102 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 103 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 104 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 105 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 106 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 107 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 108 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 109 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 110 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 111 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 112 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 113 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 114 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 115 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 116 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 117 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 118 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 119 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 120 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 121 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 122 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 123 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 124 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 125 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 126 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 127 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 128 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 129 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 130 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 131 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 132 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 133 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 134 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 135 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 136 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 137 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 138 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 139 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 140 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 141 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 142 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 143 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 144 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 145 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 146 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 147 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 148 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 149 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 150 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 151 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 152 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 153 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 154 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 155 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 156 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 157 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 158 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 159 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 160 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 161 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 162 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 163 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 164 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 165 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 166 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 167 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 168 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 169 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
##     27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49
## 2                                                                       
## 3                                                                       
## 4                                                                       
## 5                                                                       
## 6                                                                       
## 7                                                                       
## 8                                                                       
## 9                                                                       
## 10                                                                      
## 11                                                                      
## 12                                                                      
## 13                                                                      
## 14                                                                      
## 15                                                                      
## 16                                                                      
## 17                                                                      
## 18                                                                      
## 19                                                                      
## 20                                                                      
## 21                                                                      
## 22                                                                      
## 23                                                                      
## 24                                                                      
## 25                                                                      
## 26                                                                      
## 27                                                                      
## 28  0                                                                   
## 29  0  0                                                                
## 30  0  0  0                                                             
## 31  0  0  0  0                                                          
## 32  0  0  0  0  0                                                       
## 33  0  0  0  0  0  0                                                    
## 34  0  0  0  0  0  0  0                                                 
## 35  0  0  0  0  0  0  0  0                                              
## 36  0  0  0  0  0  0  0  0  0                                           
## 37  0  0  0  0  0  0  0  0  0  0                                        
## 38  0  0  0  0  0  0  0  0  0  0  0                                     
## 39  0  0  0  0  0  0  0  0  0  0  0  0                                  
## 40  0  0  0  0  0  0  0  0  0  0  0  0  0                               
## 41  0  0  0  0  0  0  0  0  0  0  0  0  0  0                            
## 42  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                         
## 43  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                      
## 44  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                   
## 45  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                
## 46  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0             
## 47  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0          
## 48  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0       
## 49  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
## 50  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 51  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 52  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 53  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 54  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 55  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 56  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 57  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 58  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 59  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 60  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 61  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 63  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 64  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 65  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 66  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 67  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 68  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 69  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 70  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 71  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 72  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 73  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 74  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 75  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 76  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 77  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 78  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 79  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 80  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 81  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 82  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 83  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 84  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 85  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 86  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 87  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 88  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 89  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 90  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 91  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 92  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 93  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 94  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 95  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 96  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 97  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 98  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 99  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 100 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 101 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 102 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 103 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 104 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 105 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 106 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 107 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 108 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 109 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 110 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 111 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 112 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 113 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 114 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 115 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 116 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 117 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 118 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 119 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 120 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 121 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 122 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 123 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 124 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 125 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 126 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 127 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 128 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 129 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 130 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 131 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 132 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 133 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 134 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 135 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 136 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 137 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 138 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 139 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 140 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 141 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 142 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 143 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 144 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 145 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 146 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 147 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 148 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 149 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 150 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 151 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 152 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 153 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 154 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 155 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 156 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 157 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 158 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 159 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 160 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 161 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 162 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 163 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 164 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 165 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 166 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 167 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 168 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 169 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
##     50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72
## 2                                                                       
## 3                                                                       
## 4                                                                       
## 5                                                                       
## 6                                                                       
## 7                                                                       
## 8                                                                       
## 9                                                                       
## 10                                                                      
## 11                                                                      
## 12                                                                      
## 13                                                                      
## 14                                                                      
## 15                                                                      
## 16                                                                      
## 17                                                                      
## 18                                                                      
## 19                                                                      
## 20                                                                      
## 21                                                                      
## 22                                                                      
## 23                                                                      
## 24                                                                      
## 25                                                                      
## 26                                                                      
## 27                                                                      
## 28                                                                      
## 29                                                                      
## 30                                                                      
## 31                                                                      
## 32                                                                      
## 33                                                                      
## 34                                                                      
## 35                                                                      
## 36                                                                      
## 37                                                                      
## 38                                                                      
## 39                                                                      
## 40                                                                      
## 41                                                                      
## 42                                                                      
## 43                                                                      
## 44                                                                      
## 45                                                                      
## 46                                                                      
## 47                                                                      
## 48                                                                      
## 49                                                                      
## 50                                                                      
## 51  0                                                                   
## 52  0  0                                                                
## 53  0  0  0                                                             
## 54  0  0  0  0                                                          
## 55  0  0  0  0  0                                                       
## 56  0  0  0  0  0  0                                                    
## 57  0  0  0  0  0  0  0                                                 
## 58  0  0  0  0  0  0  0  0                                              
## 59  0  0  0  0  0  0  0  0  0                                           
## 60  0  0  0  0  0  0  0  0  0  0                                        
## 61  0  0  0  0  0  0  0  0  0  0  0                                     
## 62  0  0  0  0  0  0  0  0  0  0  0  0                                  
## 63  0  0  0  0  0  0  0  0  0  0  0  0  0                               
## 64  0  0  0  0  0  0  0  0  0  0  0  0  0  0                            
## 65  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                         
## 66  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                      
## 67  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                   
## 68  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                
## 69  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0             
## 70  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0          
## 71  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0       
## 72  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
## 73  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 74  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 75  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 76  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 77  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 78  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 79  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 80  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 81  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 82  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 83  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 84  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 85  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 86  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 87  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 88  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 89  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 90  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 91  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 92  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 93  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 94  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 95  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 96  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 97  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 98  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 99  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 100 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 101 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 102 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 103 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 104 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 105 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 106 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 107 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 108 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 109 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 110 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 111 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 112 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 113 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 114 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 115 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 116 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 117 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 118 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 119 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 120 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 121 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 122 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 123 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 124 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 125 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 126 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 127 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 128 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 129 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 130 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 131 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 132 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 133 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 134 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 135 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 136 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 137 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 138 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 139 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 140 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 141 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 142 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 143 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 144 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 145 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 146 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 147 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 148 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 149 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 150 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 151 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 152 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 153 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 154 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 155 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 156 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 157 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 158 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 159 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 160 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 161 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 162 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 163 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 164 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 165 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 166 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 167 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 168 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 169 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
##     73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95
## 2                                                                       
## 3                                                                       
## 4                                                                       
## 5                                                                       
## 6                                                                       
## 7                                                                       
## 8                                                                       
## 9                                                                       
## 10                                                                      
## 11                                                                      
## 12                                                                      
## 13                                                                      
## 14                                                                      
## 15                                                                      
## 16                                                                      
## 17                                                                      
## 18                                                                      
## 19                                                                      
## 20                                                                      
## 21                                                                      
## 22                                                                      
## 23                                                                      
## 24                                                                      
## 25                                                                      
## 26                                                                      
## 27                                                                      
## 28                                                                      
## 29                                                                      
## 30                                                                      
## 31                                                                      
## 32                                                                      
## 33                                                                      
## 34                                                                      
## 35                                                                      
## 36                                                                      
## 37                                                                      
## 38                                                                      
## 39                                                                      
## 40                                                                      
## 41                                                                      
## 42                                                                      
## 43                                                                      
## 44                                                                      
## 45                                                                      
## 46                                                                      
## 47                                                                      
## 48                                                                      
## 49                                                                      
## 50                                                                      
## 51                                                                      
## 52                                                                      
## 53                                                                      
## 54                                                                      
## 55                                                                      
## 56                                                                      
## 57                                                                      
## 58                                                                      
## 59                                                                      
## 60                                                                      
## 61                                                                      
## 62                                                                      
## 63                                                                      
## 64                                                                      
## 65                                                                      
## 66                                                                      
## 67                                                                      
## 68                                                                      
## 69                                                                      
## 70                                                                      
## 71                                                                      
## 72                                                                      
## 73                                                                      
## 74  0                                                                   
## 75  0  0                                                                
## 76  0  0  0                                                             
## 77  0  0  0  0                                                          
## 78  0  0  0  0  0                                                       
## 79  0  0  0  0  0  0                                                    
## 80  0  0  0  0  0  0  0                                                 
## 81  0  0  0  0  0  0  0  0                                              
## 82  0  0  0  0  0  0  0  0  0                                           
## 83  0  0  0  0  0  0  0  0  0  0                                        
## 84  0  0  0  0  0  0  0  0  0  0  0                                     
## 85  0  0  0  0  0  0  0  0  0  0  0  0                                  
## 86  0  0  0  0  0  0  0  0  0  0  0  0  0                               
## 87  0  0  0  0  0  0  0  0  0  0  0  0  0  0                            
## 88  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                         
## 89  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                      
## 90  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                   
## 91  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0                
## 92  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0             
## 93  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0          
## 94  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0       
## 95  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0    
## 96  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 97  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 98  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 99  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 100 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 101 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 102 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 103 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 104 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 105 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 106 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 107 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 108 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 109 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 110 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 111 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 112 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 113 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 114 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 115 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 116 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 117 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 118 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 119 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 120 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 121 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 122 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 123 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 124 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 125 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 126 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 127 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 128 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 129 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 130 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 131 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 132 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 133 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 134 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 135 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 136 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 137 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 138 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 139 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 140 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 141 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 142 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 143 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 144 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 145 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 146 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 147 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 148 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 149 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 150 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 151 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 152 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 153 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 154 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 155 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 156 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 157 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 158 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 159 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 160 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 161 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 162 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 163 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 164 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 165 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 166 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 167 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 168 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
## 169 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
##     96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113
## 2                                                                      
## 3                                                                      
## 4                                                                      
## 5                                                                      
## 6                                                                      
## 7                                                                      
## 8                                                                      
## 9                                                                      
## 10                                                                     
## 11                                                                     
## 12                                                                     
## 13                                                                     
## 14                                                                     
## 15                                                                     
## 16                                                                     
## 17                                                                     
## 18                                                                     
## 19                                                                     
## 20                                                                     
## 21                                                                     
## 22                                                                     
## 23                                                                     
## 24                                                                     
## 25                                                                     
## 26                                                                     
## 27                                                                     
## 28                                                                     
## 29                                                                     
## 30                                                                     
## 31                                                                     
## 32                                                                     
## 33                                                                     
## 34                                                                     
## 35                                                                     
## 36                                                                     
## 37                                                                     
## 38                                                                     
## 39                                                                     
## 40                                                                     
## 41                                                                     
## 42                                                                     
## 43                                                                     
## 44                                                                     
## 45                                                                     
## 46                                                                     
## 47                                                                     
## 48                                                                     
## 49                                                                     
## 50                                                                     
## 51                                                                     
## 52                                                                     
## 53                                                                     
## 54                                                                     
## 55                                                                     
## 56                                                                     
## 57                                                                     
## 58                                                                     
## 59                                                                     
## 60                                                                     
## 61                                                                     
## 62                                                                     
## 63                                                                     
## 64                                                                     
## 65                                                                     
## 66                                                                     
## 67                                                                     
## 68                                                                     
## 69                                                                     
## 70                                                                     
## 71                                                                     
## 72                                                                     
## 73                                                                     
## 74                                                                     
## 75                                                                     
## 76                                                                     
## 77                                                                     
## 78                                                                     
## 79                                                                     
## 80                                                                     
## 81                                                                     
## 82                                                                     
## 83                                                                     
## 84                                                                     
## 85                                                                     
## 86                                                                     
## 87                                                                     
## 88                                                                     
## 89                                                                     
## 90                                                                     
## 91                                                                     
## 92                                                                     
## 93                                                                     
## 94                                                                     
## 95                                                                     
## 96                                                                     
## 97  0                                                                  
## 98  0  0                                                               
## 99  0  0  0                                                            
## 100 0  0  0  0                                                         
## 101 0  0  0  0  0                                                      
## 102 0  0  0  0  0   0                                                  
## 103 0  0  0  0  0   0   0                                              
## 104 0  0  0  0  0   0   0   0                                          
## 105 0  0  0  0  0   0   0   0   0                                      
## 106 0  0  0  0  0   0   0   0   0   0                                  
## 107 0  0  0  0  0   0   0   0   0   0   0                              
## 108 0  0  0  0  0   0   0   0   0   0   0   0                          
## 109 0  0  0  0  0   0   0   0   0   0   0   0   0                      
## 110 0  0  0  0  0   0   0   0   0   0   0   0   0   0                  
## 111 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0              
## 112 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0          
## 113 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0      
## 114 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 115 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 116 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 117 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 118 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 119 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 120 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 121 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 122 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 123 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 124 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 125 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 126 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 127 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 128 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 129 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 130 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 131 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 132 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 133 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 134 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 135 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 136 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 137 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 138 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 139 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 140 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 141 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 142 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 143 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 144 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 145 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 146 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 147 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 148 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 149 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 150 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 151 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 152 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 153 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 154 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 155 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 156 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 157 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 158 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 159 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 160 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 161 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 162 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 163 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 164 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 165 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 166 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 167 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 168 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 169 0  0  0  0  0   0   0   0   0   0   0   0   0   0   0   0   0   0  
##     114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130
## 2                                                                      
## 3                                                                      
## 4                                                                      
## 5                                                                      
## 6                                                                      
## 7                                                                      
## 8                                                                      
## 9                                                                      
## 10                                                                     
## 11                                                                     
## 12                                                                     
## 13                                                                     
## 14                                                                     
## 15                                                                     
## 16                                                                     
## 17                                                                     
## 18                                                                     
## 19                                                                     
## 20                                                                     
## 21                                                                     
## 22                                                                     
## 23                                                                     
## 24                                                                     
## 25                                                                     
## 26                                                                     
## 27                                                                     
## 28                                                                     
## 29                                                                     
## 30                                                                     
## 31                                                                     
## 32                                                                     
## 33                                                                     
## 34                                                                     
## 35                                                                     
## 36                                                                     
## 37                                                                     
## 38                                                                     
## 39                                                                     
## 40                                                                     
## 41                                                                     
## 42                                                                     
## 43                                                                     
## 44                                                                     
## 45                                                                     
## 46                                                                     
## 47                                                                     
## 48                                                                     
## 49                                                                     
## 50                                                                     
## 51                                                                     
## 52                                                                     
## 53                                                                     
## 54                                                                     
## 55                                                                     
## 56                                                                     
## 57                                                                     
## 58                                                                     
## 59                                                                     
## 60                                                                     
## 61                                                                     
## 62                                                                     
## 63                                                                     
## 64                                                                     
## 65                                                                     
## 66                                                                     
## 67                                                                     
## 68                                                                     
## 69                                                                     
## 70                                                                     
## 71                                                                     
## 72                                                                     
## 73                                                                     
## 74                                                                     
## 75                                                                     
## 76                                                                     
## 77                                                                     
## 78                                                                     
## 79                                                                     
## 80                                                                     
## 81                                                                     
## 82                                                                     
## 83                                                                     
## 84                                                                     
## 85                                                                     
## 86                                                                     
## 87                                                                     
## 88                                                                     
## 89                                                                     
## 90                                                                     
## 91                                                                     
## 92                                                                     
## 93                                                                     
## 94                                                                     
## 95                                                                     
## 96                                                                     
## 97                                                                     
## 98                                                                     
## 99                                                                     
## 100                                                                    
## 101                                                                    
## 102                                                                    
## 103                                                                    
## 104                                                                    
## 105                                                                    
## 106                                                                    
## 107                                                                    
## 108                                                                    
## 109                                                                    
## 110                                                                    
## 111                                                                    
## 112                                                                    
## 113                                                                    
## 114                                                                    
## 115 0                                                                  
## 116 0   0                                                              
## 117 0   0   0                                                          
## 118 0   0   0   0                                                      
## 119 0   0   0   0   0                                                  
## 120 0   0   0   0   0   0                                              
## 121 0   0   0   0   0   0   0                                          
## 122 0   0   0   0   0   0   0   0                                      
## 123 0   0   0   0   0   0   0   0   0                                  
## 124 0   0   0   0   0   0   0   0   0   0                              
## 125 0   0   0   0   0   0   0   0   0   0   0                          
## 126 0   0   0   0   0   0   0   0   0   0   0   0                      
## 127 0   0   0   0   0   0   0   0   0   0   0   0   0                  
## 128 0   0   0   0   0   0   0   0   0   0   0   0   0   0              
## 129 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0          
## 130 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0      
## 131 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 132 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 133 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 134 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 135 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 136 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 137 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 138 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 139 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 140 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 141 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 142 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 143 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 144 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 145 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 146 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 147 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 148 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 149 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 150 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 151 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 152 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 153 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 154 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 155 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 156 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 157 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 158 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 159 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 160 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 161 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 162 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 163 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 164 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 165 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 166 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 167 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 168 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 169 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
##     131 132 133 134 135 136 137 138 139 140 141 142 143 144 145 146 147
## 2                                                                      
## 3                                                                      
## 4                                                                      
## 5                                                                      
## 6                                                                      
## 7                                                                      
## 8                                                                      
## 9                                                                      
## 10                                                                     
## 11                                                                     
## 12                                                                     
## 13                                                                     
## 14                                                                     
## 15                                                                     
## 16                                                                     
## 17                                                                     
## 18                                                                     
## 19                                                                     
## 20                                                                     
## 21                                                                     
## 22                                                                     
## 23                                                                     
## 24                                                                     
## 25                                                                     
## 26                                                                     
## 27                                                                     
## 28                                                                     
## 29                                                                     
## 30                                                                     
## 31                                                                     
## 32                                                                     
## 33                                                                     
## 34                                                                     
## 35                                                                     
## 36                                                                     
## 37                                                                     
## 38                                                                     
## 39                                                                     
## 40                                                                     
## 41                                                                     
## 42                                                                     
## 43                                                                     
## 44                                                                     
## 45                                                                     
## 46                                                                     
## 47                                                                     
## 48                                                                     
## 49                                                                     
## 50                                                                     
## 51                                                                     
## 52                                                                     
## 53                                                                     
## 54                                                                     
## 55                                                                     
## 56                                                                     
## 57                                                                     
## 58                                                                     
## 59                                                                     
## 60                                                                     
## 61                                                                     
## 62                                                                     
## 63                                                                     
## 64                                                                     
## 65                                                                     
## 66                                                                     
## 67                                                                     
## 68                                                                     
## 69                                                                     
## 70                                                                     
## 71                                                                     
## 72                                                                     
## 73                                                                     
## 74                                                                     
## 75                                                                     
## 76                                                                     
## 77                                                                     
## 78                                                                     
## 79                                                                     
## 80                                                                     
## 81                                                                     
## 82                                                                     
## 83                                                                     
## 84                                                                     
## 85                                                                     
## 86                                                                     
## 87                                                                     
## 88                                                                     
## 89                                                                     
## 90                                                                     
## 91                                                                     
## 92                                                                     
## 93                                                                     
## 94                                                                     
## 95                                                                     
## 96                                                                     
## 97                                                                     
## 98                                                                     
## 99                                                                     
## 100                                                                    
## 101                                                                    
## 102                                                                    
## 103                                                                    
## 104                                                                    
## 105                                                                    
## 106                                                                    
## 107                                                                    
## 108                                                                    
## 109                                                                    
## 110                                                                    
## 111                                                                    
## 112                                                                    
## 113                                                                    
## 114                                                                    
## 115                                                                    
## 116                                                                    
## 117                                                                    
## 118                                                                    
## 119                                                                    
## 120                                                                    
## 121                                                                    
## 122                                                                    
## 123                                                                    
## 124                                                                    
## 125                                                                    
## 126                                                                    
## 127                                                                    
## 128                                                                    
## 129                                                                    
## 130                                                                    
## 131                                                                    
## 132 0                                                                  
## 133 0   0                                                              
## 134 0   0   0                                                          
## 135 0   0   0   0                                                      
## 136 0   0   0   0   0                                                  
## 137 0   0   0   0   0   0                                              
## 138 0   0   0   0   0   0   0                                          
## 139 0   0   0   0   0   0   0   0                                      
## 140 0   0   0   0   0   0   0   0   0                                  
## 141 0   0   0   0   0   0   0   0   0   0                              
## 142 0   0   0   0   0   0   0   0   0   0   0                          
## 143 0   0   0   0   0   0   0   0   0   0   0   0                      
## 144 0   0   0   0   0   0   0   0   0   0   0   0   0                  
## 145 0   0   0   0   0   0   0   0   0   0   0   0   0   0              
## 146 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0          
## 147 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0      
## 148 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 149 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 150 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 151 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 152 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 153 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 154 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 155 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 156 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 157 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 158 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 159 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 160 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 161 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 162 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 163 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 164 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 165 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 166 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 167 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 168 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 169 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
##     148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164
## 2                                                                      
## 3                                                                      
## 4                                                                      
## 5                                                                      
## 6                                                                      
## 7                                                                      
## 8                                                                      
## 9                                                                      
## 10                                                                     
## 11                                                                     
## 12                                                                     
## 13                                                                     
## 14                                                                     
## 15                                                                     
## 16                                                                     
## 17                                                                     
## 18                                                                     
## 19                                                                     
## 20                                                                     
## 21                                                                     
## 22                                                                     
## 23                                                                     
## 24                                                                     
## 25                                                                     
## 26                                                                     
## 27                                                                     
## 28                                                                     
## 29                                                                     
## 30                                                                     
## 31                                                                     
## 32                                                                     
## 33                                                                     
## 34                                                                     
## 35                                                                     
## 36                                                                     
## 37                                                                     
## 38                                                                     
## 39                                                                     
## 40                                                                     
## 41                                                                     
## 42                                                                     
## 43                                                                     
## 44                                                                     
## 45                                                                     
## 46                                                                     
## 47                                                                     
## 48                                                                     
## 49                                                                     
## 50                                                                     
## 51                                                                     
## 52                                                                     
## 53                                                                     
## 54                                                                     
## 55                                                                     
## 56                                                                     
## 57                                                                     
## 58                                                                     
## 59                                                                     
## 60                                                                     
## 61                                                                     
## 62                                                                     
## 63                                                                     
## 64                                                                     
## 65                                                                     
## 66                                                                     
## 67                                                                     
## 68                                                                     
## 69                                                                     
## 70                                                                     
## 71                                                                     
## 72                                                                     
## 73                                                                     
## 74                                                                     
## 75                                                                     
## 76                                                                     
## 77                                                                     
## 78                                                                     
## 79                                                                     
## 80                                                                     
## 81                                                                     
## 82                                                                     
## 83                                                                     
## 84                                                                     
## 85                                                                     
## 86                                                                     
## 87                                                                     
## 88                                                                     
## 89                                                                     
## 90                                                                     
## 91                                                                     
## 92                                                                     
## 93                                                                     
## 94                                                                     
## 95                                                                     
## 96                                                                     
## 97                                                                     
## 98                                                                     
## 99                                                                     
## 100                                                                    
## 101                                                                    
## 102                                                                    
## 103                                                                    
## 104                                                                    
## 105                                                                    
## 106                                                                    
## 107                                                                    
## 108                                                                    
## 109                                                                    
## 110                                                                    
## 111                                                                    
## 112                                                                    
## 113                                                                    
## 114                                                                    
## 115                                                                    
## 116                                                                    
## 117                                                                    
## 118                                                                    
## 119                                                                    
## 120                                                                    
## 121                                                                    
## 122                                                                    
## 123                                                                    
## 124                                                                    
## 125                                                                    
## 126                                                                    
## 127                                                                    
## 128                                                                    
## 129                                                                    
## 130                                                                    
## 131                                                                    
## 132                                                                    
## 133                                                                    
## 134                                                                    
## 135                                                                    
## 136                                                                    
## 137                                                                    
## 138                                                                    
## 139                                                                    
## 140                                                                    
## 141                                                                    
## 142                                                                    
## 143                                                                    
## 144                                                                    
## 145                                                                    
## 146                                                                    
## 147                                                                    
## 148                                                                    
## 149 0                                                                  
## 150 0   0                                                              
## 151 0   0   0                                                          
## 152 0   0   0   0                                                      
## 153 0   0   0   0   0                                                  
## 154 0   0   0   0   0   0                                              
## 155 0   0   0   0   0   0   0                                          
## 156 0   0   0   0   0   0   0   0                                      
## 157 0   0   0   0   0   0   0   0   0                                  
## 158 0   0   0   0   0   0   0   0   0   0                              
## 159 0   0   0   0   0   0   0   0   0   0   0                          
## 160 0   0   0   0   0   0   0   0   0   0   0   0                      
## 161 0   0   0   0   0   0   0   0   0   0   0   0   0                  
## 162 0   0   0   0   0   0   0   0   0   0   0   0   0   0              
## 163 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0          
## 164 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0      
## 165 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 166 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 167 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 168 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
## 169 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  
##     165 166 167 168
## 2                  
## 3                  
## 4                  
## 5                  
## 6                  
## 7                  
## 8                  
## 9                  
## 10                 
## 11                 
## 12                 
## 13                 
## 14                 
## 15                 
## 16                 
## 17                 
## 18                 
## 19                 
## 20                 
## 21                 
## 22                 
## 23                 
## 24                 
## 25                 
## 26                 
## 27                 
## 28                 
## 29                 
## 30                 
## 31                 
## 32                 
## 33                 
## 34                 
## 35                 
## 36                 
## 37                 
## 38                 
## 39                 
## 40                 
## 41                 
## 42                 
## 43                 
## 44                 
## 45                 
## 46                 
## 47                 
## 48                 
## 49                 
## 50                 
## 51                 
## 52                 
## 53                 
## 54                 
## 55                 
## 56                 
## 57                 
## 58                 
## 59                 
## 60                 
## 61                 
## 62                 
## 63                 
## 64                 
## 65                 
## 66                 
## 67                 
## 68                 
## 69                 
## 70                 
## 71                 
## 72                 
## 73                 
## 74                 
## 75                 
## 76                 
## 77                 
## 78                 
## 79                 
## 80                 
## 81                 
## 82                 
## 83                 
## 84                 
## 85                 
## 86                 
## 87                 
## 88                 
## 89                 
## 90                 
## 91                 
## 92                 
## 93                 
## 94                 
## 95                 
## 96                 
## 97                 
## 98                 
## 99                 
## 100                
## 101                
## 102                
## 103                
## 104                
## 105                
## 106                
## 107                
## 108                
## 109                
## 110                
## 111                
## 112                
## 113                
## 114                
## 115                
## 116                
## 117                
## 118                
## 119                
## 120                
## 121                
## 122                
## 123                
## 124                
## 125                
## 126                
## 127                
## 128                
## 129                
## 130                
## 131                
## 132                
## 133                
## 134                
## 135                
## 136                
## 137                
## 138                
## 139                
## 140                
## 141                
## 142                
## 143                
## 144                
## 145                
## 146                
## 147                
## 148                
## 149                
## 150                
## 151                
## 152                
## 153                
## 154                
## 155                
## 156                
## 157                
## 158                
## 159                
## 160                
## 161                
## 162                
## 163                
## 164                
## 165                
## 166 0              
## 167 0   0          
## 168 0   0   0      
## 169 0   0   0   0  
## Fixed effects: list(fixed) 
##                                                   Value  Std.Error   DF
## (Intercept)                                  -1.3072088 0.04005958 7388
## genotypePos/Pos                              -0.1382761 0.05665281   40
## chowYellow (C8)                              -0.1224580 0.04616532   40
## metabolite_typeAmino Acids                    1.4396406 0.05359902 7388
## metabolite_typeorganic acids                  1.0233612 0.06037107 7388
## genotypePos/Pos:chowYellow (C8)               0.0580403 0.06528762   40
## genotypePos/Pos:metabolite_typeAmino Acids    0.0753990 0.07580046 7388
## genotypePos/Pos:metabolite_typeorganic acids  0.1235755 0.08537759 7388
##                                                t-value p-value
## (Intercept)                                  -32.63161  0.0000
## genotypePos/Pos                               -2.44076  0.0192
## chowYellow (C8)                               -2.65260  0.0114
## metabolite_typeAmino Acids                    26.85946  0.0000
## metabolite_typeorganic acids                  16.95118  0.0000
## genotypePos/Pos:chowYellow (C8)                0.88899  0.3793
## genotypePos/Pos:metabolite_typeAmino Acids     0.99470  0.3199
## genotypePos/Pos:metabolite_typeorganic acids   1.44740  0.1478
##  Correlation: 
##                                              (Intr) gntP/P cY(C8) mtb_AA
## genotypePos/Pos                              -0.707                     
## chowYellow (C8)                              -0.576  0.407              
## metabolite_typeAmino Acids                   -0.499  0.353  0.000       
## metabolite_typeorganic acids                 -0.443  0.313  0.000  0.331
## genotypePos/Pos:chowYellow (C8)               0.407 -0.576 -0.707  0.000
## genotypePos/Pos:metabolite_typeAmino Acids    0.353 -0.499  0.000 -0.707
## genotypePos/Pos:metabolite_typeorganic acids  0.313 -0.443  0.000 -0.234
##                                              mtbl_a gP/P:( gP/P:A
## genotypePos/Pos                                                  
## chowYellow (C8)                                                  
## metabolite_typeAmino Acids                                       
## metabolite_typeorganic acids                                     
## genotypePos/Pos:chowYellow (C8)               0.000              
## genotypePos/Pos:metabolite_typeAmino Acids   -0.234  0.000       
## genotypePos/Pos:metabolite_typeorganic acids -0.707  0.000  0.331
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -5.03533870 -0.47485373  0.05099419  0.50934516  4.29841962 
## 
## Number of Observations: 7436
## Number of Groups: 44
```
