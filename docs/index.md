---
title: "Metabolomics of very long-chain aclCoA dehydrogenase knockout mice"
date: "2016-12-13 12:25:46"
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
L1 <- importDataToList("../data/raw/Ultragenyx Aim 1 Z score.xlsx")
```

```
## Warning in `levels<-`(`*tmp*`, value = if (nl == nL) as.character(labels)
## else paste0(labels, : duplicated levels in factors are deprecated

## Warning in `levels<-`(`*tmp*`, value = if (nl == nL) as.character(labels)
## else paste0(labels, : duplicated levels in factors are deprecated

## Warning in `levels<-`(`*tmp*`, value = if (nl == nL) as.character(labels)
## else paste0(labels, : duplicated levels in factors are deprecated

## Warning in `levels<-`(`*tmp*`, value = if (nl == nL) as.character(labels)
## else paste0(labels, : duplicated levels in factors are deprecated
```

```r
L1[["file"]]
```

```
## [1] "../data/raw/Ultragenyx Aim 1 Z score.xlsx"
```

```r
L1[["dim"]]
```

```
## [1] 622  11
```

```r
L1[["names"]]
```

```
##  [1] "id"              "genotype"        "activity"       
##  [4] "chow"            "metabolite_type" "metabolite"     
##  [7] "value"           "logValue"        "zValue"         
## [10] "zLogValue"       "important"
```

```r
L1[["head"]]
```

```
## # A tibble: 6 × 11
##      id genotype activity    chow metabolite_type metabolite   value
##   <chr>   <fctr>   <fctr>  <fctr>          <fctr>      <chr>   <dbl>
## 1  1170       WT     Rest Regular  Acylcarnitines LCAC total  4.8986
## 2  1171       WT     Rest Regular  Acylcarnitines LCAC total  5.8048
## 3  1172       WT     Rest Regular  Acylcarnitines LCAC total  5.7869
## 4  1201       WT     Rest Regular  Acylcarnitines LCAC total 10.1177
## 5  1202       WT     Rest Regular  Acylcarnitines LCAC total  6.7489
## 6  1209       WT     Rest Regular  Acylcarnitines LCAC total  7.6464
## # ... with 4 more variables: logValue <dbl>, zValue <dbl>,
## #   zLogValue <dbl>, important <lgl>
```

```r
D1 <- L1[["data"]]
```

Import Aim 2: Diet data.
Metabolite `LCAC total` is the sum of `LC even AC total` and `LC odd AC total`; remove `LCAC total` data.


```r
L2 <- importDataToList("../data/raw/Ultragenyx Aim 2 Z score.xlsx")
```

```
## Warning in read_xlsx_(path, sheet, col_names = col_names, col_types =
## col_types, : [4173, 23]: expecting numeric: got ' ctrl mean val'
```

```
## Warning in read_xlsx_(path, sheet, col_names = col_names, col_types =
## col_types, : [4173, 24]: expecting numeric: got 'stdev ctrl val'
```

```
## Warning in read_xlsx_(path, sheet, col_names = col_names, col_types =
## col_types, : [4173, 25]: expecting numeric: got 'mean log ctrl val'
```

```
## Warning in read_xlsx_(path, sheet, col_names = col_names, col_types =
## col_types, : [4173, 26]: expecting numeric: got 'stdev log ctrl val'
```

```
## Warning in `levels<-`(`*tmp*`, value = if (nl == nL) as.character(labels)
## else paste0(labels, : duplicated levels in factors are deprecated

## Warning in `levels<-`(`*tmp*`, value = if (nl == nL) as.character(labels)
## else paste0(labels, : duplicated levels in factors are deprecated

## Warning in `levels<-`(`*tmp*`, value = if (nl == nL) as.character(labels)
## else paste0(labels, : duplicated levels in factors are deprecated

## Warning in `levels<-`(`*tmp*`, value = if (nl == nL) as.character(labels)
## else paste0(labels, : duplicated levels in factors are deprecated
```

```r
L2[["file"]]
```

```
## [1] "../data/raw/Ultragenyx Aim 2 Z score.xlsx"
```

```r
L2[["dim"]]
```

```
## [1] 688  11
```

```r
L2[["names"]]
```

```
##  [1] "id"              "genotype"        "activity"       
##  [4] "chow"            "metabolite_type" "metabolite"     
##  [7] "value"           "logValue"        "zValue"         
## [10] "zLogValue"       "important"
```

```r
L2[["head"]]
```

```
## # A tibble: 6 × 11
##      id genotype activity        chow metabolite_type       metabolite
##   <chr>   <fctr>   <fctr>      <fctr>          <fctr>            <chr>
## 1  1101       WT Exercise Yellow (C8)  Acylcarnitines LC even AC total
## 2  1102       WT Exercise Yellow (C8)  Acylcarnitines LC even AC total
## 3  1103       WT Exercise Yellow (C8)  Acylcarnitines LC even AC total
## 4  1184       WT Exercise Yellow (C8)  Acylcarnitines LC even AC total
## 5  1185       WT Exercise Yellow (C8)  Acylcarnitines LC even AC total
## 6  1186       WT Exercise Yellow (C8)  Acylcarnitines LC even AC total
## # ... with 5 more variables: value <dbl>, logValue <dbl>, zValue <dbl>,
## #   zLogValue <dbl>, important <lgl>
```

```r
D2 <- L2[["data"]] %>% filter(metabolite != "LCAC total")
```


## Check data

Check the `value` and `logValue`.


```r
kable(summarizeOutcome(D1))
```



|              |       Min.| X1st.Qu.| Median|    Mean| X3rd.Qu.|    Max.|
|:-------------|----------:|--------:|------:|-------:|--------:|-------:|
|nominal       |  0.0000001|  1.06900| 3.0360| 28.4200|   12.070| 438.700|
|log-transform | -6.9550000|  0.02907| 0.4823|  0.4017|    1.082|   2.642|

```r
ggplot(D1) +
  aes(x = logValue, color = metabolite_type, fill = metabolite_type) +
  geom_density(alpha = 1/3) +
  scale_x_continuous("log10 scale") +
  facet_grid(genotype ~ activity) +
  scale_color_brewer("Metabolite type", palette = "Set1") +
  scale_fill_brewer("Metabolite type", palette = "Set1") +
  theme_bw()
```

![plot of chunk densitiesAim1](../figures/densitiesAim1-1.png)


```r
kable(summarizeOutcome(D2))
```



|              |       Min.| X1st.Qu.| Median|    Mean| X3rd.Qu.|    Max.|
|:-------------|----------:|--------:|------:|-------:|--------:|-------:|
|nominal       |  0.0000001|  0.84630| 3.2620| 24.7700|   12.070| 408.800|
|log-transform | -6.9550000| -0.07246| 0.5134|  0.4957|    1.082|   2.612|

```r
ggplot(D2) +
  aes(x = logValue, color = metabolite_type, fill = metabolite_type) +
  geom_density(alpha = 1/3) +
  scale_x_continuous("log10 scale") +
  facet_grid(genotype ~ chow) +
  scale_color_brewer("Metabolite type", palette = "Set1") +
  scale_fill_brewer("Metabolite type", palette = "Set1") +
  theme_bw()
```

![plot of chunk densitiesAim2](../figures/densitiesAim2-1.png)


Check fixed effects factors.


```r
D1 %>% group_by(genotype, activity) %>% tally
```

```
## Source: local data frame [4 x 3]
## Groups: genotype [?]
## 
##   genotype activity     n
##     <fctr>   <fctr> <int>
## 1       WT     Rest   142
## 2       WT Exercise   165
## 3       KO     Rest   150
## 4       KO Exercise   165
```

```r
D1 %>% group_by(genotype, metabolite) %>% tally %>% print(n = nrow(.))
```

```
## Source: local data frame [30 x 3]
## Groups: genotype [?]
## 
##    genotype       metabolite     n
##      <fctr>            <chr> <int>
## 1        WT 3-HYDROXYBUTYRIC    20
## 2        WT         arginine    21
## 3        WT           CITRIC    20
## 4        WT          FUMARIC    20
## 5        WT        glutamine    21
## 6        WT       isoleucine    21
## 7        WT           LACTIC    20
## 8        WT       LCAC total    21
## 9        WT          leucine    21
## 10       WT            MALIC    20
## 11       WT       MCAC Total    21
## 12       WT   METHYLSUCCINIC    20
## 13       WT      PYRUVIC_P2P    20
## 14       WT       SUCCINIC-2    20
## 15       WT           valine    21
## 16       KO 3-HYDROXYBUTYRIC    21
## 17       KO         arginine    21
## 18       KO           CITRIC    21
## 19       KO          FUMARIC    21
## 20       KO        glutamine    21
## 21       KO       isoleucine    21
## 22       KO           LACTIC    21
## 23       KO       LCAC total    21
## 24       KO          leucine    21
## 25       KO            MALIC    21
## 26       KO       MCAC Total    21
## 27       KO   METHYLSUCCINIC    21
## 28       KO      PYRUVIC_P2P    21
## 29       KO       SUCCINIC-2    21
## 30       KO           valine    21
```

```r
D1 %>% group_by(activity, metabolite) %>% tally %>% print(n = nrow(.))
```

```
## Source: local data frame [30 x 3]
## Groups: activity [?]
## 
##    activity       metabolite     n
##      <fctr>            <chr> <int>
## 1      Rest 3-HYDROXYBUTYRIC    19
## 2      Rest         arginine    20
## 3      Rest           CITRIC    19
## 4      Rest          FUMARIC    19
## 5      Rest        glutamine    20
## 6      Rest       isoleucine    20
## 7      Rest           LACTIC    19
## 8      Rest       LCAC total    20
## 9      Rest          leucine    20
## 10     Rest            MALIC    19
## 11     Rest       MCAC Total    20
## 12     Rest   METHYLSUCCINIC    19
## 13     Rest      PYRUVIC_P2P    19
## 14     Rest       SUCCINIC-2    19
## 15     Rest           valine    20
## 16 Exercise 3-HYDROXYBUTYRIC    22
## 17 Exercise         arginine    22
## 18 Exercise           CITRIC    22
## 19 Exercise          FUMARIC    22
## 20 Exercise        glutamine    22
## 21 Exercise       isoleucine    22
## 22 Exercise           LACTIC    22
## 23 Exercise       LCAC total    22
## 24 Exercise          leucine    22
## 25 Exercise            MALIC    22
## 26 Exercise       MCAC Total    22
## 27 Exercise   METHYLSUCCINIC    22
## 28 Exercise      PYRUVIC_P2P    22
## 29 Exercise       SUCCINIC-2    22
## 30 Exercise           valine    22
```


```r
D2 %>% group_by(genotype, chow) %>% tally
```

```
## Source: local data frame [4 x 3]
## Groups: genotype [?]
## 
##   genotype        chow     n
##     <fctr>      <fctr> <int>
## 1       WT  White (C7)   165
## 2       WT Yellow (C8)   165
## 3       KO  White (C7)   165
## 4       KO Yellow (C8)   150
```

```r
D2 %>% group_by(genotype, metabolite) %>% tally %>% print(n = nrow(.))
```

```
## Source: local data frame [30 x 3]
## Groups: genotype [?]
## 
##    genotype       metabolite     n
##      <fctr>            <chr> <int>
## 1        WT 3-HYDROXYBUTYRIC    22
## 2        WT         arginine    22
## 3        WT           CITRIC    22
## 4        WT          FUMARIC    22
## 5        WT        glutamine    22
## 6        WT       isoleucine    22
## 7        WT           LACTIC    22
## 8        WT LC even AC total    22
## 9        WT  LC odd AC total    22
## 10       WT          leucine    22
## 11       WT            MALIC    22
## 12       WT       MCAC total    22
## 13       WT   METHYLSUCCINIC    22
## 14       WT       SUCCINIC-2    22
## 15       WT           valine    22
## 16       KO 3-HYDROXYBUTYRIC    21
## 17       KO         arginine    21
## 18       KO           CITRIC    21
## 19       KO          FUMARIC    21
## 20       KO        glutamine    21
## 21       KO       isoleucine    21
## 22       KO           LACTIC    21
## 23       KO LC even AC total    21
## 24       KO  LC odd AC total    21
## 25       KO          leucine    21
## 26       KO            MALIC    21
## 27       KO       MCAC total    21
## 28       KO   METHYLSUCCINIC    21
## 29       KO       SUCCINIC-2    21
## 30       KO           valine    21
```

```r
D2 %>% group_by(chow, metabolite) %>% tally %>% print(n = nrow(.))
```

```
## Source: local data frame [30 x 3]
## Groups: chow [?]
## 
##           chow       metabolite     n
##         <fctr>            <chr> <int>
## 1   White (C7) 3-HYDROXYBUTYRIC    22
## 2   White (C7)         arginine    22
## 3   White (C7)           CITRIC    22
## 4   White (C7)          FUMARIC    22
## 5   White (C7)        glutamine    22
## 6   White (C7)       isoleucine    22
## 7   White (C7)           LACTIC    22
## 8   White (C7) LC even AC total    22
## 9   White (C7)  LC odd AC total    22
## 10  White (C7)          leucine    22
## 11  White (C7)            MALIC    22
## 12  White (C7)       MCAC total    22
## 13  White (C7)   METHYLSUCCINIC    22
## 14  White (C7)       SUCCINIC-2    22
## 15  White (C7)           valine    22
## 16 Yellow (C8) 3-HYDROXYBUTYRIC    21
## 17 Yellow (C8)         arginine    21
## 18 Yellow (C8)           CITRIC    21
## 19 Yellow (C8)          FUMARIC    21
## 20 Yellow (C8)        glutamine    21
## 21 Yellow (C8)       isoleucine    21
## 22 Yellow (C8)           LACTIC    21
## 23 Yellow (C8) LC even AC total    21
## 24 Yellow (C8)  LC odd AC total    21
## 25 Yellow (C8)          leucine    21
## 26 Yellow (C8)            MALIC    21
## 27 Yellow (C8)       MCAC total    21
## 28 Yellow (C8)   METHYLSUCCINIC    21
## 29 Yellow (C8)       SUCCINIC-2    21
## 30 Yellow (C8)           valine    21
```

---

# Model


## Aim 1: Exercise

**Rest versus exhaustive exercise**

Estimate model.
Specify the correlation structure using `cs`.
Use `corSymm`, *general correlation matrix, with no additional structure*.


```r
fixed <- formula(logValue ~
                   genotype +
                   activity +
                   metabolite +
                   genotype * activity +
                   genotype * metabolite +
                   activity * metabolite)
random <- formula(~ 1 | id)
ctrl <- lmeControl(opt = "optim",
                   maxIter = 500, msMaxIter = 500)
cs <-
  corSymm(form = random, fixed = FALSE) %>%
  Initialize(data = D1)
Dim(cs)
```

```
## $N
## [1] 622
## 
## $M
## [1] 42
## 
## $maxLen
## [1] 15
## 
## $sumLenSq
## [1] 9274
## 
## $len
## groups
##   1170   1171   1172   1201   1202   1209   1210   1211   1212 1208B  
##     15     15     15     15      7     15     15     15     15     15 
##   1028   1029   1030   1033   1034   1046   1090   1091   1092   1094 
##     15     15     15     15     15     15     15     15     15     15 
##   1095   1134   1135   1150   1151   1163   1164 1204B  1205B  1206B  
##     15     15     15     15     15     15     15     15     15     15 
## 1207B    1010   1012   1014   1017   1018   1019   1060   1066   1073 
##     15     15     15     15     15     15     15     15     15     15 
##   1076   1077 
##     15     15 
## 
## $start
##  [1]  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
## [24] 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41
```

```r
M <- lme(fixed, data = D1, random = random, correlation = cs, control = ctrl)
plot(M)
```

![plot of chunk lmeDiagnosticAim1](../figures/lmeDiagnosticAim1-1.png)

```r
anova(M)
```

```
##                     numDF denDF  F-value p-value
## (Intercept)             1   538 21296925  <.0001
## genotype                1    38 14177385  <.0001
## activity                1    38  9785253  <.0001
## metabolite             14   538 22936360  <.0001
## genotype:activity       1    38        0  0.6035
## genotype:metabolite    14   538        7  <.0001
## activity:metabolite    14   538        1  0.8918
```

```r
summary(M)
```

```
## Linear mixed-effects model fit by REML
##  Data: D1 
##        AIC      BIC    logLik
##   1694.892 2361.376 -694.4458
## 
## Random effects:
##  Formula: ~1 | id
##         (Intercept)  Residual
## StdDev:  0.06222017 0.7252289
## 
## Correlation Structure: General
##  Formula: ~1 | id 
##  Parameter estimate(s):
##  Correlation: 
##    1      2      3      4      5      6      7      8      9      10    
## 2  -0.086                                                               
## 3   0.025  0.080                                                        
## 4   0.038 -0.043  0.014                                                 
## 5   0.038 -0.052  0.056  0.014                                          
## 6   0.042  0.021  0.027  0.086  0.034                                   
## 7  -0.076  0.162  0.035  0.007 -0.151  0.002                            
## 8   0.020  0.047  0.237  0.169  0.068  0.177  0.031                     
## 9   0.019  0.008 -0.089  0.010 -0.003 -0.077 -0.245  0.025              
## 10 -0.164  0.076  0.031 -0.024  0.016  0.008  0.289  0.145 -0.058       
## 11  0.027 -0.009  0.083  0.060  0.155  0.023  0.017  0.318 -0.052  0.061
## 12  0.108  0.057  0.029  0.111 -0.443  0.074  0.062 -0.145 -0.039 -0.333
## 13 -0.154  0.119  0.053 -0.023 -0.117  0.058  0.265  0.020 -0.050  0.187
## 14  0.100 -0.017  0.074 -0.007  0.052 -0.062  0.000 -0.128 -0.018 -0.058
## 15  0.006  0.021  0.044  0.172 -0.033  0.070  0.036  0.176 -0.011  0.054
##    11     12     13     14    
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
## 12 -0.388                     
## 13 -0.017 -0.006              
## 14  0.061  0.069 -0.060       
## 15  0.026  0.137  0.098  0.010
## Fixed effects: list(fixed) 
##                                                Value Std.Error  DF
## (Intercept)                                0.8435123 0.1993774 538
## genotypeKO                                -0.0096514 0.2070457  38
## activityExercise                           0.1620288 0.2094355  38
## metabolitearginine                        -0.4253956 0.2719001 538
## metaboliteCITRIC                           0.2235985 0.2760201 538
## metaboliteFUMARIC                         -1.7593045 0.2801625 538
## metaboliteglutamine                        1.2263403 0.2739290 538
## metaboliteisoleucine                      -0.5665789 0.2694938 538
## metaboliteLACTIC                           1.6101043 0.2800638 538
## metaboliteLCAC total                       0.0580709 0.2729833 538
## metaboliteleucine                         -0.2426250 0.2711761 538
## metaboliteMALIC                            0.1188729 0.2799638 538
## metaboliteMCAC Total                      -1.2910636 0.2760959 538
## metaboliteMETHYLSUCCINIC                  -2.1973093 0.2813740 538
## metabolitePYRUVIC_P2P                     -1.4600854 0.2748223 538
## metaboliteSUCCINIC-2                       0.3582403 0.2798113 538
## metabolitevaline                          -0.4935056 0.2715735 538
## genotypeKO:activityExercise                0.0743726 0.1205622  38
## genotypeKO:metabolitearginine              0.0531687 0.2755158 538
## genotypeKO:metaboliteCITRIC               -0.0929074 0.2741748 538
## genotypeKO:metaboliteFUMARIC              -0.0496999 0.2889446 538
## genotypeKO:metaboliteglutamine            -0.0404563 0.2829928 538
## genotypeKO:metaboliteisoleucine           -0.1019044 0.2658592 538
## genotypeKO:metaboliteLACTIC               -0.0887323 0.2885991 538
## genotypeKO:metaboliteLCAC total            0.3701918 0.2795042 538
## genotypeKO:metaboliteleucine              -0.4141598 0.2725178 538
## genotypeKO:metaboliteMALIC                -0.1397596 0.2882483 538
## genotypeKO:metaboliteMCAC Total           -0.2073685 0.2910631 538
## genotypeKO:metaboliteMETHYLSUCCINIC       -0.1269491 0.2931634 538
## genotypeKO:metabolitePYRUVIC_P2P          -1.5588673 0.2697961 538
## genotypeKO:metaboliteSUCCINIC-2           -0.1630092 0.2877133 538
## genotypeKO:metabolitevaline                0.0464627 0.2741423 538
## activityExercise:metabolitearginine       -0.3857341 0.2801195 538
## activityExercise:metaboliteCITRIC         -0.1193411 0.2792636 538
## activityExercise:metaboliteFUMARIC        -0.1783012 0.2925823 538
## activityExercise:metaboliteglutamine      -0.2868221 0.2868500 538
## activityExercise:metaboliteisoleucine     -0.2323430 0.2714854 538
## activityExercise:metaboliteLACTIC         -0.1391883 0.2922699 538
## activityExercise:metaboliteLCAC total     -0.3102560 0.2837081 538
## activityExercise:metaboliteleucine        -0.5660946 0.2774401 538
## activityExercise:metaboliteMALIC          -0.1077386 0.2919527 538
## activityExercise:metaboliteMCAC Total     -0.3244149 0.2941258 538
## activityExercise:metaboliteMETHYLSUCCINIC -0.2171200 0.2964005 538
## activityExercise:metabolitePYRUVIC_P2P    -0.3803011 0.2753307 538
## activityExercise:metaboliteSUCCINIC-2     -0.3209592 0.2914690 538
## activityExercise:metabolitevaline         -0.0771800 0.2788925 538
##                                             t-value p-value
## (Intercept)                                4.230732  0.0000
## genotypeKO                                -0.046615  0.9631
## activityExercise                           0.773646  0.4439
## metabolitearginine                        -1.564529  0.1183
## metaboliteCITRIC                           0.810081  0.4183
## metaboliteFUMARIC                         -6.279586  0.0000
## metaboliteglutamine                        4.476855  0.0000
## metaboliteisoleucine                      -2.102382  0.0360
## metaboliteLACTIC                           5.749062  0.0000
## metaboliteLCAC total                       0.212727  0.8316
## metaboliteleucine                         -0.894714  0.3713
## metaboliteMALIC                            0.424601  0.6713
## metaboliteMCAC Total                      -4.676141  0.0000
## metaboliteMETHYLSUCCINIC                  -7.809213  0.0000
## metabolitePYRUVIC_P2P                     -5.312834  0.0000
## metaboliteSUCCINIC-2                       1.280293  0.2010
## metabolitevaline                          -1.817209  0.0697
## genotypeKO:activityExercise                0.616882  0.5410
## genotypeKO:metabolitearginine              0.192979  0.8470
## genotypeKO:metaboliteCITRIC               -0.338862  0.7348
## genotypeKO:metaboliteFUMARIC              -0.172005  0.8635
## genotypeKO:metaboliteglutamine            -0.142959  0.8864
## genotypeKO:metaboliteisoleucine           -0.383302  0.7016
## genotypeKO:metaboliteLACTIC               -0.307459  0.7586
## genotypeKO:metaboliteLCAC total            1.324459  0.1859
## genotypeKO:metaboliteleucine              -1.519753  0.1292
## genotypeKO:metaboliteMALIC                -0.484858  0.6280
## genotypeKO:metaboliteMCAC Total           -0.712452  0.4765
## genotypeKO:metaboliteMETHYLSUCCINIC       -0.433032  0.6652
## genotypeKO:metabolitePYRUVIC_P2P          -5.777947  0.0000
## genotypeKO:metaboliteSUCCINIC-2           -0.566568  0.5712
## genotypeKO:metabolitevaline                0.169484  0.8655
## activityExercise:metabolitearginine       -1.377034  0.1691
## activityExercise:metaboliteCITRIC         -0.427342  0.6693
## activityExercise:metaboliteFUMARIC        -0.609405  0.5425
## activityExercise:metaboliteglutamine      -0.999903  0.3178
## activityExercise:metaboliteisoleucine     -0.855821  0.3925
## activityExercise:metaboliteLACTIC         -0.476232  0.6341
## activityExercise:metaboliteLCAC total     -1.093574  0.2746
## activityExercise:metaboliteleucine        -2.040421  0.0418
## activityExercise:metaboliteMALIC          -0.369027  0.7123
## activityExercise:metaboliteMCAC Total     -1.102980  0.2705
## activityExercise:metaboliteMETHYLSUCCINIC -0.732522  0.4642
## activityExercise:metabolitePYRUVIC_P2P    -1.381252  0.1678
## activityExercise:metaboliteSUCCINIC-2     -1.101178  0.2713
## activityExercise:metabolitevaline         -0.276738  0.7821
##  Correlation: 
##                                           (Intr) gntyKO actvtE mtbltr
## genotypeKO                                -0.546                     
## activityExercise                          -0.578 -0.165              
## metabolitearginine                        -0.693  0.333  0.360       
## metaboliteCITRIC                          -0.684  0.331  0.356  0.501
## metaboliteFUMARIC                         -0.690  0.355  0.379  0.507
## metaboliteglutamine                       -0.699  0.352  0.377  0.500
## metaboliteisoleucine                      -0.700  0.338  0.364  0.502
## metaboliteLACTIC                          -0.687  0.351  0.375  0.495
## metaboliteLCAC total                      -0.697  0.344  0.370  0.492
## metaboliteleucine                         -0.698  0.340  0.366  0.502
## metaboliteMALIC                           -0.690  0.356  0.379  0.494
## metaboliteMCAC Total                      -0.700  0.360  0.384  0.498
## metaboliteMETHYLSUCCINIC                  -0.688  0.356  0.379  0.498
## metabolitePYRUVIC_P2P                     -0.686  0.330  0.356  0.493
## metaboliteSUCCINIC-2                      -0.689  0.353  0.376  0.500
## metabolitevaline                          -0.696  0.337  0.363  0.494
## genotypeKO:activityExercise                0.185 -0.340 -0.321 -0.014
## genotypeKO:metabolitearginine              0.338 -0.619  0.243 -0.519
## genotypeKO:metaboliteCITRIC                0.346 -0.633  0.234 -0.252
## genotypeKO:metaboliteFUMARIC               0.358 -0.655  0.171 -0.265
## genotypeKO:metaboliteglutamine             0.351 -0.643  0.199 -0.235
## genotypeKO:metaboliteisoleucine            0.352 -0.644  0.250 -0.238
## genotypeKO:metaboliteLACTIC                0.354 -0.647  0.179 -0.242
## genotypeKO:metaboliteLCAC total            0.346 -0.634  0.217 -0.220
## genotypeKO:metaboliteleucine               0.348 -0.636  0.236 -0.237
## genotypeKO:metaboliteMALIC                 0.359 -0.656  0.171 -0.240
## genotypeKO:metaboliteMCAC Total            0.353 -0.646  0.173 -0.233
## genotypeKO:metaboliteMETHYLSUCCINIC        0.354 -0.648  0.165 -0.248
## genotypeKO:metabolitePYRUVIC_P2P           0.349 -0.639  0.241 -0.238
## genotypeKO:metaboliteSUCCINIC-2            0.356 -0.651  0.177 -0.252
## genotypeKO:metabolitevaline                0.343 -0.628  0.239 -0.222
## activityExercise:metabolitearginine        0.364  0.242 -0.629 -0.552
## activityExercise:metaboliteCITRIC          0.370  0.232 -0.641 -0.270
## activityExercise:metaboliteFUMARIC         0.381  0.170 -0.660 -0.282
## activityExercise:metaboliteglutamine       0.376  0.198 -0.651 -0.254
## activityExercise:metaboliteisoleucine      0.376  0.247 -0.651 -0.257
## activityExercise:metaboliteLACTIC          0.377  0.179 -0.653 -0.260
## activityExercise:metaboliteLCAC total      0.371  0.216 -0.642 -0.240
## activityExercise:metaboliteleucine         0.373  0.234 -0.645 -0.257
## activityExercise:metaboliteMALIC           0.382  0.171 -0.661 -0.258
## activityExercise:metaboliteMCAC Total      0.377  0.173 -0.653 -0.252
## activityExercise:metaboliteMETHYLSUCCINIC  0.378  0.165 -0.654 -0.266
## activityExercise:metabolitePYRUVIC_P2P     0.374  0.239 -0.647 -0.256
## activityExercise:metaboliteSUCCINIC-2      0.380  0.177 -0.657 -0.270
## activityExercise:metabolitevaline          0.368  0.237 -0.638 -0.242
##                                           mCITRI mFUMAR mtbltg mtblts
## genotypeKO                                                           
## activityExercise                                                     
## metabolitearginine                                                   
## metaboliteCITRIC                                                     
## metaboliteFUMARIC                          0.488                     
## metaboliteglutamine                        0.498  0.497              
## metaboliteisoleucine                       0.502  0.495  0.503       
## metaboliteLACTIC                           0.489  0.490  0.506  0.506
## metaboliteLCAC total                       0.502  0.505  0.502  0.504
## metaboliteleucine                          0.493  0.486  0.502  0.516
## metaboliteMALIC                            0.491  0.501  0.502  0.507
## metaboliteMCAC Total                       0.495  0.507  0.511  0.511
## metaboliteMETHYLSUCCINIC                   0.490  0.501  0.499  0.493
## metabolitePYRUVIC_P2P                      0.496  0.498  0.501  0.506
## metaboliteSUCCINIC-2                       0.491  0.507  0.505  0.499
## metabolitevaline                           0.486  0.499  0.513  0.508
## genotypeKO:activityExercise               -0.008 -0.001 -0.022  0.002
## genotypeKO:metabolitearginine             -0.247 -0.269 -0.238 -0.231
## genotypeKO:metaboliteCITRIC               -0.523 -0.249 -0.252 -0.248
## genotypeKO:metaboliteFUMARIC              -0.240 -0.543 -0.251 -0.236
## genotypeKO:metaboliteglutamine            -0.243 -0.251 -0.529 -0.234
## genotypeKO:metaboliteisoleucine           -0.249 -0.246 -0.245 -0.506
## genotypeKO:metaboliteLACTIC               -0.241 -0.252 -0.268 -0.257
## genotypeKO:metaboliteLCAC total           -0.249 -0.265 -0.243 -0.237
## genotypeKO:metaboliteleucine              -0.233 -0.227 -0.242 -0.259
## genotypeKO:metaboliteMALIC                -0.244 -0.273 -0.260 -0.259
## genotypeKO:metaboliteMCAC Total           -0.238 -0.269 -0.261 -0.251
## genotypeKO:metaboliteMETHYLSUCCINIC       -0.244 -0.273 -0.256 -0.234
## genotypeKO:metabolitePYRUVIC_P2P          -0.253 -0.267 -0.257 -0.256
## genotypeKO:metaboliteSUCCINIC-2           -0.245 -0.285 -0.266 -0.243
## genotypeKO:metabolitevaline               -0.218 -0.253 -0.266 -0.244
## activityExercise:metabolitearginine       -0.265 -0.286 -0.257 -0.251
## activityExercise:metaboliteCITRIC         -0.556 -0.266 -0.270 -0.266
## activityExercise:metaboliteFUMARIC        -0.258 -0.574 -0.269 -0.255
## activityExercise:metaboliteglutamine      -0.261 -0.268 -0.562 -0.254
## activityExercise:metaboliteisoleucine     -0.267 -0.263 -0.263 -0.541
## activityExercise:metaboliteLACTIC         -0.259 -0.269 -0.285 -0.275
## activityExercise:metaboliteLCAC total     -0.267 -0.282 -0.261 -0.256
## activityExercise:metaboliteleucine        -0.252 -0.246 -0.261 -0.278
## activityExercise:metaboliteMALIC          -0.262 -0.289 -0.277 -0.276
## activityExercise:metaboliteMCAC Total     -0.256 -0.286 -0.278 -0.269
## activityExercise:metaboliteMETHYLSUCCINIC -0.262 -0.288 -0.273 -0.253
## activityExercise:metabolitePYRUVIC_P2P    -0.270 -0.283 -0.275 -0.274
## activityExercise:metaboliteSUCCINIC-2     -0.262 -0.300 -0.283 -0.261
## activityExercise:metabolitevaline         -0.238 -0.270 -0.283 -0.263
##                                           mLACTI mLCACt mtbltl mMALIC
## genotypeKO                                                           
## activityExercise                                                     
## metabolitearginine                                                   
## metaboliteCITRIC                                                     
## metaboliteFUMARIC                                                    
## metaboliteglutamine                                                  
## metaboliteisoleucine                                                 
## metaboliteLACTIC                                                     
## metaboliteLCAC total                       0.498                     
## metaboliteleucine                          0.505  0.508              
## metaboliteMALIC                            0.495  0.486  0.500       
## metaboliteMCAC Total                       0.501  0.506  0.508  0.503
## metaboliteMETHYLSUCCINIC                   0.491  0.504  0.502  0.511
## metabolitePYRUVIC_P2P                      0.495  0.505  0.489  0.495
## metaboliteSUCCINIC-2                       0.499  0.503  0.502  0.495
## metabolitevaline                           0.485  0.501  0.512  0.506
## genotypeKO:activityExercise                0.000 -0.020 -0.020 -0.006
## genotypeKO:metabolitearginine             -0.245 -0.222 -0.235 -0.242
## genotypeKO:metaboliteCITRIC               -0.250 -0.257 -0.236 -0.253
## genotypeKO:metaboliteFUMARIC              -0.252 -0.263 -0.223 -0.272
## genotypeKO:metaboliteglutamine            -0.268 -0.241 -0.236 -0.259
## genotypeKO:metaboliteisoleucine           -0.268 -0.245 -0.264 -0.270
## genotypeKO:metaboliteLACTIC               -0.542 -0.250 -0.259 -0.262
## genotypeKO:metaboliteLCAC total           -0.251 -0.525 -0.248 -0.229
## genotypeKO:metaboliteleucine              -0.265 -0.253 -0.515 -0.256
## genotypeKO:metaboliteMALIC                -0.262 -0.229 -0.250 -0.542
## genotypeKO:metaboliteMCAC Total           -0.258 -0.250 -0.250 -0.262
## genotypeKO:metaboliteMETHYLSUCCINIC       -0.254 -0.261 -0.254 -0.292
## genotypeKO:metabolitePYRUVIC_P2P          -0.263 -0.263 -0.228 -0.262
## genotypeKO:metaboliteSUCCINIC-2           -0.270 -0.260 -0.254 -0.261
## genotypeKO:metabolitevaline               -0.225 -0.238 -0.255 -0.267
## activityExercise:metabolitearginine       -0.263 -0.242 -0.255 -0.260
## activityExercise:metaboliteCITRIC         -0.267 -0.274 -0.255 -0.270
## activityExercise:metaboliteFUMARIC        -0.269 -0.280 -0.243 -0.288
## activityExercise:metaboliteglutamine      -0.284 -0.260 -0.255 -0.276
## activityExercise:metaboliteisoleucine     -0.285 -0.264 -0.282 -0.286
## activityExercise:metaboliteLACTIC         -0.574 -0.268 -0.276 -0.278
## activityExercise:metaboliteLCAC total     -0.268 -0.557 -0.267 -0.247
## activityExercise:metaboliteleucine        -0.282 -0.271 -0.549 -0.273
## activityExercise:metaboliteMALIC          -0.278 -0.248 -0.268 -0.574
## activityExercise:metaboliteMCAC Total     -0.275 -0.268 -0.268 -0.278
## activityExercise:metaboliteMETHYLSUCCINIC -0.271 -0.278 -0.272 -0.307
## activityExercise:metabolitePYRUVIC_P2P    -0.279 -0.280 -0.248 -0.278
## activityExercise:metaboliteSUCCINIC-2     -0.286 -0.278 -0.272 -0.278
## activityExercise:metabolitevaline         -0.244 -0.257 -0.274 -0.283
##                                           mMCACT mMETHY mPYRUV mSUCCI
## genotypeKO                                                           
## activityExercise                                                     
## metabolitearginine                                                   
## metaboliteCITRIC                                                     
## metaboliteFUMARIC                                                    
## metaboliteglutamine                                                  
## metaboliteisoleucine                                                 
## metaboliteLACTIC                                                     
## metaboliteLCAC total                                                 
## metaboliteleucine                                                    
## metaboliteMALIC                                                      
## metaboliteMCAC Total                                                 
## metaboliteMETHYLSUCCINIC                   0.493                     
## metabolitePYRUVIC_P2P                      0.501  0.490              
## metaboliteSUCCINIC-2                       0.497  0.486  0.483       
## metabolitevaline                           0.507  0.497  0.493  0.500
## genotypeKO:activityExercise               -0.027 -0.018  0.003 -0.006
## genotypeKO:metabolitearginine             -0.241 -0.254 -0.229 -0.256
## genotypeKO:metaboliteCITRIC               -0.252 -0.256 -0.250 -0.253
## genotypeKO:metaboliteFUMARIC              -0.275 -0.275 -0.255 -0.284
## genotypeKO:metaboliteglutamine            -0.266 -0.258 -0.245 -0.265
## genotypeKO:metaboliteisoleucine           -0.267 -0.245 -0.255 -0.252
## genotypeKO:metaboliteLACTIC               -0.263 -0.257 -0.250 -0.269
## genotypeKO:metaboliteLCAC total           -0.256 -0.266 -0.253 -0.262
## genotypeKO:metaboliteleucine              -0.261 -0.264 -0.222 -0.259
## genotypeKO:metaboliteMALIC                -0.267 -0.295 -0.250 -0.261
## genotypeKO:metaboliteMCAC Total           -0.540 -0.246 -0.245 -0.250
## genotypeKO:metaboliteMETHYLSUCCINIC       -0.249 -0.548 -0.241 -0.245
## genotypeKO:metabolitePYRUVIC_P2P          -0.263 -0.255 -0.517 -0.238
## genotypeKO:metaboliteSUCCINIC-2           -0.257 -0.248 -0.228 -0.541
## genotypeKO:metabolitevaline               -0.259 -0.253 -0.228 -0.256
## activityExercise:metabolitearginine       -0.260 -0.271 -0.249 -0.273
## activityExercise:metaboliteCITRIC         -0.270 -0.273 -0.268 -0.270
## activityExercise:metaboliteFUMARIC        -0.291 -0.291 -0.272 -0.299
## activityExercise:metaboliteglutamine      -0.283 -0.275 -0.263 -0.282
## activityExercise:metaboliteisoleucine     -0.284 -0.263 -0.272 -0.269
## activityExercise:metaboliteLACTIC         -0.280 -0.273 -0.268 -0.285
## activityExercise:metaboliteLCAC total     -0.274 -0.282 -0.271 -0.278
## activityExercise:metaboliteleucine        -0.278 -0.280 -0.241 -0.276
## activityExercise:metaboliteMALIC          -0.284 -0.310 -0.267 -0.277
## activityExercise:metaboliteMCAC Total     -0.571 -0.263 -0.263 -0.268
## activityExercise:metaboliteMETHYLSUCCINIC -0.267 -0.579 -0.259 -0.262
## activityExercise:metabolitePYRUVIC_P2P    -0.280 -0.272 -0.551 -0.256
## activityExercise:metaboliteSUCCINIC-2     -0.274 -0.265 -0.246 -0.573
## activityExercise:metabolitevaline         -0.276 -0.270 -0.248 -0.273
##                                           mtbltv gnKO:E gntypKO:mtbltr
## genotypeKO                                                            
## activityExercise                                                      
## metabolitearginine                                                    
## metaboliteCITRIC                                                      
## metaboliteFUMARIC                                                     
## metaboliteglutamine                                                   
## metaboliteisoleucine                                                  
## metaboliteLACTIC                                                      
## metaboliteLCAC total                                                  
## metaboliteleucine                                                     
## metaboliteMALIC                                                       
## metaboliteMCAC Total                                                  
## metaboliteMETHYLSUCCINIC                                              
## metabolitePYRUVIC_P2P                                                 
## metaboliteSUCCINIC-2                                                  
## metabolitevaline                                                      
## genotypeKO:activityExercise               -0.017                      
## genotypeKO:metabolitearginine             -0.221  0.015               
## genotypeKO:metaboliteCITRIC               -0.224  0.014  0.473        
## genotypeKO:metaboliteFUMARIC              -0.248  0.002  0.496        
## genotypeKO:metaboliteglutamine            -0.260  0.030  0.438        
## genotypeKO:metaboliteisoleucine           -0.249 -0.016  0.444        
## genotypeKO:metaboliteLACTIC               -0.222  0.000  0.452        
## genotypeKO:metaboliteLCAC total           -0.235  0.026  0.409        
## genotypeKO:metaboliteleucine              -0.257  0.026  0.443        
## genotypeKO:metaboliteMALIC                -0.261  0.011  0.447        
## genotypeKO:metaboliteMCAC Total           -0.249  0.039  0.434        
## genotypeKO:metaboliteMETHYLSUCCINIC       -0.246  0.034  0.463        
## genotypeKO:metabolitePYRUVIC_P2P          -0.236 -0.005  0.444        
## genotypeKO:metaboliteSUCCINIC-2           -0.252  0.011  0.472        
## genotypeKO:metabolitevaline               -0.517  0.022  0.413        
## activityExercise:metabolitearginine       -0.241  0.015 -0.299        
## activityExercise:metaboliteCITRIC         -0.243  0.014 -0.177        
## activityExercise:metaboliteFUMARIC        -0.266  0.002 -0.124        
## activityExercise:metaboliteglutamine      -0.278  0.029 -0.192        
## activityExercise:metaboliteisoleucine     -0.268 -0.015 -0.224        
## activityExercise:metaboliteLACTIC         -0.242  0.000 -0.166        
## activityExercise:metaboliteLCAC total     -0.254  0.025 -0.226        
## activityExercise:metaboliteleucine        -0.275  0.025 -0.209        
## activityExercise:metaboliteMALIC          -0.279  0.011 -0.172        
## activityExercise:metaboliteMCAC Total     -0.268  0.037 -0.179        
## activityExercise:metaboliteMETHYLSUCCINIC -0.264  0.032 -0.147        
## activityExercise:metabolitePYRUVIC_P2P    -0.255 -0.005 -0.214        
## activityExercise:metaboliteSUCCINIC-2     -0.269  0.011 -0.149        
## activityExercise:metabolitevaline         -0.551  0.021 -0.234        
##                                           gKO:CI gKO:FU gntypKO:mtbltg
## genotypeKO                                                            
## activityExercise                                                      
## metabolitearginine                                                    
## metaboliteCITRIC                                                      
## metaboliteFUMARIC                                                     
## metaboliteglutamine                                                   
## metaboliteisoleucine                                                  
## metaboliteLACTIC                                                      
## metaboliteLCAC total                                                  
## metaboliteleucine                                                     
## metaboliteMALIC                                                       
## metaboliteMCAC Total                                                  
## metaboliteMETHYLSUCCINIC                                              
## metabolitePYRUVIC_P2P                                                 
## metaboliteSUCCINIC-2                                                  
## metabolitevaline                                                      
## genotypeKO:activityExercise                                           
## genotypeKO:metabolitearginine                                         
## genotypeKO:metaboliteCITRIC                                           
## genotypeKO:metaboliteFUMARIC               0.460                      
## genotypeKO:metaboliteglutamine             0.464  0.462               
## genotypeKO:metaboliteisoleucine            0.476  0.452  0.448        
## genotypeKO:metaboliteLACTIC                0.462  0.465  0.494        
## genotypeKO:metaboliteLCAC total            0.477  0.489  0.446        
## genotypeKO:metaboliteleucine               0.445  0.418  0.444        
## genotypeKO:metaboliteMALIC                 0.467  0.503  0.478        
## genotypeKO:metaboliteMCAC Total            0.455  0.496  0.480        
## genotypeKO:metaboliteMETHYLSUCCINIC        0.468  0.502  0.470        
## genotypeKO:metabolitePYRUVIC_P2P           0.484  0.493  0.474        
## genotypeKO:metaboliteSUCCINIC-2            0.468  0.525  0.490        
## genotypeKO:metabolitevaline                0.418  0.466  0.490        
## activityExercise:metabolitearginine       -0.177 -0.124 -0.191        
## activityExercise:metaboliteCITRIC         -0.309 -0.161 -0.169        
## activityExercise:metaboliteFUMARIC        -0.162 -0.186 -0.141        
## activityExercise:metaboliteglutamine      -0.170 -0.141 -0.235        
## activityExercise:metaboliteisoleucine     -0.196 -0.185 -0.202        
## activityExercise:metaboliteLACTIC         -0.161 -0.127 -0.111        
## activityExercise:metaboliteLCAC total     -0.165 -0.123 -0.175        
## activityExercise:metaboliteleucine        -0.210 -0.203 -0.191        
## activityExercise:metaboliteMALIC          -0.156 -0.092 -0.127        
## activityExercise:metaboliteMCAC Total     -0.162 -0.093 -0.119        
## activityExercise:metaboliteMETHYLSUCCINIC -0.146 -0.083 -0.125        
## activityExercise:metabolitePYRUVIC_P2P    -0.180 -0.139 -0.169        
## activityExercise:metaboliteSUCCINIC-2     -0.156 -0.071 -0.117        
## activityExercise:metabolitevaline         -0.232 -0.155 -0.145        
##                                           gntypKO:mtblts gKO:LA gKO:Lt
## genotypeKO                                                            
## activityExercise                                                      
## metabolitearginine                                                    
## metaboliteCITRIC                                                      
## metaboliteFUMARIC                                                     
## metaboliteglutamine                                                   
## metaboliteisoleucine                                                  
## metaboliteLACTIC                                                      
## metaboliteLCAC total                                                  
## metaboliteleucine                                                     
## metaboliteMALIC                                                       
## metaboliteMCAC Total                                                  
## metaboliteMETHYLSUCCINIC                                              
## metabolitePYRUVIC_P2P                                                 
## metaboliteSUCCINIC-2                                                  
## metabolitevaline                                                      
## genotypeKO:activityExercise                                           
## genotypeKO:metabolitearginine                                         
## genotypeKO:metaboliteCITRIC                                           
## genotypeKO:metaboliteFUMARIC                                          
## genotypeKO:metaboliteglutamine                                        
## genotypeKO:metaboliteisoleucine                                       
## genotypeKO:metaboliteLACTIC                0.495                      
## genotypeKO:metaboliteLCAC total            0.454          0.463       
## genotypeKO:metaboliteleucine               0.500          0.489  0.468
## genotypeKO:metaboliteMALIC                 0.498          0.483  0.422
## genotypeKO:metaboliteMCAC Total            0.482          0.476  0.463
## genotypeKO:metaboliteMETHYLSUCCINIC        0.447          0.469  0.485
## genotypeKO:metabolitePYRUVIC_P2P           0.493          0.485  0.489
## genotypeKO:metaboliteSUCCINIC-2            0.466          0.497  0.483
## genotypeKO:metabolitevaline                0.468          0.415  0.440
## activityExercise:metabolitearginine       -0.225         -0.166 -0.225
## activityExercise:metaboliteCITRIC         -0.197         -0.160 -0.165
## activityExercise:metaboliteFUMARIC        -0.187         -0.127 -0.123
## activityExercise:metaboliteglutamine      -0.204         -0.111 -0.176
## activityExercise:metaboliteisoleucine     -0.389         -0.147 -0.205
## activityExercise:metaboliteLACTIC         -0.148         -0.189 -0.148
## activityExercise:metaboliteLCAC total     -0.206         -0.148 -0.264
## activityExercise:metaboliteleucine        -0.178         -0.138 -0.176
## activityExercise:metaboliteMALIC          -0.146         -0.111 -0.187
## activityExercise:metaboliteMCAC Total     -0.155         -0.113 -0.144
## activityExercise:metaboliteMETHYLSUCCINIC -0.183         -0.115 -0.118
## activityExercise:metabolitePYRUVIC_P2P    -0.191         -0.147 -0.162
## activityExercise:metaboliteSUCCINIC-2     -0.177         -0.099 -0.130
## activityExercise:metabolitevaline         -0.204         -0.204 -0.199
##                                           gntypKO:mtbltl gKO:MA gKO:MT
## genotypeKO                                                            
## activityExercise                                                      
## metabolitearginine                                                    
## metaboliteCITRIC                                                      
## metaboliteFUMARIC                                                     
## metaboliteglutamine                                                   
## metaboliteisoleucine                                                  
## metaboliteLACTIC                                                      
## metaboliteLCAC total                                                  
## metaboliteleucine                                                     
## metaboliteMALIC                                                       
## metaboliteMCAC Total                                                  
## metaboliteMETHYLSUCCINIC                                              
## metabolitePYRUVIC_P2P                                                 
## metaboliteSUCCINIC-2                                                  
## metabolitevaline                                                      
## genotypeKO:activityExercise                                           
## genotypeKO:metabolitearginine                                         
## genotypeKO:metaboliteCITRIC                                           
## genotypeKO:metaboliteFUMARIC                                          
## genotypeKO:metaboliteglutamine                                        
## genotypeKO:metaboliteisoleucine                                       
## genotypeKO:metaboliteLACTIC                                           
## genotypeKO:metaboliteLCAC total                                       
## genotypeKO:metaboliteleucine                                          
## genotypeKO:metaboliteMALIC                 0.472                      
## genotypeKO:metaboliteMCAC Total            0.471          0.483       
## genotypeKO:metaboliteMETHYLSUCCINIC        0.481          0.539  0.449
## genotypeKO:metabolitePYRUVIC_P2P           0.429          0.483  0.475
## genotypeKO:metaboliteSUCCINIC-2            0.479          0.482  0.463
## genotypeKO:metabolitevaline                0.483          0.492  0.468
## activityExercise:metabolitearginine       -0.209         -0.171 -0.177
## activityExercise:metaboliteCITRIC         -0.210         -0.155 -0.161
## activityExercise:metaboliteFUMARIC        -0.204         -0.092 -0.092
## activityExercise:metaboliteglutamine      -0.192         -0.127 -0.119
## activityExercise:metaboliteisoleucine     -0.178         -0.144 -0.153
## activityExercise:metaboliteLACTIC         -0.138         -0.111 -0.113
## activityExercise:metaboliteLCAC total     -0.177         -0.187 -0.143
## activityExercise:metaboliteleucine        -0.326         -0.155 -0.150
## activityExercise:metaboliteMALIC          -0.155         -0.192 -0.107
## activityExercise:metaboliteMCAC Total     -0.151         -0.107 -0.171
## activityExercise:metaboliteMETHYLSUCCINIC -0.137         -0.050 -0.130
## activityExercise:metabolitePYRUVIC_P2P    -0.235         -0.150 -0.151
## activityExercise:metaboliteSUCCINIC-2     -0.149         -0.114 -0.126
## activityExercise:metabolitevaline         -0.175         -0.132 -0.149
##                                           gKO:ME gKO:PY gKO:SU
## genotypeKO                                                    
## activityExercise                                              
## metabolitearginine                                            
## metaboliteCITRIC                                              
## metaboliteFUMARIC                                             
## metaboliteglutamine                                           
## metaboliteisoleucine                                          
## metaboliteLACTIC                                              
## metaboliteLCAC total                                          
## metaboliteleucine                                             
## metaboliteMALIC                                               
## metaboliteMCAC Total                                          
## metaboliteMETHYLSUCCINIC                                      
## metabolitePYRUVIC_P2P                                         
## metaboliteSUCCINIC-2                                          
## metabolitevaline                                              
## genotypeKO:activityExercise                                   
## genotypeKO:metabolitearginine                                 
## genotypeKO:metaboliteCITRIC                                   
## genotypeKO:metaboliteFUMARIC                                  
## genotypeKO:metaboliteglutamine                                
## genotypeKO:metaboliteisoleucine                               
## genotypeKO:metaboliteLACTIC                                   
## genotypeKO:metaboliteLCAC total                               
## genotypeKO:metaboliteleucine                                  
## genotypeKO:metaboliteMALIC                                    
## genotypeKO:metaboliteMCAC Total                               
## genotypeKO:metaboliteMETHYLSUCCINIC                           
## genotypeKO:metabolitePYRUVIC_P2P           0.466              
## genotypeKO:metaboliteSUCCINIC-2            0.453  0.440       
## genotypeKO:metabolitevaline                0.462  0.442  0.473
## activityExercise:metabolitearginine       -0.146 -0.215 -0.149
## activityExercise:metaboliteCITRIC         -0.145 -0.180 -0.155
## activityExercise:metaboliteFUMARIC        -0.083 -0.140 -0.071
## activityExercise:metaboliteglutamine      -0.125 -0.171 -0.117
## activityExercise:metaboliteisoleucine     -0.182 -0.191 -0.176
## activityExercise:metaboliteLACTIC         -0.115 -0.148 -0.099
## activityExercise:metaboliteLCAC total     -0.118 -0.163 -0.130
## activityExercise:metaboliteleucine        -0.136 -0.235 -0.149
## activityExercise:metaboliteMALIC          -0.050 -0.151 -0.114
## activityExercise:metaboliteMCAC Total     -0.130 -0.153 -0.127
## activityExercise:metaboliteMETHYLSUCCINIC -0.154 -0.157 -0.132
## activityExercise:metabolitePYRUVIC_P2P    -0.156 -0.350 -0.191
## activityExercise:metaboliteSUCCINIC-2     -0.132 -0.192 -0.196
## activityExercise:metabolitevaline         -0.150 -0.220 -0.151
##                                           gntypKO:mtbltv
## genotypeKO                                              
## activityExercise                                        
## metabolitearginine                                      
## metaboliteCITRIC                                        
## metaboliteFUMARIC                                       
## metaboliteglutamine                                     
## metaboliteisoleucine                                    
## metaboliteLACTIC                                        
## metaboliteLCAC total                                    
## metaboliteleucine                                       
## metaboliteMALIC                                         
## metaboliteMCAC Total                                    
## metaboliteMETHYLSUCCINIC                                
## metabolitePYRUVIC_P2P                                   
## metaboliteSUCCINIC-2                                    
## metabolitevaline                                        
## genotypeKO:activityExercise                             
## genotypeKO:metabolitearginine                           
## genotypeKO:metaboliteCITRIC                             
## genotypeKO:metaboliteFUMARIC                            
## genotypeKO:metaboliteglutamine                          
## genotypeKO:metaboliteisoleucine                         
## genotypeKO:metaboliteLACTIC                             
## genotypeKO:metaboliteLCAC total                         
## genotypeKO:metaboliteleucine                            
## genotypeKO:metaboliteMALIC                              
## genotypeKO:metaboliteMCAC Total                         
## genotypeKO:metaboliteMETHYLSUCCINIC                     
## genotypeKO:metabolitePYRUVIC_P2P                        
## genotypeKO:metaboliteSUCCINIC-2                         
## genotypeKO:metabolitevaline                             
## activityExercise:metabolitearginine       -0.234        
## activityExercise:metaboliteCITRIC         -0.232        
## activityExercise:metaboliteFUMARIC        -0.156        
## activityExercise:metaboliteglutamine      -0.146        
## activityExercise:metaboliteisoleucine     -0.204        
## activityExercise:metaboliteLACTIC         -0.205        
## activityExercise:metaboliteLCAC total     -0.200        
## activityExercise:metaboliteleucine        -0.175        
## activityExercise:metaboliteMALIC          -0.133        
## activityExercise:metaboliteMCAC Total     -0.150        
## activityExercise:metaboliteMETHYLSUCCINIC -0.151        
## activityExercise:metabolitePYRUVIC_P2P    -0.219        
## activityExercise:metaboliteSUCCINIC-2     -0.152        
## activityExercise:metabolitevaline         -0.311        
##                                           actvtyExrcs:mtbltr aE:CIT aE:FUM
## genotypeKO                                                                
## activityExercise                                                          
## metabolitearginine                                                        
## metaboliteCITRIC                                                          
## metaboliteFUMARIC                                                         
## metaboliteglutamine                                                       
## metaboliteisoleucine                                                      
## metaboliteLACTIC                                                          
## metaboliteLCAC total                                                      
## metaboliteleucine                                                         
## metaboliteMALIC                                                           
## metaboliteMCAC Total                                                      
## metaboliteMETHYLSUCCINIC                                                  
## metabolitePYRUVIC_P2P                                                     
## metaboliteSUCCINIC-2                                                      
## metabolitevaline                                                          
## genotypeKO:activityExercise                                               
## genotypeKO:metabolitearginine                                             
## genotypeKO:metaboliteCITRIC                                               
## genotypeKO:metaboliteFUMARIC                                              
## genotypeKO:metaboliteglutamine                                            
## genotypeKO:metaboliteisoleucine                                           
## genotypeKO:metaboliteLACTIC                                               
## genotypeKO:metaboliteLCAC total                                           
## genotypeKO:metaboliteleucine                                              
## genotypeKO:metaboliteMALIC                                                
## genotypeKO:metaboliteMCAC Total                                           
## genotypeKO:metaboliteMETHYLSUCCINIC                                       
## genotypeKO:metabolitePYRUVIC_P2P                                          
## genotypeKO:metaboliteSUCCINIC-2                                           
## genotypeKO:metabolitevaline                                               
## activityExercise:metabolitearginine                                       
## activityExercise:metaboliteCITRIC          0.477                          
## activityExercise:metaboliteFUMARIC         0.497              0.464       
## activityExercise:metaboliteglutamine       0.446              0.469  0.466
## activityExercise:metaboliteisoleucine      0.452              0.479  0.458
## activityExercise:metaboliteLACTIC          0.458              0.466  0.469
## activityExercise:metaboliteLCAC total      0.421              0.480  0.491
## activityExercise:metaboliteleucine         0.451              0.452  0.428
## activityExercise:metaboliteMALIC           0.454              0.471  0.502
## activityExercise:metaboliteMCAC Total      0.443              0.461  0.497
## activityExercise:metaboliteMETHYLSUCCINIC  0.468              0.471  0.502
## activityExercise:metabolitePYRUVIC_P2P     0.451              0.486  0.493
## activityExercise:metaboliteSUCCINIC-2      0.476              0.472  0.523
## activityExercise:metabolitevaline          0.424              0.428  0.470
##                                           actvtyExrcs:mtbltg
## genotypeKO                                                  
## activityExercise                                            
## metabolitearginine                                          
## metaboliteCITRIC                                            
## metaboliteFUMARIC                                           
## metaboliteglutamine                                         
## metaboliteisoleucine                                        
## metaboliteLACTIC                                            
## metaboliteLCAC total                                        
## metaboliteleucine                                           
## metaboliteMALIC                                             
## metaboliteMCAC Total                                        
## metaboliteMETHYLSUCCINIC                                    
## metabolitePYRUVIC_P2P                                       
## metaboliteSUCCINIC-2                                        
## metabolitevaline                                            
## genotypeKO:activityExercise                                 
## genotypeKO:metabolitearginine                               
## genotypeKO:metaboliteCITRIC                                 
## genotypeKO:metaboliteFUMARIC                                
## genotypeKO:metaboliteglutamine                              
## genotypeKO:metaboliteisoleucine                             
## genotypeKO:metaboliteLACTIC                                 
## genotypeKO:metaboliteLCAC total                             
## genotypeKO:metaboliteleucine                                
## genotypeKO:metaboliteMALIC                                  
## genotypeKO:metaboliteMCAC Total                             
## genotypeKO:metaboliteMETHYLSUCCINIC                         
## genotypeKO:metabolitePYRUVIC_P2P                            
## genotypeKO:metaboliteSUCCINIC-2                             
## genotypeKO:metabolitevaline                                 
## activityExercise:metabolitearginine                         
## activityExercise:metaboliteCITRIC                           
## activityExercise:metaboliteFUMARIC                          
## activityExercise:metaboliteglutamine                        
## activityExercise:metaboliteisoleucine      0.456            
## activityExercise:metaboliteLACTIC          0.495            
## activityExercise:metaboliteLCAC total      0.453            
## activityExercise:metaboliteleucine         0.452            
## activityExercise:metaboliteMALIC           0.481            
## activityExercise:metaboliteMCAC Total      0.484            
## activityExercise:metaboliteMETHYLSUCCINIC  0.474            
## activityExercise:metabolitePYRUVIC_P2P     0.477            
## activityExercise:metaboliteSUCCINIC-2      0.492            
## activityExercise:metabolitevaline          0.492            
##                                           actvtyExrcs:mtblts aE:LAC aE:LCt
## genotypeKO                                                                
## activityExercise                                                          
## metabolitearginine                                                        
## metaboliteCITRIC                                                          
## metaboliteFUMARIC                                                         
## metaboliteglutamine                                                       
## metaboliteisoleucine                                                      
## metaboliteLACTIC                                                          
## metaboliteLCAC total                                                      
## metaboliteleucine                                                         
## metaboliteMALIC                                                           
## metaboliteMCAC Total                                                      
## metaboliteMETHYLSUCCINIC                                                  
## metabolitePYRUVIC_P2P                                                     
## metaboliteSUCCINIC-2                                                      
## metabolitevaline                                                          
## genotypeKO:activityExercise                                               
## genotypeKO:metabolitearginine                                             
## genotypeKO:metaboliteCITRIC                                               
## genotypeKO:metaboliteFUMARIC                                              
## genotypeKO:metaboliteglutamine                                            
## genotypeKO:metaboliteisoleucine                                           
## genotypeKO:metaboliteLACTIC                                               
## genotypeKO:metaboliteLCAC total                                           
## genotypeKO:metaboliteleucine                                              
## genotypeKO:metaboliteMALIC                                                
## genotypeKO:metaboliteMCAC Total                                           
## genotypeKO:metaboliteMETHYLSUCCINIC                                       
## genotypeKO:metabolitePYRUVIC_P2P                                          
## genotypeKO:metaboliteSUCCINIC-2                                           
## genotypeKO:metabolitevaline                                               
## activityExercise:metabolitearginine                                       
## activityExercise:metaboliteCITRIC                                         
## activityExercise:metaboliteFUMARIC                                        
## activityExercise:metaboliteglutamine                                      
## activityExercise:metaboliteisoleucine                                     
## activityExercise:metaboliteLACTIC          0.496                          
## activityExercise:metaboliteLCAC total      0.460              0.467       
## activityExercise:metaboliteleucine         0.501              0.491  0.473
## activityExercise:metaboliteMALIC           0.499              0.485  0.431
## activityExercise:metaboliteMCAC Total      0.485              0.479  0.468
## activityExercise:metaboliteMETHYLSUCCINIC  0.454              0.472  0.487
## activityExercise:metabolitePYRUVIC_P2P     0.494              0.486  0.491
## activityExercise:metaboliteSUCCINIC-2      0.470              0.498  0.486
## activityExercise:metabolitevaline          0.473              0.425  0.448
##                                           actvtyExrcs:mtbltl aE:MAL aE:MCT
## genotypeKO                                                                
## activityExercise                                                          
## metabolitearginine                                                        
## metaboliteCITRIC                                                          
## metaboliteFUMARIC                                                         
## metaboliteglutamine                                                       
## metaboliteisoleucine                                                      
## metaboliteLACTIC                                                          
## metaboliteLCAC total                                                      
## metaboliteleucine                                                         
## metaboliteMALIC                                                           
## metaboliteMCAC Total                                                      
## metaboliteMETHYLSUCCINIC                                                  
## metabolitePYRUVIC_P2P                                                     
## metaboliteSUCCINIC-2                                                      
## metabolitevaline                                                          
## genotypeKO:activityExercise                                               
## genotypeKO:metabolitearginine                                             
## genotypeKO:metaboliteCITRIC                                               
## genotypeKO:metaboliteFUMARIC                                              
## genotypeKO:metaboliteglutamine                                            
## genotypeKO:metaboliteisoleucine                                           
## genotypeKO:metaboliteLACTIC                                               
## genotypeKO:metaboliteLCAC total                                           
## genotypeKO:metaboliteleucine                                              
## genotypeKO:metaboliteMALIC                                                
## genotypeKO:metaboliteMCAC Total                                           
## genotypeKO:metaboliteMETHYLSUCCINIC                                       
## genotypeKO:metabolitePYRUVIC_P2P                                          
## genotypeKO:metaboliteSUCCINIC-2                                           
## genotypeKO:metabolitevaline                                               
## activityExercise:metabolitearginine                                       
## activityExercise:metaboliteCITRIC                                         
## activityExercise:metaboliteFUMARIC                                        
## activityExercise:metaboliteglutamine                                      
## activityExercise:metaboliteisoleucine                                     
## activityExercise:metaboliteLACTIC                                         
## activityExercise:metaboliteLCAC total                                     
## activityExercise:metaboliteleucine                                        
## activityExercise:metaboliteMALIC           0.475                          
## activityExercise:metaboliteMCAC Total      0.475              0.485       
## activityExercise:metaboliteMETHYLSUCCINIC  0.483              0.535  0.455
## activityExercise:metabolitePYRUVIC_P2P     0.438              0.485  0.478
## activityExercise:metaboliteSUCCINIC-2      0.482              0.484  0.467
## activityExercise:metabolitevaline          0.486              0.493  0.473
##                                           aE:MET aE:PYR aE:SUC
## genotypeKO                                                    
## activityExercise                                              
## metabolitearginine                                            
## metaboliteCITRIC                                              
## metaboliteFUMARIC                                             
## metaboliteglutamine                                           
## metaboliteisoleucine                                          
## metaboliteLACTIC                                              
## metaboliteLCAC total                                          
## metaboliteleucine                                             
## metaboliteMALIC                                               
## metaboliteMCAC Total                                          
## metaboliteMETHYLSUCCINIC                                      
## metabolitePYRUVIC_P2P                                         
## metaboliteSUCCINIC-2                                          
## metabolitevaline                                              
## genotypeKO:activityExercise                                   
## genotypeKO:metabolitearginine                                 
## genotypeKO:metaboliteCITRIC                                   
## genotypeKO:metaboliteFUMARIC                                  
## genotypeKO:metaboliteglutamine                                
## genotypeKO:metaboliteisoleucine                               
## genotypeKO:metaboliteLACTIC                                   
## genotypeKO:metaboliteLCAC total                               
## genotypeKO:metaboliteleucine                                  
## genotypeKO:metaboliteMALIC                                    
## genotypeKO:metaboliteMCAC Total                               
## genotypeKO:metaboliteMETHYLSUCCINIC                           
## genotypeKO:metabolitePYRUVIC_P2P                              
## genotypeKO:metaboliteSUCCINIC-2                               
## genotypeKO:metabolitevaline                                   
## activityExercise:metabolitearginine                           
## activityExercise:metaboliteCITRIC                             
## activityExercise:metaboliteFUMARIC                            
## activityExercise:metaboliteglutamine                          
## activityExercise:metaboliteisoleucine                         
## activityExercise:metaboliteLACTIC                             
## activityExercise:metaboliteLCAC total                         
## activityExercise:metaboliteleucine                            
## activityExercise:metaboliteMALIC                              
## activityExercise:metaboliteMCAC Total                         
## activityExercise:metaboliteMETHYLSUCCINIC                     
## activityExercise:metabolitePYRUVIC_P2P     0.469              
## activityExercise:metaboliteSUCCINIC-2      0.458  0.447       
## activityExercise:metabolitevaline          0.467  0.449  0.476
## 
## Standardized Within-Group Residuals:
##        Min         Q1        Med         Q3        Max 
## -9.3276555 -0.1084469  0.0128430  0.1374870  3.4032215 
## 
## Number of Observations: 622
## Number of Groups: 42
```


## Aim 2: Chow

**MCT chow versus experimental chow with triheptanoin**

Estimate model.
Specify the correlation structure using `cs`.
Use `corSymm`, *general correlation matrix, with no additional structure*.


```r
rm(M)
fixed <- formula(logValue ~
                   genotype +
                   chow +
                   metabolite +
                   genotype * chow +
                   genotype * metabolite +
                   chow * metabolite)
cs <-
  corSymm(form = random, fixed = FALSE) %>%
  Initialize(data = D2)
Dim(cs)
```

```
## $N
## [1] 645
## 
## $M
## [1] 43
## 
## $maxLen
## [1] 15
## 
## $sumLenSq
## [1] 9675
## 
## $len
## groups
##       1101       1102       1103       1184       1185       1186 
##         15         15         15         15         15         15 
##       1195       1196       1197       1203 1192/1198B       1176 
##         15         15         15         15         15         15 
##       1177       1179       1180       1190       1191       1193 
##         15         15         15         15         15         15 
##       1194       1199       1200 1192A/1198       1107       1113 
##         15         15         15         15         15         15 
##       1114       1115       1117       1118       1119       1204 
##         15         15         15         15         15         15 
##       1205       1206       1120       1126       1127       1128 
##         15         15         15         15         15         15 
##       1142       1143       1144       1145       1146       1158 
##         15         15         15         15         15         15 
##       1159 
##         15 
## 
## $start
##  [1]  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
## [24] 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42
```

```r
M <- lme(fixed, data = D2, random = random, correlation = cs, control = ctrl)
plot(M)
```

![plot of chunk lmeDiagnosticAim2](../figures/lmeDiagnosticAim2-1.png)

```r
anova(M)
```

```
##                     numDF denDF      F-value p-value
## (Intercept)             1   560 1.495067e+06  <.0001
## genotype                1    39 2.360281e+13  <.0001
## chow                    1    39 3.881290e+11  <.0001
## metabolite             14   560 6.164543e+12  <.0001
## genotype:chow           1    39 0.000000e+00  0.5853
## genotype:metabolite    14   560 6.000000e+00  <.0001
## chow:metabolite        14   560 3.000000e+00  0.0001
```

```r
summary(M)
```

```
## Linear mixed-effects model fit by REML
##  Data: D2 
##        AIC      BIC    logLik
##   1377.217 2049.692 -535.6084
## 
## Random effects:
##  Formula: ~1 | id
##         (Intercept)  Residual
## StdDev:   0.1011938 0.5108908
## 
## Correlation Structure: General
##  Formula: ~1 | id 
##  Parameter estimate(s):
##  Correlation: 
##    1      2      3      4      5      6      7      8      9      10    
## 2   0.047                                                               
## 3   0.152  0.075                                                        
## 4   0.063  0.029 -0.032                                                 
## 5  -0.021  0.073  0.071 -0.016                                          
## 6   0.054  0.014  0.012 -0.090  0.032                                   
## 7   0.138  0.139  0.055 -0.045  0.021  0.004                            
## 8   0.082  0.076  0.077 -0.054  0.013  0.002  0.010                     
## 9  -0.046  0.039  0.051 -0.004  0.047 -0.034  0.013 -0.047              
## 10 -0.004 -0.050 -0.007 -0.253  0.043 -0.219 -0.025 -0.030  0.083       
## 11 -0.078 -0.014 -0.001 -0.017  0.025  0.060  0.003  0.018  0.021  0.053
## 12  0.070  0.083  0.028 -0.044 -0.001 -0.033  0.060  0.038  0.043  0.010
## 13 -0.103 -0.017  0.001 -0.285  0.006 -0.176 -0.003 -0.061  0.006 -0.050
## 14 -0.247 -0.174 -0.058 -0.039 -0.010  0.044 -0.127 -0.025 -0.033  0.137
## 15  0.318  0.238  0.180  0.097  0.128 -0.013  0.112  0.287  0.060 -0.215
##    11     12     13     14    
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
## 12  0.042                     
## 13 -0.016 -0.079              
## 14  0.379  0.038 -0.154       
## 15  0.023  0.048 -0.060 -0.290
## Fixed effects: list(fixed) 
##                                                 Value Std.Error  DF
## (Intercept)                                 0.3804304 0.1338167 560
## genotypeKO                                  0.4132004 0.1493635  39
## chowYellow (C8)                             0.2559378 0.1466806  39
## metabolitearginine                         -0.0866796 0.1837287 560
## metaboliteCITRIC                           -0.0210452 0.1818259 560
## metaboliteFUMARIC                          -1.1859464 0.1832035 560
## metaboliteglutamine                         1.6520904 0.1820680 560
## metaboliteisoleucine                       -0.2258870 0.1834171 560
## metaboliteLACTIC                            2.0540483 0.1825703 560
## metaboliteLC even AC total                  0.2802441 0.1828464 560
## metaboliteLC odd AC total                  -0.5825606 0.1847896 560
## metaboliteleucine                          -0.1073972 0.1838510 560
## metaboliteMALIC                             0.6388126 0.1817697 560
## metaboliteMCAC total                       -0.8465304 0.1850690 560
## metaboliteMETHYLSUCCINIC                   -1.6865104 0.1847424 560
## metaboliteSUCCINIC-2                        0.8553226 0.1834661 560
## metabolitevaline                            0.0382429 0.1815215 560
## genotypeKO:chowYellow (C8)                  0.0546717 0.0963259  39
## genotypeKO:metabolitearginine              -0.2529326 0.2002809 560
## genotypeKO:metaboliteCITRIC                 0.1301356 0.1932101 560
## genotypeKO:metaboliteFUMARIC               -0.6107086 0.1983471 560
## genotypeKO:metaboliteglutamine             -0.5019595 0.1941198 560
## genotypeKO:metaboliteisoleucine            -0.3282557 0.1991353 560
## genotypeKO:metaboliteLACTIC                -0.6285505 0.1959977 560
## genotypeKO:metaboliteLC even AC total       0.1445959 0.1970245 560
## genotypeKO:metaboliteLC odd AC total        0.1196305 0.2041478 560
## genotypeKO:metaboliteleucine               -0.2930849 0.2007294 560
## genotypeKO:metaboliteMALIC                 -0.5853944 0.1929982 560
## genotypeKO:metaboliteMCAC total            -0.2831915 0.2051575 560
## genotypeKO:metaboliteMETHYLSUCCINIC        -0.5405962 0.2039766 560
## genotypeKO:metaboliteSUCCINIC-2            -0.5507818 0.1993155 560
## genotypeKO:metabolitevaline                -0.4998860 0.1920618 560
## chowYellow (C8):metabolitearginine         -0.5971881 0.1984089 560
## chowYellow (C8):metaboliteCITRIC           -0.1146199 0.1905589 560
## chowYellow (C8):metaboliteFUMARIC          -0.3666861 0.1935204 560
## chowYellow (C8):metaboliteglutamine        -0.2605192 0.1925106 560
## chowYellow (C8):metaboliteisoleucine       -0.4832110 0.1995500 560
## chowYellow (C8):metaboliteLACTIC           -0.3800258 0.1934636 560
## chowYellow (C8):metaboliteLC even AC total -0.2521383 0.1961389 560
## chowYellow (C8):metaboliteLC odd AC total  -0.7688645 0.1985949 560
## chowYellow (C8):metaboliteleucine          -0.4609447 0.2005267 560
## chowYellow (C8):metaboliteMALIC            -0.4063418 0.1910524 560
## chowYellow (C8):metaboliteMCAC total       -0.3194601 0.1987676 560
## chowYellow (C8):metaboliteMETHYLSUCCINIC    0.2928661 0.2017272 560
## chowYellow (C8):metaboliteSUCCINIC-2       -0.3319469 0.1967629 560
## chowYellow (C8):metabolitevaline           -0.0846998 0.1907487 560
##                                              t-value p-value
## (Intercept)                                 2.842923  0.0046
## genotypeKO                                  2.766407  0.0086
## chowYellow (C8)                             1.744864  0.0889
## metabolitearginine                         -0.471780  0.6373
## metaboliteCITRIC                           -0.115743  0.9079
## metaboliteFUMARIC                          -6.473383  0.0000
## metaboliteglutamine                         9.074030  0.0000
## metaboliteisoleucine                       -1.231548  0.2186
## metaboliteLACTIC                           11.250726  0.0000
## metaboliteLC even AC total                  1.532675  0.1259
## metaboliteLC odd AC total                  -3.152561  0.0017
## metaboliteleucine                          -0.584153  0.5594
## metaboliteMALIC                             3.514408  0.0005
## metaboliteMCAC total                       -4.574134  0.0000
## metaboliteMETHYLSUCCINIC                   -9.128986  0.0000
## metaboliteSUCCINIC-2                        4.662021  0.0000
## metabolitevaline                            0.210679  0.8332
## genotypeKO:chowYellow (C8)                  0.567570  0.5736
## genotypeKO:metabolitearginine              -1.262890  0.2072
## genotypeKO:metaboliteCITRIC                 0.673544  0.5009
## genotypeKO:metaboliteFUMARIC               -3.078989  0.0022
## genotypeKO:metaboliteglutamine             -2.585824  0.0100
## genotypeKO:metaboliteisoleucine            -1.648406  0.0998
## genotypeKO:metaboliteLACTIC                -3.206927  0.0014
## genotypeKO:metaboliteLC even AC total       0.733898  0.4633
## genotypeKO:metaboliteLC odd AC total        0.586000  0.5581
## genotypeKO:metaboliteleucine               -1.460100  0.1448
## genotypeKO:metaboliteMALIC                 -3.033160  0.0025
## genotypeKO:metaboliteMCAC total            -1.380361  0.1680
## genotypeKO:metaboliteMETHYLSUCCINIC        -2.650286  0.0083
## genotypeKO:metaboliteSUCCINIC-2            -2.763366  0.0059
## genotypeKO:metabolitevaline                -2.602734  0.0095
## chowYellow (C8):metabolitearginine         -3.009885  0.0027
## chowYellow (C8):metaboliteCITRIC           -0.601493  0.5478
## chowYellow (C8):metaboliteFUMARIC          -1.894819  0.0586
## chowYellow (C8):metaboliteglutamine        -1.353272  0.1765
## chowYellow (C8):metaboliteisoleucine       -2.421504  0.0158
## chowYellow (C8):metaboliteLACTIC           -1.964327  0.0500
## chowYellow (C8):metaboliteLC even AC total -1.285509  0.1991
## chowYellow (C8):metaboliteLC odd AC total  -3.871523  0.0001
## chowYellow (C8):metaboliteleucine          -2.298670  0.0219
## chowYellow (C8):metaboliteMALIC            -2.126861  0.0339
## chowYellow (C8):metaboliteMCAC total       -1.607204  0.1086
## chowYellow (C8):metaboliteMETHYLSUCCINIC    1.451793  0.1471
## chowYellow (C8):metaboliteSUCCINIC-2       -1.687040  0.0922
## chowYellow (C8):metabolitevaline           -0.444039  0.6572
##  Correlation: 
##                                            (Intr) gntyKO chY(C8) mtbltr
## genotypeKO                                 -0.558                      
## chowYellow (C8)                            -0.829  0.359               
## metabolitearginine                         -0.676  0.346  0.521        
## metaboliteCITRIC                           -0.667  0.321  0.563   0.477
## metaboliteFUMARIC                          -0.676  0.343  0.553   0.484
## metaboliteglutamine                        -0.669  0.326  0.573   0.498
## metaboliteisoleucine                       -0.679  0.350  0.525   0.499
## metaboliteLACTIC                           -0.676  0.340  0.554   0.500
## metaboliteLC even AC total                 -0.670  0.332  0.567   0.490
## metaboliteLC odd AC total                  -0.684  0.366  0.512   0.509
## metaboliteleucine                          -0.666  0.329  0.545   0.481
## metaboliteMALIC                            -0.670  0.327  0.567   0.492
## metaboliteMCAC total                       -0.677  0.354  0.532   0.495
## metaboliteMETHYLSUCCINIC                   -0.672  0.344  0.537   0.492
## metaboliteSUCCINIC-2                       -0.671  0.336  0.558   0.491
## metabolitevaline                           -0.667  0.320  0.574   0.494
## genotypeKO:chowYellow (C8)                  0.178 -0.320 -0.302   0.040
## genotypeKO:metabolitearginine               0.355 -0.635 -0.147  -0.545
## genotypeKO:metaboliteCITRIC                 0.338 -0.605 -0.222  -0.229
## genotypeKO:metaboliteFUMARIC                0.354 -0.634 -0.205  -0.242
## genotypeKO:metaboliteglutamine              0.341 -0.611 -0.241  -0.268
## genotypeKO:metaboliteisoleucine             0.360 -0.645 -0.155  -0.271
## genotypeKO:metaboliteLACTIC                 0.354 -0.634 -0.207  -0.272
## genotypeKO:metaboliteLC even AC total       0.344 -0.616 -0.231  -0.254
## genotypeKO:metaboliteLC odd AC total        0.370 -0.663 -0.135  -0.289
## genotypeKO:metaboliteleucine                0.337 -0.604 -0.192  -0.238
## genotypeKO:metaboliteMALIC                  0.343 -0.615 -0.229  -0.259
## genotypeKO:metaboliteMCAC total             0.356 -0.638 -0.170  -0.263
## genotypeKO:metaboliteMETHYLSUCCINIC         0.348 -0.623 -0.179  -0.259
## genotypeKO:metaboliteSUCCINIC-2             0.346 -0.619 -0.215  -0.257
## genotypeKO:metabolitevaline                 0.337 -0.605 -0.243  -0.262
## chowYellow (C8):metabolitearginine          0.551 -0.187 -0.667  -0.786
## chowYellow (C8):metaboliteCITRIC            0.593 -0.228 -0.599  -0.443
## chowYellow (C8):metaboliteFUMARIC           0.575 -0.209 -0.607  -0.434
## chowYellow (C8):metaboliteglutamine         0.594 -0.239 -0.639  -0.404
## chowYellow (C8):metaboliteisoleucine        0.531 -0.155 -0.679  -0.375
## chowYellow (C8):metaboliteLACTIC            0.577 -0.213 -0.641  -0.404
## chowYellow (C8):metaboliteLC even AC total  0.570 -0.211 -0.609  -0.391
## chowYellow (C8):metaboliteLC odd AC total   0.543 -0.173 -0.631  -0.377
## chowYellow (C8):metaboliteleucine           0.561 -0.213 -0.638  -0.413
## chowYellow (C8):metaboliteMALIC             0.601 -0.245 -0.625  -0.427
## chowYellow (C8):metaboliteMCAC total        0.556 -0.197 -0.633  -0.399
## chowYellow (C8):metaboliteMETHYLSUCCINIC    0.545 -0.189 -0.633  -0.390
## chowYellow (C8):metaboliteSUCCINIC-2        0.570 -0.213 -0.615  -0.408
## chowYellow (C8):metabolitevaline            0.598 -0.239 -0.605  -0.415
##                                            mCITRI mFUMAR mtbltg mtblts
## genotypeKO                                                            
## chowYellow (C8)                                                       
## metabolitearginine                                                    
## metaboliteCITRIC                                                      
## metaboliteFUMARIC                           0.493                     
## metaboliteglutamine                         0.484  0.494              
## metaboliteisoleucine                        0.499  0.507  0.496       
## metaboliteLACTIC                            0.497  0.485  0.491  0.493
## metaboliteLC even AC total                  0.478  0.494  0.473  0.485
## metaboliteLC odd AC total                   0.482  0.512  0.507  0.511
## metaboliteleucine                           0.477  0.484  0.471  0.490
## metaboliteMALIC                             0.484  0.492  0.488  0.486
## metaboliteMCAC total                        0.492  0.493  0.480  0.497
## metaboliteMETHYLSUCCINIC                    0.483  0.493  0.485  0.492
## metaboliteSUCCINIC-2                        0.483  0.488  0.487  0.493
## metabolitevaline                            0.483  0.483  0.477  0.488
## genotypeKO:chowYellow (C8)                  0.010 -0.001  0.007  0.030
## genotypeKO:metabolitearginine              -0.223 -0.240 -0.262 -0.270
## genotypeKO:metaboliteCITRIC                -0.531 -0.258 -0.236 -0.269
## genotypeKO:metaboliteFUMARIC               -0.254 -0.541 -0.255 -0.284
## genotypeKO:metaboliteglutamine             -0.236 -0.259 -0.533 -0.263
## genotypeKO:metaboliteisoleucine            -0.263 -0.283 -0.259 -0.543
## genotypeKO:metaboliteLACTIC                -0.261 -0.243 -0.250 -0.259
## genotypeKO:metaboliteLC even AC total      -0.226 -0.259 -0.217 -0.243
## genotypeKO:metaboliteLC odd AC total       -0.233 -0.293 -0.279 -0.291
## genotypeKO:metaboliteleucine               -0.223 -0.241 -0.214 -0.253
## genotypeKO:metaboliteMALIC                 -0.235 -0.257 -0.245 -0.245
## genotypeKO:metaboliteMCAC total            -0.252 -0.258 -0.231 -0.266
## genotypeKO:metaboliteMETHYLSUCCINIC        -0.235 -0.258 -0.239 -0.258
## genotypeKO:metaboliteSUCCINIC-2            -0.234 -0.248 -0.243 -0.259
## genotypeKO:metabolitevaline                -0.233 -0.239 -0.223 -0.250
## chowYellow (C8):metabolitearginine         -0.433 -0.424 -0.405 -0.391
## chowYellow (C8):metaboliteCITRIC           -0.839 -0.411 -0.454 -0.408
## chowYellow (C8):metaboliteFUMARIC          -0.414 -0.821 -0.431 -0.397
## chowYellow (C8):metaboliteglutamine        -0.441 -0.423 -0.863 -0.412
## chowYellow (C8):metaboliteisoleucine       -0.379 -0.367 -0.408 -0.770
## chowYellow (C8):metaboliteLACTIC           -0.415 -0.440 -0.436 -0.420
## chowYellow (C8):metaboliteLC even AC total -0.439 -0.416 -0.458 -0.416
## chowYellow (C8):metaboliteLC odd AC total  -0.435 -0.365 -0.410 -0.388
## chowYellow (C8):metaboliteleucine          -0.433 -0.413 -0.442 -0.399
## chowYellow (C8):metaboliteMALIC            -0.451 -0.433 -0.451 -0.444
## chowYellow (C8):metaboliteMCAC total       -0.414 -0.401 -0.442 -0.410
## chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.414 -0.394 -0.419 -0.393
## chowYellow (C8):metaboliteSUCCINIC-2       -0.431 -0.421 -0.432 -0.408
## chowYellow (C8):metabolitevaline           -0.444 -0.446 -0.467 -0.433
##                                            mLACTI mLCeAt mLCoAt mtbltl
## genotypeKO                                                            
## chowYellow (C8)                                                       
## metabolitearginine                                                    
## metaboliteCITRIC                                                      
## metaboliteFUMARIC                                                     
## metaboliteglutamine                                                   
## metaboliteisoleucine                                                  
## metaboliteLACTIC                                                      
## metaboliteLC even AC total                  0.490                     
## metaboliteLC odd AC total                   0.501  0.504              
## metaboliteleucine                           0.489  0.483  0.499       
## metaboliteMALIC                             0.487  0.498  0.488  0.486
## metaboliteMCAC total                        0.493  0.486  0.500  0.486
## metaboliteMETHYLSUCCINIC                    0.490  0.487  0.493  0.492
## metaboliteSUCCINIC-2                        0.501  0.490  0.490  0.475
## metabolitevaline                            0.486  0.487  0.497  0.483
## genotypeKO:chowYellow (C8)                  0.010  0.000  0.016  0.022
## genotypeKO:metabolitearginine              -0.268 -0.252 -0.293 -0.239
## genotypeKO:metaboliteCITRIC                -0.264 -0.229 -0.242 -0.230
## genotypeKO:metaboliteFUMARIC               -0.241 -0.258 -0.299 -0.243
## genotypeKO:metaboliteglutamine             -0.252 -0.219 -0.289 -0.219
## genotypeKO:metaboliteisoleucine            -0.256 -0.241 -0.296 -0.254
## genotypeKO:metaboliteLACTIC                -0.537 -0.252 -0.278 -0.253
## genotypeKO:metaboliteLC even AC total      -0.251 -0.539 -0.285 -0.241
## genotypeKO:metaboliteLC odd AC total       -0.270 -0.278 -0.552 -0.270
## genotypeKO:metaboliteleucine               -0.248 -0.238 -0.274 -0.546
## genotypeKO:metaboliteMALIC                 -0.243 -0.265 -0.255 -0.246
## genotypeKO:metaboliteMCAC total            -0.256 -0.244 -0.276 -0.248
## genotypeKO:metaboliteMETHYLSUCCINIC        -0.250 -0.245 -0.264 -0.258
## genotypeKO:metaboliteSUCCINIC-2            -0.270 -0.251 -0.258 -0.228
## genotypeKO:metabolitevaline                -0.241 -0.244 -0.271 -0.241
## chowYellow (C8):metabolitearginine         -0.402 -0.401 -0.364 -0.421
## chowYellow (C8):metaboliteCITRIC           -0.425 -0.453 -0.428 -0.452
## chowYellow (C8):metaboliteFUMARIC          -0.448 -0.426 -0.361 -0.431
## chowYellow (C8):metaboliteglutamine        -0.433 -0.462 -0.400 -0.456
## chowYellow (C8):metaboliteisoleucine       -0.399 -0.404 -0.340 -0.393
## chowYellow (C8):metaboliteLACTIC           -0.834 -0.431 -0.396 -0.423
## chowYellow (C8):metaboliteLC even AC total -0.428 -0.811 -0.385 -0.423
## chowYellow (C8):metaboliteLC odd AC total  -0.410 -0.399 -0.733 -0.404
## chowYellow (C8):metaboliteleucine          -0.409 -0.417 -0.387 -0.763
## chowYellow (C8):metaboliteMALIC            -0.452 -0.426 -0.431 -0.438
## chowYellow (C8):metaboliteMCAC total       -0.414 -0.430 -0.391 -0.415
## chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.407 -0.410 -0.383 -0.388
## chowYellow (C8):metaboliteSUCCINIC-2       -0.402 -0.420 -0.406 -0.436
## chowYellow (C8):metabolitevaline           -0.446 -0.439 -0.405 -0.441
##                                            mMALIC mMCACt mMETHY mSUCCI
## genotypeKO                                                            
## chowYellow (C8)                                                       
## metabolitearginine                                                    
## metaboliteCITRIC                                                      
## metaboliteFUMARIC                                                     
## metaboliteglutamine                                                   
## metaboliteisoleucine                                                  
## metaboliteLACTIC                                                      
## metaboliteLC even AC total                                            
## metaboliteLC odd AC total                                             
## metaboliteleucine                                                     
## metaboliteMALIC                                                       
## metaboliteMCAC total                        0.496                     
## metaboliteMETHYLSUCCINIC                    0.483  0.495              
## metaboliteSUCCINIC-2                        0.485  0.496  0.495       
## metabolitevaline                            0.480  0.499  0.478  0.477
## genotypeKO:chowYellow (C8)                  0.015  0.003  0.019  0.007
## genotypeKO:metabolitearginine              -0.252 -0.268 -0.262 -0.256
## genotypeKO:metaboliteCITRIC                -0.235 -0.263 -0.244 -0.240
## genotypeKO:metaboliteFUMARIC               -0.252 -0.264 -0.263 -0.249
## genotypeKO:metaboliteglutamine             -0.244 -0.240 -0.248 -0.248
## genotypeKO:metaboliteisoleucine            -0.239 -0.272 -0.262 -0.259
## genotypeKO:metaboliteLACTIC                -0.240 -0.264 -0.257 -0.273
## genotypeKO:metaboliteLC even AC total      -0.261 -0.251 -0.251 -0.253
## genotypeKO:metaboliteLC odd AC total       -0.245 -0.277 -0.264 -0.254
## genotypeKO:metaboliteleucine               -0.239 -0.252 -0.261 -0.227
## genotypeKO:metaboliteMALIC                 -0.531 -0.269 -0.245 -0.244
## genotypeKO:metaboliteMCAC total            -0.258 -0.554 -0.267 -0.264
## genotypeKO:metaboliteMETHYLSUCCINIC        -0.235 -0.268 -0.552 -0.262
## genotypeKO:metaboliteSUCCINIC-2            -0.238 -0.269 -0.266 -0.543
## genotypeKO:metabolitevaline                -0.228 -0.276 -0.235 -0.227
## chowYellow (C8):metabolitearginine         -0.416 -0.391 -0.395 -0.406
## chowYellow (C8):metaboliteCITRIC           -0.451 -0.416 -0.428 -0.441
## chowYellow (C8):metaboliteFUMARIC          -0.432 -0.410 -0.410 -0.428
## chowYellow (C8):metaboliteglutamine        -0.445 -0.440 -0.436 -0.440
## chowYellow (C8):metaboliteisoleucine       -0.396 -0.383 -0.389 -0.398
## chowYellow (C8):metaboliteLACTIC           -0.445 -0.410 -0.420 -0.408
## chowYellow (C8):metaboliteLC even AC total -0.407 -0.413 -0.411 -0.417
## chowYellow (C8):metaboliteLC odd AC total  -0.420 -0.398 -0.383 -0.400
## chowYellow (C8):metaboliteleucine          -0.420 -0.400 -0.395 -0.432
## chowYellow (C8):metaboliteMALIC            -0.864 -0.415 -0.443 -0.448
## chowYellow (C8):metaboliteMCAC total       -0.422 -0.760 -0.392 -0.399
## chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.423 -0.378 -0.748 -0.391
## chowYellow (C8):metaboliteSUCCINIC-2       -0.436 -0.393 -0.399 -0.799
## chowYellow (C8):metabolitevaline           -0.463 -0.406 -0.443 -0.454
##                                            mtbltv gKO:Y( gntypKO:mtbltr
## genotypeKO                                                             
## chowYellow (C8)                                                        
## metabolitearginine                                                     
## metaboliteCITRIC                                                       
## metaboliteFUMARIC                                                      
## metaboliteglutamine                                                    
## metaboliteisoleucine                                                   
## metaboliteLACTIC                                                       
## metaboliteLC even AC total                                             
## metaboliteLC odd AC total                                              
## metaboliteleucine                                                      
## metaboliteMALIC                                                        
## metaboliteMCAC total                                                   
## metaboliteMETHYLSUCCINIC                                               
## metaboliteSUCCINIC-2                                                   
## metabolitevaline                                                       
## genotypeKO:chowYellow (C8)                  0.005                      
## genotypeKO:metabolitearginine              -0.254 -0.073               
## genotypeKO:metaboliteCITRIC                -0.232 -0.019  0.419        
## genotypeKO:metaboliteFUMARIC               -0.233  0.002  0.444        
## genotypeKO:metaboliteglutamine             -0.221 -0.013  0.492        
## genotypeKO:metaboliteisoleucine            -0.243 -0.055  0.498        
## genotypeKO:metaboliteLACTIC                -0.238 -0.019  0.500        
## genotypeKO:metaboliteLC even AC total      -0.240 -0.001  0.467        
## genotypeKO:metaboliteLC odd AC total       -0.259 -0.028  0.530        
## genotypeKO:metaboliteleucine               -0.234 -0.040  0.438        
## genotypeKO:metaboliteMALIC                 -0.228 -0.029  0.474        
## genotypeKO:metaboliteMCAC total            -0.263 -0.006  0.483        
## genotypeKO:metaboliteMETHYLSUCCINIC        -0.225 -0.035  0.475        
## genotypeKO:metaboliteSUCCINIC-2            -0.222 -0.012  0.471        
## genotypeKO:metabolitevaline                -0.529 -0.009  0.480        
## chowYellow (C8):metabolitearginine         -0.408  0.041  0.247        
## chowYellow (C8):metaboliteCITRIC           -0.449 -0.025  0.190        
## chowYellow (C8):metaboliteFUMARIC          -0.454 -0.058  0.184        
## chowYellow (C8):metaboliteglutamine        -0.457  0.027  0.126        
## chowYellow (C8):metaboliteisoleucine       -0.420  0.029  0.095        
## chowYellow (C8):metaboliteLACTIC           -0.443  0.002  0.128        
## chowYellow (C8):metaboliteLC even AC total -0.428 -0.031  0.113        
## chowYellow (C8):metaboliteLC odd AC total  -0.424 -0.050  0.096        
## chowYellow (C8):metaboliteleucine          -0.423  0.049  0.167        
## chowYellow (C8):metaboliteMALIC            -0.464  0.006  0.164        
## chowYellow (C8):metaboliteMCAC total       -0.402  0.005  0.137        
## chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.430  0.014  0.128        
## chowYellow (C8):metaboliteSUCCINIC-2       -0.450 -0.011  0.147        
## chowYellow (C8):metabolitevaline           -0.856 -0.003  0.141        
##                                            gKO:CI gKO:FU gntypKO:mtbltg
## genotypeKO                                                             
## chowYellow (C8)                                                        
## metabolitearginine                                                     
## metaboliteCITRIC                                                       
## metaboliteFUMARIC                                                      
## metaboliteglutamine                                                    
## metaboliteisoleucine                                                   
## metaboliteLACTIC                                                       
## metaboliteLC even AC total                                             
## metaboliteLC odd AC total                                              
## metaboliteleucine                                                      
## metaboliteMALIC                                                        
## metaboliteMCAC total                                                   
## metaboliteMETHYLSUCCINIC                                               
## metaboliteSUCCINIC-2                                                   
## metabolitevaline                                                       
## genotypeKO:chowYellow (C8)                                             
## genotypeKO:metabolitearginine                                          
## genotypeKO:metaboliteCITRIC                                            
## genotypeKO:metaboliteFUMARIC                0.477                      
## genotypeKO:metaboliteglutamine              0.443  0.478               
## genotypeKO:metaboliteisoleucine             0.495  0.523  0.485        
## genotypeKO:metaboliteLACTIC                 0.491  0.448  0.469        
## genotypeKO:metaboliteLC even AC total       0.425  0.479  0.407        
## genotypeKO:metaboliteLC odd AC total        0.438  0.541  0.524        
## genotypeKO:metaboliteleucine                0.421  0.446  0.402        
## genotypeKO:metaboliteMALIC                  0.442  0.474  0.459        
## genotypeKO:metaboliteMCAC total             0.475  0.476  0.433        
## genotypeKO:metaboliteMETHYLSUCCINIC         0.442  0.476  0.449        
## genotypeKO:metaboliteSUCCINIC-2             0.441  0.459  0.456        
## genotypeKO:metabolitevaline                 0.439  0.441  0.417        
## chowYellow (C8):metabolitearginine          0.196  0.181  0.143        
## chowYellow (C8):metaboliteCITRIC            0.291  0.132  0.210        
## chowYellow (C8):metaboliteFUMARIC           0.144  0.280  0.177        
## chowYellow (C8):metaboliteglutamine         0.191  0.159  0.349        
## chowYellow (C8):metaboliteisoleucine        0.097  0.078  0.152        
## chowYellow (C8):metaboliteLACTIC            0.147  0.194  0.186        
## chowYellow (C8):metaboliteLC even AC total  0.200  0.159  0.236        
## chowYellow (C8):metaboliteLC odd AC total   0.200  0.072  0.154        
## chowYellow (C8):metaboliteleucine           0.203  0.166  0.220        
## chowYellow (C8):metaboliteMALIC             0.207  0.174  0.207        
## chowYellow (C8):metaboliteMCAC total        0.161  0.139  0.215        
## chowYellow (C8):metaboliteMETHYLSUCCINIC    0.171  0.136  0.180        
## chowYellow (C8):metaboliteSUCCINIC-2        0.186  0.169  0.189        
## chowYellow (C8):metabolitevaline            0.191  0.197  0.235        
##                                            gntypKO:mtblts gKO:LA gKOeAt
## genotypeKO                                                             
## chowYellow (C8)                                                        
## metabolitearginine                                                     
## metaboliteCITRIC                                                       
## metaboliteFUMARIC                                                      
## metaboliteglutamine                                                    
## metaboliteisoleucine                                                   
## metaboliteLACTIC                                                       
## metaboliteLC even AC total                                             
## metaboliteLC odd AC total                                              
## metaboliteleucine                                                      
## metaboliteMALIC                                                        
## metaboliteMCAC total                                                   
## metaboliteMETHYLSUCCINIC                                               
## metaboliteSUCCINIC-2                                                   
## metabolitevaline                                                       
## genotypeKO:chowYellow (C8)                                             
## genotypeKO:metabolitearginine                                          
## genotypeKO:metaboliteCITRIC                                            
## genotypeKO:metaboliteFUMARIC                                           
## genotypeKO:metaboliteglutamine                                         
## genotypeKO:metaboliteisoleucine                                        
## genotypeKO:metaboliteLACTIC                 0.477                      
## genotypeKO:metaboliteLC even AC total       0.447          0.467       
## genotypeKO:metaboliteLC odd AC total        0.535          0.504  0.515
## genotypeKO:metaboliteleucine                0.466          0.463  0.442
## genotypeKO:metaboliteMALIC                  0.451          0.453  0.492
## genotypeKO:metaboliteMCAC total             0.491          0.477  0.452
## genotypeKO:metaboliteMETHYLSUCCINIC         0.475          0.465  0.455
## genotypeKO:metaboliteSUCCINIC-2             0.477          0.502  0.465
## genotypeKO:metabolitevaline                 0.460          0.450  0.453
## chowYellow (C8):metabolitearginine          0.119          0.138  0.138
## chowYellow (C8):metaboliteCITRIC            0.126          0.156  0.208
## chowYellow (C8):metaboliteFUMARIC           0.115          0.210  0.169
## chowYellow (C8):metaboliteglutamine         0.140          0.177  0.232
## chowYellow (C8):metaboliteisoleucine        0.225          0.137  0.146
## chowYellow (C8):metaboliteLACTIC            0.157          0.303  0.177
## chowYellow (C8):metaboliteLC even AC total  0.159          0.180  0.277
## chowYellow (C8):metaboliteLC odd AC total   0.114          0.154  0.134
## chowYellow (C8):metaboliteleucine           0.141          0.158  0.173
## chowYellow (C8):metaboliteMALIC             0.195          0.208  0.161
## chowYellow (C8):metaboliteMCAC total        0.156          0.161  0.192
## chowYellow (C8):metaboliteMETHYLSUCCINIC    0.133          0.158  0.164
## chowYellow (C8):metaboliteSUCCINIC-2        0.146          0.134  0.168
## chowYellow (C8):metabolitevaline            0.174          0.197  0.184
##                                            gKOoAt gntypKO:mtbltl gKO:MA
## genotypeKO                                                             
## chowYellow (C8)                                                        
## metabolitearginine                                                     
## metaboliteCITRIC                                                       
## metaboliteFUMARIC                                                      
## metaboliteglutamine                                                    
## metaboliteisoleucine                                                   
## metaboliteLACTIC                                                       
## metaboliteLC even AC total                                             
## metaboliteLC odd AC total                                              
## metaboliteleucine                                                      
## metaboliteMALIC                                                        
## metaboliteMCAC total                                                   
## metaboliteMETHYLSUCCINIC                                               
## metaboliteSUCCINIC-2                                                   
## metabolitevaline                                                       
## genotypeKO:chowYellow (C8)                                             
## genotypeKO:metabolitearginine                                          
## genotypeKO:metaboliteCITRIC                                            
## genotypeKO:metaboliteFUMARIC                                           
## genotypeKO:metaboliteglutamine                                         
## genotypeKO:metaboliteisoleucine                                        
## genotypeKO:metaboliteLACTIC                                            
## genotypeKO:metaboliteLC even AC total                                  
## genotypeKO:metaboliteLC odd AC total                                   
## genotypeKO:metaboliteleucine                0.495                      
## genotypeKO:metaboliteMALIC                  0.461  0.451               
## genotypeKO:metaboliteMCAC total             0.499  0.454          0.486
## genotypeKO:metaboliteMETHYLSUCCINIC         0.478  0.473          0.443
## genotypeKO:metaboliteSUCCINIC-2             0.468  0.417          0.449
## genotypeKO:metabolitevaline                 0.490  0.442          0.430
## chowYellow (C8):metabolitearginine          0.073  0.175          0.165
## chowYellow (C8):metaboliteCITRIC            0.164  0.208          0.204
## chowYellow (C8):metaboliteFUMARIC           0.052  0.179          0.179
## chowYellow (C8):metaboliteglutamine         0.121  0.221          0.200
## chowYellow (C8):metaboliteisoleucine        0.032  0.127          0.129
## chowYellow (C8):metaboliteLACTIC            0.117  0.163          0.203
## chowYellow (C8):metaboliteLC even AC total  0.105  0.172          0.140
## chowYellow (C8):metaboliteLC odd AC total   0.157  0.145          0.173
## chowYellow (C8):metaboliteleucine           0.121  0.219          0.178
## chowYellow (C8):metaboliteMALIC             0.171  0.183          0.341
## chowYellow (C8):metaboliteMCAC total        0.122  0.165          0.177
## chowYellow (C8):metaboliteMETHYLSUCCINIC    0.118  0.124          0.187
## chowYellow (C8):metaboliteSUCCINIC-2        0.144  0.199          0.196
## chowYellow (C8):metabolitevaline            0.124  0.189          0.227
##                                            gKO:Mt gKO:ME gKO:SU
## genotypeKO                                                     
## chowYellow (C8)                                                
## metabolitearginine                                             
## metaboliteCITRIC                                               
## metaboliteFUMARIC                                              
## metaboliteglutamine                                            
## metaboliteisoleucine                                           
## metaboliteLACTIC                                               
## metaboliteLC even AC total                                     
## metaboliteLC odd AC total                                      
## metaboliteleucine                                              
## metaboliteMALIC                                                
## metaboliteMCAC total                                           
## metaboliteMETHYLSUCCINIC                                       
## metaboliteSUCCINIC-2                                           
## metabolitevaline                                               
## genotypeKO:chowYellow (C8)                                     
## genotypeKO:metabolitearginine                                  
## genotypeKO:metaboliteCITRIC                                    
## genotypeKO:metaboliteFUMARIC                                   
## genotypeKO:metaboliteglutamine                                 
## genotypeKO:metaboliteisoleucine                                
## genotypeKO:metaboliteLACTIC                                    
## genotypeKO:metaboliteLC even AC total                          
## genotypeKO:metaboliteLC odd AC total                           
## genotypeKO:metaboliteleucine                                   
## genotypeKO:metaboliteMALIC                                     
## genotypeKO:metaboliteMCAC total                                
## genotypeKO:metaboliteMETHYLSUCCINIC         0.484              
## genotypeKO:metaboliteSUCCINIC-2             0.486  0.482       
## genotypeKO:metabolitevaline                 0.497  0.425  0.419
## chowYellow (C8):metabolitearginine          0.123  0.128  0.147
## chowYellow (C8):metaboliteCITRIC            0.143  0.165  0.188
## chowYellow (C8):metaboliteFUMARIC           0.142  0.141  0.173
## chowYellow (C8):metaboliteglutamine         0.194  0.186  0.191
## chowYellow (C8):metaboliteisoleucine        0.111  0.121  0.135
## chowYellow (C8):metaboliteLACTIC            0.141  0.159  0.136
## chowYellow (C8):metaboliteLC even AC total  0.155  0.152  0.161
## chowYellow (C8):metaboliteLC odd AC total   0.136  0.108  0.137
## chowYellow (C8):metaboliteleucine           0.144  0.135  0.202
## chowYellow (C8):metaboliteMALIC             0.143  0.194  0.202
## chowYellow (C8):metaboliteMCAC total        0.207  0.125  0.135
## chowYellow (C8):metaboliteMETHYLSUCCINIC    0.109  0.201  0.130
## chowYellow (C8):metaboliteSUCCINIC-2        0.122  0.131  0.260
## chowYellow (C8):metabolitevaline            0.126  0.193  0.213
##                                            gntypKO:mtbltv
## genotypeKO                                               
## chowYellow (C8)                                          
## metabolitearginine                                       
## metaboliteCITRIC                                         
## metaboliteFUMARIC                                        
## metaboliteglutamine                                      
## metaboliteisoleucine                                     
## metaboliteLACTIC                                         
## metaboliteLC even AC total                               
## metaboliteLC odd AC total                                
## metaboliteleucine                                        
## metaboliteMALIC                                          
## metaboliteMCAC total                                     
## metaboliteMETHYLSUCCINIC                                 
## metaboliteSUCCINIC-2                                     
## metabolitevaline                                         
## genotypeKO:chowYellow (C8)                               
## genotypeKO:metabolitearginine                            
## genotypeKO:metaboliteCITRIC                              
## genotypeKO:metaboliteFUMARIC                             
## genotypeKO:metaboliteglutamine                           
## genotypeKO:metaboliteisoleucine                          
## genotypeKO:metaboliteLACTIC                              
## genotypeKO:metaboliteLC even AC total                    
## genotypeKO:metaboliteLC odd AC total                     
## genotypeKO:metaboliteleucine                             
## genotypeKO:metaboliteMALIC                               
## genotypeKO:metaboliteMCAC total                          
## genotypeKO:metaboliteMETHYLSUCCINIC                      
## genotypeKO:metaboliteSUCCINIC-2                          
## genotypeKO:metabolitevaline                              
## chowYellow (C8):metabolitearginine          0.148        
## chowYellow (C8):metaboliteCITRIC            0.200        
## chowYellow (C8):metaboliteFUMARIC           0.219        
## chowYellow (C8):metaboliteglutamine         0.222        
## chowYellow (C8):metaboliteisoleucine        0.175        
## chowYellow (C8):metaboliteLACTIC            0.199        
## chowYellow (C8):metaboliteLC even AC total  0.180        
## chowYellow (C8):metaboliteLC odd AC total   0.179        
## chowYellow (C8):metaboliteleucine           0.184        
## chowYellow (C8):metaboliteMALIC             0.231        
## chowYellow (C8):metaboliteMCAC total        0.138        
## chowYellow (C8):metaboliteMETHYLSUCCINIC    0.200        
## chowYellow (C8):metaboliteSUCCINIC-2        0.222        
## chowYellow (C8):metabolitevaline            0.322        
##                                            chwYllw(C8):mtbltr cY(C8):C
## genotypeKO                                                            
## chowYellow (C8)                                                       
## metabolitearginine                                                    
## metaboliteCITRIC                                                      
## metaboliteFUMARIC                                                     
## metaboliteglutamine                                                   
## metaboliteisoleucine                                                  
## metaboliteLACTIC                                                      
## metaboliteLC even AC total                                            
## metaboliteLC odd AC total                                             
## metaboliteleucine                                                     
## metaboliteMALIC                                                       
## metaboliteMCAC total                                                  
## metaboliteMETHYLSUCCINIC                                              
## metaboliteSUCCINIC-2                                                  
## metabolitevaline                                                      
## genotypeKO:chowYellow (C8)                                            
## genotypeKO:metabolitearginine                                         
## genotypeKO:metaboliteCITRIC                                           
## genotypeKO:metaboliteFUMARIC                                          
## genotypeKO:metaboliteglutamine                                        
## genotypeKO:metaboliteisoleucine                                       
## genotypeKO:metaboliteLACTIC                                           
## genotypeKO:metaboliteLC even AC total                                 
## genotypeKO:metaboliteLC odd AC total                                  
## genotypeKO:metaboliteleucine                                          
## genotypeKO:metaboliteMALIC                                            
## genotypeKO:metaboliteMCAC total                                       
## genotypeKO:metaboliteMETHYLSUCCINIC                                   
## genotypeKO:metaboliteSUCCINIC-2                                       
## genotypeKO:metabolitevaline                                           
## chowYellow (C8):metabolitearginine                                    
## chowYellow (C8):metaboliteCITRIC            0.408                     
## chowYellow (C8):metaboliteFUMARIC           0.422              0.470  
## chowYellow (C8):metaboliteglutamine         0.512              0.454  
## chowYellow (C8):metaboliteisoleucine        0.502              0.500  
## chowYellow (C8):metaboliteLACTIC            0.492              0.486  
## chowYellow (C8):metaboliteLC even AC total  0.485              0.414  
## chowYellow (C8):metaboliteLC odd AC total   0.504              0.396  
## chowYellow (C8):metaboliteleucine           0.436              0.424  
## chowYellow (C8):metaboliteMALIC             0.471              0.440  
## chowYellow (C8):metaboliteMCAC total        0.459              0.447  
## chowYellow (C8):metaboliteMETHYLSUCCINIC    0.476              0.447  
## chowYellow (C8):metaboliteSUCCINIC-2        0.471              0.438  
## chowYellow (C8):metabolitevaline            0.481              0.435  
##                                            cY(C8):F chwYllw(C8):mtbltg
## genotypeKO                                                            
## chowYellow (C8)                                                       
## metabolitearginine                                                    
## metaboliteCITRIC                                                      
## metaboliteFUMARIC                                                     
## metaboliteglutamine                                                   
## metaboliteisoleucine                                                  
## metaboliteLACTIC                                                      
## metaboliteLC even AC total                                            
## metaboliteLC odd AC total                                             
## metaboliteleucine                                                     
## metaboliteMALIC                                                       
## metaboliteMCAC total                                                  
## metaboliteMETHYLSUCCINIC                                              
## metaboliteSUCCINIC-2                                                  
## metabolitevaline                                                      
## genotypeKO:chowYellow (C8)                                            
## genotypeKO:metabolitearginine                                         
## genotypeKO:metaboliteCITRIC                                           
## genotypeKO:metaboliteFUMARIC                                          
## genotypeKO:metaboliteglutamine                                        
## genotypeKO:metaboliteisoleucine                                       
## genotypeKO:metaboliteLACTIC                                           
## genotypeKO:metaboliteLC even AC total                                 
## genotypeKO:metaboliteLC odd AC total                                  
## genotypeKO:metaboliteleucine                                          
## genotypeKO:metaboliteMALIC                                            
## genotypeKO:metaboliteMCAC total                                       
## genotypeKO:metaboliteMETHYLSUCCINIC                                   
## genotypeKO:metaboliteSUCCINIC-2                                       
## genotypeKO:metabolitevaline                                           
## chowYellow (C8):metabolitearginine                                    
## chowYellow (C8):metaboliteCITRIC                                      
## chowYellow (C8):metaboliteFUMARIC                                     
## chowYellow (C8):metaboliteglutamine         0.487                     
## chowYellow (C8):metaboliteisoleucine        0.513    0.499            
## chowYellow (C8):metaboliteLACTIC            0.434    0.478            
## chowYellow (C8):metaboliteLC even AC total  0.460    0.407            
## chowYellow (C8):metaboliteLC odd AC total   0.528    0.519            
## chowYellow (C8):metaboliteleucine           0.447    0.405            
## chowYellow (C8):metaboliteMALIC             0.472    0.463            
## chowYellow (C8):metaboliteMCAC total        0.465    0.436            
## chowYellow (C8):metaboliteMETHYLSUCCINIC    0.477    0.442            
## chowYellow (C8):metaboliteSUCCINIC-2        0.459    0.452            
## chowYellow (C8):metabolitevaline            0.422    0.426            
##                                            chwYllw(C8):mtblts cY(C8):L
## genotypeKO                                                            
## chowYellow (C8)                                                       
## metabolitearginine                                                    
## metaboliteCITRIC                                                      
## metaboliteFUMARIC                                                     
## metaboliteglutamine                                                   
## metaboliteisoleucine                                                  
## metaboliteLACTIC                                                      
## metaboliteLC even AC total                                            
## metaboliteLC odd AC total                                             
## metaboliteleucine                                                     
## metaboliteMALIC                                                       
## metaboliteMCAC total                                                  
## metaboliteMETHYLSUCCINIC                                              
## metaboliteSUCCINIC-2                                                  
## metabolitevaline                                                      
## genotypeKO:chowYellow (C8)                                            
## genotypeKO:metabolitearginine                                         
## genotypeKO:metaboliteCITRIC                                           
## genotypeKO:metaboliteFUMARIC                                          
## genotypeKO:metaboliteglutamine                                        
## genotypeKO:metaboliteisoleucine                                       
## genotypeKO:metaboliteLACTIC                                           
## genotypeKO:metaboliteLC even AC total                                 
## genotypeKO:metaboliteLC odd AC total                                  
## genotypeKO:metaboliteleucine                                          
## genotypeKO:metaboliteMALIC                                            
## genotypeKO:metaboliteMCAC total                                       
## genotypeKO:metaboliteMETHYLSUCCINIC                                   
## genotypeKO:metaboliteSUCCINIC-2                                       
## genotypeKO:metabolitevaline                                           
## chowYellow (C8):metabolitearginine                                    
## chowYellow (C8):metaboliteCITRIC                                      
## chowYellow (C8):metaboliteFUMARIC                                     
## chowYellow (C8):metaboliteglutamine                                   
## chowYellow (C8):metaboliteisoleucine                                  
## chowYellow (C8):metaboliteLACTIC            0.484                     
## chowYellow (C8):metaboliteLC even AC total  0.475              0.455  
## chowYellow (C8):metaboliteLC odd AC total   0.488              0.478  
## chowYellow (C8):metaboliteleucine           0.478              0.468  
## chowYellow (C8):metaboliteMALIC             0.498              0.450  
## chowYellow (C8):metaboliteMCAC total        0.489              0.467  
## chowYellow (C8):metaboliteMETHYLSUCCINIC    0.472              0.461  
## chowYellow (C8):metaboliteSUCCINIC-2        0.472              0.501  
## chowYellow (C8):metabolitevaline            0.447              0.444  
##                                            cY(eAt cY(oAt
## genotypeKO                                              
## chowYellow (C8)                                         
## metabolitearginine                                      
## metaboliteCITRIC                                        
## metaboliteFUMARIC                                       
## metaboliteglutamine                                     
## metaboliteisoleucine                                    
## metaboliteLACTIC                                        
## metaboliteLC even AC total                              
## metaboliteLC odd AC total                               
## metaboliteleucine                                       
## metaboliteMALIC                                         
## metaboliteMCAC total                                    
## metaboliteMETHYLSUCCINIC                                
## metaboliteSUCCINIC-2                                    
## metabolitevaline                                        
## genotypeKO:chowYellow (C8)                              
## genotypeKO:metabolitearginine                           
## genotypeKO:metaboliteCITRIC                             
## genotypeKO:metaboliteFUMARIC                            
## genotypeKO:metaboliteglutamine                          
## genotypeKO:metaboliteisoleucine                         
## genotypeKO:metaboliteLACTIC                             
## genotypeKO:metaboliteLC even AC total                   
## genotypeKO:metaboliteLC odd AC total                    
## genotypeKO:metaboliteleucine                            
## genotypeKO:metaboliteMALIC                              
## genotypeKO:metaboliteMCAC total                         
## genotypeKO:metaboliteMETHYLSUCCINIC                     
## genotypeKO:metaboliteSUCCINIC-2                         
## genotypeKO:metabolitevaline                             
## chowYellow (C8):metabolitearginine                      
## chowYellow (C8):metaboliteCITRIC                        
## chowYellow (C8):metaboliteFUMARIC                       
## chowYellow (C8):metaboliteglutamine                     
## chowYellow (C8):metaboliteisoleucine                    
## chowYellow (C8):metaboliteLACTIC                        
## chowYellow (C8):metaboliteLC even AC total              
## chowYellow (C8):metaboliteLC odd AC total   0.487       
## chowYellow (C8):metaboliteleucine           0.444  0.483
## chowYellow (C8):metaboliteMALIC             0.505  0.464
## chowYellow (C8):metaboliteMCAC total        0.436  0.463
## chowYellow (C8):metaboliteMETHYLSUCCINIC    0.457  0.498
## chowYellow (C8):metaboliteSUCCINIC-2        0.463  0.483
## chowYellow (C8):metabolitevaline            0.455  0.448
##                                            chwYllw(C8):mtbltl cY(C8):MA
## genotypeKO                                                             
## chowYellow (C8)                                                        
## metabolitearginine                                                     
## metaboliteCITRIC                                                       
## metaboliteFUMARIC                                                      
## metaboliteglutamine                                                    
## metaboliteisoleucine                                                   
## metaboliteLACTIC                                                       
## metaboliteLC even AC total                                             
## metaboliteLC odd AC total                                              
## metaboliteleucine                                                      
## metaboliteMALIC                                                        
## metaboliteMCAC total                                                   
## metaboliteMETHYLSUCCINIC                                               
## metaboliteSUCCINIC-2                                                   
## metabolitevaline                                                       
## genotypeKO:chowYellow (C8)                                             
## genotypeKO:metabolitearginine                                          
## genotypeKO:metaboliteCITRIC                                            
## genotypeKO:metaboliteFUMARIC                                           
## genotypeKO:metaboliteglutamine                                         
## genotypeKO:metaboliteisoleucine                                        
## genotypeKO:metaboliteLACTIC                                            
## genotypeKO:metaboliteLC even AC total                                  
## genotypeKO:metaboliteLC odd AC total                                   
## genotypeKO:metaboliteleucine                                           
## genotypeKO:metaboliteMALIC                                             
## genotypeKO:metaboliteMCAC total                                        
## genotypeKO:metaboliteMETHYLSUCCINIC                                    
## genotypeKO:metaboliteSUCCINIC-2                                        
## genotypeKO:metabolitevaline                                            
## chowYellow (C8):metabolitearginine                                     
## chowYellow (C8):metaboliteCITRIC                                       
## chowYellow (C8):metaboliteFUMARIC                                      
## chowYellow (C8):metaboliteglutamine                                    
## chowYellow (C8):metaboliteisoleucine                                   
## chowYellow (C8):metaboliteLACTIC                                       
## chowYellow (C8):metaboliteLC even AC total                             
## chowYellow (C8):metaboliteLC odd AC total                              
## chowYellow (C8):metaboliteleucine                                      
## chowYellow (C8):metaboliteMALIC             0.452                      
## chowYellow (C8):metaboliteMCAC total        0.447              0.458   
## chowYellow (C8):metaboliteMETHYLSUCCINIC    0.464              0.437   
## chowYellow (C8):metaboliteSUCCINIC-2        0.408              0.444   
## chowYellow (C8):metabolitevaline            0.443              0.420   
##                                            cY(C8t cY(C8):ME cY(C8):S
## genotypeKO                                                          
## chowYellow (C8)                                                     
## metabolitearginine                                                  
## metaboliteCITRIC                                                    
## metaboliteFUMARIC                                                   
## metaboliteglutamine                                                 
## metaboliteisoleucine                                                
## metaboliteLACTIC                                                    
## metaboliteLC even AC total                                          
## metaboliteLC odd AC total                                           
## metaboliteleucine                                                   
## metaboliteMALIC                                                     
## metaboliteMCAC total                                                
## metaboliteMETHYLSUCCINIC                                            
## metaboliteSUCCINIC-2                                                
## metabolitevaline                                                    
## genotypeKO:chowYellow (C8)                                          
## genotypeKO:metabolitearginine                                       
## genotypeKO:metaboliteCITRIC                                         
## genotypeKO:metaboliteFUMARIC                                        
## genotypeKO:metaboliteglutamine                                      
## genotypeKO:metaboliteisoleucine                                     
## genotypeKO:metaboliteLACTIC                                         
## genotypeKO:metaboliteLC even AC total                               
## genotypeKO:metaboliteLC odd AC total                                
## genotypeKO:metaboliteleucine                                        
## genotypeKO:metaboliteMALIC                                          
## genotypeKO:metaboliteMCAC total                                     
## genotypeKO:metaboliteMETHYLSUCCINIC                                 
## genotypeKO:metaboliteSUCCINIC-2                                     
## genotypeKO:metabolitevaline                                         
## chowYellow (C8):metabolitearginine                                  
## chowYellow (C8):metaboliteCITRIC                                    
## chowYellow (C8):metaboliteFUMARIC                                   
## chowYellow (C8):metaboliteglutamine                                 
## chowYellow (C8):metaboliteisoleucine                                
## chowYellow (C8):metaboliteLACTIC                                    
## chowYellow (C8):metaboliteLC even AC total                          
## chowYellow (C8):metaboliteLC odd AC total                           
## chowYellow (C8):metaboliteleucine                                   
## chowYellow (C8):metaboliteMALIC                                     
## chowYellow (C8):metaboliteMCAC total                                
## chowYellow (C8):metaboliteMETHYLSUCCINIC    0.478                   
## chowYellow (C8):metaboliteSUCCINIC-2        0.483  0.478            
## chowYellow (C8):metabolitevaline            0.494  0.421     0.416  
## 
## Standardized Within-Group Residuals:
##          Min           Q1          Med           Q3          Max 
## -14.11233569  -0.11922815   0.02994209   0.21328345   4.72986981 
## 
## Number of Observations: 645
## Number of Groups: 43
```
