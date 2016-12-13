---
title: "Metabolomics of very long-chain aclCoA dehydrogenase knockout mice"
date: "2016-12-13 15:19:51"
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
  aes(x = logValue, y = metabolite, color = activity, fill = activity) +
  geom_jitter(alpha = 1/2) +
  facet_wrap(~ genotype) +
  scale_x_continuous("log10 scale") +
  scale_y_discrete("Metabolite") +
  scale_color_brewer("Activity", palette = "Set1") +
  scale_fill_brewer("Activity", palette = "Set1") +
  theme_bw()
```

![plot of chunk plotDataAim1](../figures/plotDataAim1-1.png)

```r
ggsave("../figures/plotDataAim1.png", width = 10, height = 8, units = "in")
```


```r
kable(summarizeOutcome(D2))
```



|              |       Min.| X1st.Qu.| Median|    Mean| X3rd.Qu.|    Max.|
|:-------------|----------:|--------:|------:|-------:|--------:|-------:|
|nominal       |  0.0000001|  0.84630| 3.2620| 24.7700|   12.070| 408.800|
|log-transform | -6.9550000| -0.07246| 0.5134|  0.4957|    1.082|   2.612|

```r
ggplot(D2) +
  aes(x = logValue, y = metabolite, color = chow, fill = chow) +
  geom_jitter(alpha = 1/2) +
  facet_wrap(~ genotype) +
  scale_x_continuous("log10 scale") +
  scale_y_discrete("Metabolite") +
  scale_color_brewer("Chow", palette = "Set1") +
  scale_fill_brewer("Chow", palette = "Set1") +
  theme_bw()
```

![plot of chunk plotDataAim2](../figures/plotDataAim2-1.png)

```r
ggsave("../figures/plotDataAim2.png", width = 10, height = 8, units = "in")
```


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
                   activity * metabolite +
                   activity * metabolite * genotype)
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
##                              numDF denDF   F-value p-value
## (Intercept)                      1   524  149301.5  <.0001
## genotype                         1    38   98911.6  <.0001
## activity                         1    38   63134.6  <.0001
## metabolite                      14   524 2824393.4  <.0001
## genotype:activity                1    38       0.0  0.9801
## genotype:metabolite             14   524       3.4  <.0001
## activity:metabolite             14   524       0.8  0.7007
## genotype:activity:metabolite    14   524       1.8  0.0375
```

```r
summary(M)
```

```
## Linear mixed-effects model fit by REML
##  Data: D1 
##       AIC      BIC    logLik
##   1713.92 2437.281 -689.9601
## 
## Random effects:
##  Formula: ~1 | id
##         (Intercept)  Residual
## StdDev:  0.06348574 0.7260516
## 
## Correlation Structure: General
##  Formula: ~1 | id 
##  Parameter estimate(s):
##  Correlation: 
##    1      2      3      4      5      6      7      8      9      10    
## 2  -0.020                                                               
## 3   0.080 -0.014                                                        
## 4   0.126  0.007 -0.046                                                 
## 5   0.081  0.049 -0.092  0.127                                          
## 6   0.016  0.039  0.137 -0.105  0.003                                   
## 7   0.054 -0.008  0.009  0.004 -0.168  0.069                            
## 8  -0.022 -0.057 -0.091 -0.025  0.030 -0.027 -0.223                     
## 9   0.016  0.041 -0.009  0.014  0.032 -0.047 -0.201  0.099              
## 10  0.028  0.114  0.009 -0.012 -0.007 -0.004  0.072 -0.133 -0.008       
## 11  0.026 -0.095  0.072  0.034  0.105  0.037  0.025  0.235 -0.240  0.035
## 12 -0.042  0.229 -0.014 -0.008  0.147  0.121  0.060 -0.151 -0.016  0.021
## 13  0.003  0.071 -0.107  0.013  0.066  0.056  0.065 -0.043 -0.030  0.139
## 14  0.066 -0.049  0.106  0.019 -0.108  0.092  0.088 -0.187 -0.005 -0.060
## 15 -0.017  0.147  0.087  0.012  0.082  0.066 -0.042  0.018  0.107  0.024
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
## 12 -0.255                     
## 13  0.032  0.043              
## 14 -0.044  0.218 -0.085       
## 15 -0.032  0.096  0.017 -0.010
## Fixed effects: list(fixed) 
##                                                           Value Std.Error
## (Intercept)                                           0.7800250 0.2429359
## genotypeKO                                            0.1111659 0.3348671
## activityExercise                                      0.2776343 0.3275776
## metabolitearginine                                   -0.3534360 0.3336614
## metaboliteCITRIC                                      0.2699043 0.3422640
## metaboliteFUMARIC                                    -1.7086321 0.3422640
## metaboliteglutamine                                   1.2558975 0.3336614
## metaboliteisoleucine                                 -0.5368730 0.3336614
## metaboliteLACTIC                                      1.6755101 0.3422640
## metaboliteLCAC total                                  0.1169499 0.3336614
## metaboliteleucine                                    -0.3706435 0.3336614
## metaboliteMALIC                                       0.1387884 0.3422640
## metaboliteMCAC Total                                 -1.2621611 0.3336614
## metaboliteMETHYLSUCCINIC                             -2.0968883 0.3422640
## metabolitePYRUVIC_P2P                                -0.9377425 0.3422640
## metaboliteSUCCINIC-2                                  0.3515290 0.3422640
## metabolitevaline                                     -0.3999569 0.3336614
## genotypeKO:activityExercise                          -0.2725876 0.4271739
## genotypeKO:metabolitearginine                        -0.0845930 0.4655750
## genotypeKO:metaboliteCITRIC                          -0.1808884 0.4717783
## genotypeKO:metaboliteFUMARIC                         -0.1459773 0.4717783
## genotypeKO:metaboliteglutamine                       -0.0934134 0.4655750
## genotypeKO:metaboliteisoleucine                      -0.1551589 0.4655750
## genotypeKO:metaboliteLACTIC                          -0.2130032 0.4717783
## genotypeKO:metaboliteLCAC total                       0.2585912 0.4655750
## genotypeKO:metaboliteleucine                         -0.1519655 0.4655750
## genotypeKO:metaboliteMALIC                           -0.1775990 0.4717783
## genotypeKO:metaboliteMCAC Total                      -0.2590161 0.4655750
## genotypeKO:metaboliteMETHYLSUCCINIC                  -0.3177490 0.4717783
## genotypeKO:metabolitePYRUVIC_P2P                     -2.5513187 0.4717783
## genotypeKO:metaboliteSUCCINIC-2                      -0.1502577 0.4717783
## genotypeKO:metabolitevaline                          -0.1344774 0.4655750
## activityExercise:metabolitearginine                  -0.5175138 0.4551654
## activityExercise:metaboliteCITRIC                    -0.2035335 0.4615087
## activityExercise:metaboliteFUMARIC                   -0.2704327 0.4615087
## activityExercise:metaboliteglutamine                 -0.3376519 0.4551654
## activityExercise:metaboliteisoleucine                -0.2834567 0.4551654
## activityExercise:metaboliteLACTIC                    -0.2581078 0.4615087
## activityExercise:metaboliteLCAC total                -0.4170637 0.4551654
## activityExercise:metaboliteleucine                   -0.3160980 0.4551654
## activityExercise:metaboliteMALIC                     -0.1439485 0.4615087
## activityExercise:metaboliteMCAC Total                -0.3739947 0.4551654
## activityExercise:metaboliteMETHYLSUCCINIC            -0.3997036 0.4615087
## activityExercise:metabolitePYRUVIC_P2P               -1.3300154 0.4615087
## activityExercise:metaboliteSUCCINIC-2                -0.3087569 0.4615087
## activityExercise:metabolitevaline                    -0.2501754 0.4551654
## genotypeKO:activityExercise:metabolitearginine        0.2727438 0.5742045
## genotypeKO:activityExercise:metaboliteCITRIC          0.0996485 0.6048647
## genotypeKO:activityExercise:metaboliteFUMARIC         0.2290113 0.6014938
## genotypeKO:activityExercise:metaboliteglutamine       0.0838461 0.5947977
## genotypeKO:activityExercise:metaboliteisoleucine      0.2896583 0.5969625
## genotypeKO:activityExercise:metaboliteLACTIC          0.2524058 0.5989676
## genotypeKO:activityExercise:metaboliteLCAC total      0.3834262 0.6047703
## genotypeKO:activityExercise:metaboliteleucine        -0.4830162 0.6010027
## genotypeKO:activityExercise:metaboliteMALIC           0.1327393 0.6047662
## genotypeKO:activityExercise:metaboliteMCAC Total      0.2062359 0.6050828
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC  0.5375277 0.6001742
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P     2.1895025 0.6137569
## genotypeKO:activityExercise:metaboliteSUCCINIC-2      0.0678668 0.6165921
## genotypeKO:activityExercise:metabolitevaline          0.5488350 0.5997617
##                                                       DF   t-value p-value
## (Intercept)                                          524  3.210826  0.0014
## genotypeKO                                            38  0.331970  0.7417
## activityExercise                                      38  0.847538  0.4020
## metabolitearginine                                   524 -1.059266  0.2900
## metaboliteCITRIC                                     524  0.788585  0.4307
## metaboliteFUMARIC                                    524 -4.992147  0.0000
## metaboliteglutamine                                  524  3.763988  0.0002
## metaboliteisoleucine                                 524 -1.609036  0.1082
## metaboliteLACTIC                                     524  4.895374  0.0000
## metaboliteLCAC total                                 524  0.350505  0.7261
## metaboliteleucine                                    524 -1.110837  0.2671
## metaboliteMALIC                                      524  0.405501  0.6853
## metaboliteMCAC Total                                 524 -3.782761  0.0002
## metaboliteMETHYLSUCCINIC                             524 -6.126523  0.0000
## metabolitePYRUVIC_P2P                                524 -2.739822  0.0064
## metaboliteSUCCINIC-2                                 524  1.027070  0.3049
## metabolitevaline                                     524 -1.198691  0.2312
## genotypeKO:activityExercise                           38 -0.638119  0.5272
## genotypeKO:metabolitearginine                        524 -0.181696  0.8559
## genotypeKO:metaboliteCITRIC                          524 -0.383418  0.7016
## genotypeKO:metaboliteFUMARIC                         524 -0.309419  0.7571
## genotypeKO:metaboliteglutamine                       524 -0.200641  0.8411
## genotypeKO:metaboliteisoleucine                      524 -0.333263  0.7391
## genotypeKO:metaboliteLACTIC                          524 -0.451490  0.6518
## genotypeKO:metaboliteLCAC total                      524  0.555423  0.5788
## genotypeKO:metaboliteleucine                         524 -0.326404  0.7442
## genotypeKO:metaboliteMALIC                           524 -0.376446  0.7067
## genotypeKO:metaboliteMCAC Total                      524 -0.556336  0.5782
## genotypeKO:metaboliteMETHYLSUCCINIC                  524 -0.673513  0.5009
## genotypeKO:metabolitePYRUVIC_P2P                     524 -5.407876  0.0000
## genotypeKO:metaboliteSUCCINIC-2                      524 -0.318492  0.7502
## genotypeKO:metabolitevaline                          524 -0.288841  0.7728
## activityExercise:metabolitearginine                  524 -1.136980  0.2561
## activityExercise:metaboliteCITRIC                    524 -0.441018  0.6594
## activityExercise:metaboliteFUMARIC                   524 -0.585975  0.5581
## activityExercise:metaboliteglutamine                 524 -0.741822  0.4585
## activityExercise:metaboliteisoleucine                524 -0.622755  0.5337
## activityExercise:metaboliteLACTIC                    524 -0.559270  0.5762
## activityExercise:metaboliteLCAC total                524 -0.916290  0.3599
## activityExercise:metaboliteleucine                   524 -0.694469  0.4877
## activityExercise:metaboliteMALIC                     524 -0.311909  0.7552
## activityExercise:metaboliteMCAC Total                524 -0.821668  0.4116
## activityExercise:metaboliteMETHYLSUCCINIC            524 -0.866080  0.3868
## activityExercise:metabolitePYRUVIC_P2P               524 -2.881886  0.0041
## activityExercise:metaboliteSUCCINIC-2                524 -0.669016  0.5038
## activityExercise:metabolitevaline                    524 -0.549636  0.5828
## genotypeKO:activityExercise:metabolitearginine       524  0.474994  0.6350
## genotypeKO:activityExercise:metaboliteCITRIC         524  0.164745  0.8692
## genotypeKO:activityExercise:metaboliteFUMARIC        524  0.380738  0.7036
## genotypeKO:activityExercise:metaboliteglutamine      524  0.140966  0.8880
## genotypeKO:activityExercise:metaboliteisoleucine     524  0.485220  0.6277
## genotypeKO:activityExercise:metaboliteLACTIC         524  0.421401  0.6736
## genotypeKO:activityExercise:metaboliteLCAC total     524  0.634003  0.5264
## genotypeKO:activityExercise:metaboliteleucine        524 -0.803684  0.4219
## genotypeKO:activityExercise:metaboliteMALIC          524  0.219489  0.8264
## genotypeKO:activityExercise:metaboliteMCAC Total     524  0.340839  0.7334
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC 524  0.895620  0.3709
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    524  3.567377  0.0004
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     524  0.110068  0.9124
## genotypeKO:activityExercise:metabolitevaline         524  0.915088  0.3606
##  Correlation: 
##                                                      (Intr) gntyKO actvtE
## genotypeKO                                           -0.725              
## activityExercise                                     -0.742  0.538       
## metabolitearginine                                   -0.723  0.525  0.536
## metaboliteCITRIC                                     -0.704  0.511  0.522
## metaboliteFUMARIC                                    -0.704  0.511  0.522
## metaboliteglutamine                                  -0.723  0.525  0.536
## metaboliteisoleucine                                 -0.723  0.525  0.536
## metaboliteLACTIC                                     -0.704  0.511  0.522
## metaboliteLCAC total                                 -0.723  0.525  0.536
## metaboliteleucine                                    -0.723  0.525  0.536
## metaboliteMALIC                                      -0.704  0.511  0.522
## metaboliteMCAC Total                                 -0.723  0.525  0.536
## metaboliteMETHYLSUCCINIC                             -0.704  0.511  0.522
## metabolitePYRUVIC_P2P                                -0.704  0.511  0.522
## metaboliteSUCCINIC-2                                 -0.704  0.511  0.522
## metabolitevaline                                     -0.723  0.525  0.536
## genotypeKO:activityExercise                           0.569 -0.784 -0.767
## genotypeKO:metabolitearginine                         0.518 -0.714 -0.384
## genotypeKO:metaboliteCITRIC                           0.511 -0.704 -0.379
## genotypeKO:metaboliteFUMARIC                          0.511 -0.704 -0.379
## genotypeKO:metaboliteglutamine                        0.518 -0.714 -0.384
## genotypeKO:metaboliteisoleucine                       0.518 -0.714 -0.384
## genotypeKO:metaboliteLACTIC                           0.511 -0.704 -0.379
## genotypeKO:metaboliteLCAC total                       0.518 -0.714 -0.384
## genotypeKO:metaboliteleucine                          0.518 -0.714 -0.384
## genotypeKO:metaboliteMALIC                            0.511 -0.704 -0.379
## genotypeKO:metaboliteMCAC Total                       0.518 -0.714 -0.384
## genotypeKO:metaboliteMETHYLSUCCINIC                   0.511 -0.704 -0.379
## genotypeKO:metabolitePYRUVIC_P2P                      0.511 -0.704 -0.379
## genotypeKO:metaboliteSUCCINIC-2                       0.511 -0.704 -0.379
## genotypeKO:metabolitevaline                           0.518 -0.714 -0.384
## activityExercise:metabolitearginine                   0.530 -0.385 -0.715
## activityExercise:metaboliteCITRIC                     0.522 -0.379 -0.704
## activityExercise:metaboliteFUMARIC                    0.522 -0.379 -0.704
## activityExercise:metaboliteglutamine                  0.530 -0.385 -0.715
## activityExercise:metaboliteisoleucine                 0.530 -0.385 -0.715
## activityExercise:metaboliteLACTIC                     0.522 -0.379 -0.704
## activityExercise:metaboliteLCAC total                 0.530 -0.385 -0.715
## activityExercise:metaboliteleucine                    0.530 -0.385 -0.715
## activityExercise:metaboliteMALIC                      0.522 -0.379 -0.704
## activityExercise:metaboliteMCAC Total                 0.530 -0.385 -0.715
## activityExercise:metaboliteMETHYLSUCCINIC             0.522 -0.379 -0.704
## activityExercise:metabolitePYRUVIC_P2P                0.522 -0.379 -0.704
## activityExercise:metaboliteSUCCINIC-2                 0.522 -0.379 -0.704
## activityExercise:metabolitevaline                     0.530 -0.385 -0.715
## genotypeKO:activityExercise:metabolitearginine       -0.420  0.579  0.566
## genotypeKO:activityExercise:metaboliteCITRIC         -0.399  0.549  0.537
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.401  0.553  0.540
## genotypeKO:activityExercise:metaboliteglutamine      -0.406  0.559  0.547
## genotypeKO:activityExercise:metaboliteisoleucine     -0.404  0.557  0.545
## genotypeKO:activityExercise:metaboliteLACTIC         -0.403  0.555  0.543
## genotypeKO:activityExercise:metaboliteLCAC total     -0.399  0.550  0.538
## genotypeKO:activityExercise:metaboliteleucine        -0.401  0.553  0.541
## genotypeKO:activityExercise:metaboliteMALIC          -0.399  0.550  0.538
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.399  0.549  0.537
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.402  0.554  0.542
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.393  0.541  0.530
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.391  0.539  0.527
## genotypeKO:activityExercise:metabolitevaline         -0.402  0.554  0.542
##                                                      mtbltr mCITRI mFUMAR
## genotypeKO                                                               
## activityExercise                                                         
## metabolitearginine                                                       
## metaboliteCITRIC                                      0.513              
## metaboliteFUMARIC                                     0.513  0.500       
## metaboliteglutamine                                   0.526  0.513  0.513
## metaboliteisoleucine                                  0.526  0.513  0.513
## metaboliteLACTIC                                      0.513  0.500  0.500
## metaboliteLCAC total                                  0.526  0.513  0.513
## metaboliteleucine                                     0.526  0.513  0.513
## metaboliteMALIC                                       0.513  0.500  0.500
## metaboliteMCAC Total                                  0.526  0.513  0.513
## metaboliteMETHYLSUCCINIC                              0.513  0.500  0.500
## metabolitePYRUVIC_P2P                                 0.513  0.500  0.500
## metaboliteSUCCINIC-2                                  0.513  0.500  0.500
## metabolitevaline                                      0.526  0.513  0.513
## genotypeKO:activityExercise                          -0.411 -0.401 -0.401
## genotypeKO:metabolitearginine                        -0.717 -0.368 -0.368
## genotypeKO:metaboliteCITRIC                          -0.372 -0.725 -0.363
## genotypeKO:metaboliteFUMARIC                         -0.372 -0.363 -0.725
## genotypeKO:metaboliteglutamine                       -0.377 -0.368 -0.368
## genotypeKO:metaboliteisoleucine                      -0.377 -0.368 -0.368
## genotypeKO:metaboliteLACTIC                          -0.372 -0.363 -0.363
## genotypeKO:metaboliteLCAC total                      -0.377 -0.368 -0.368
## genotypeKO:metaboliteleucine                         -0.377 -0.368 -0.368
## genotypeKO:metaboliteMALIC                           -0.372 -0.363 -0.363
## genotypeKO:metaboliteMCAC Total                      -0.377 -0.368 -0.368
## genotypeKO:metaboliteMETHYLSUCCINIC                  -0.372 -0.363 -0.363
## genotypeKO:metabolitePYRUVIC_P2P                     -0.372 -0.363 -0.363
## genotypeKO:metaboliteSUCCINIC-2                      -0.372 -0.363 -0.363
## genotypeKO:metabolitevaline                          -0.377 -0.368 -0.368
## activityExercise:metabolitearginine                  -0.733 -0.376 -0.376
## activityExercise:metaboliteCITRIC                    -0.380 -0.742 -0.371
## activityExercise:metaboliteFUMARIC                   -0.380 -0.371 -0.742
## activityExercise:metaboliteglutamine                 -0.386 -0.376 -0.376
## activityExercise:metaboliteisoleucine                -0.386 -0.376 -0.376
## activityExercise:metaboliteLACTIC                    -0.380 -0.371 -0.371
## activityExercise:metaboliteLCAC total                -0.386 -0.376 -0.376
## activityExercise:metaboliteleucine                   -0.386 -0.376 -0.376
## activityExercise:metaboliteMALIC                     -0.380 -0.371 -0.371
## activityExercise:metaboliteMCAC Total                -0.386 -0.376 -0.376
## activityExercise:metaboliteMETHYLSUCCINIC            -0.380 -0.371 -0.371
## activityExercise:metabolitePYRUVIC_P2P               -0.380 -0.371 -0.371
## activityExercise:metaboliteSUCCINIC-2                -0.380 -0.371 -0.371
## activityExercise:metabolitevaline                    -0.386 -0.376 -0.376
## genotypeKO:activityExercise:metabolitearginine        0.581  0.298  0.298
## genotypeKO:activityExercise:metaboliteCITRIC          0.290  0.566  0.283
## genotypeKO:activityExercise:metaboliteFUMARIC         0.292  0.285  0.569
## genotypeKO:activityExercise:metaboliteglutamine       0.295  0.288  0.288
## genotypeKO:activityExercise:metaboliteisoleucine      0.294  0.287  0.287
## genotypeKO:activityExercise:metaboliteLACTIC          0.293  0.286  0.286
## genotypeKO:activityExercise:metaboliteLCAC total      0.290  0.283  0.283
## genotypeKO:activityExercise:metaboliteleucine         0.292  0.285  0.285
## genotypeKO:activityExercise:metaboliteMALIC           0.290  0.283  0.283
## genotypeKO:activityExercise:metaboliteMCAC Total      0.290  0.283  0.283
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC  0.292  0.285  0.285
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P     0.286  0.279  0.279
## genotypeKO:activityExercise:metaboliteSUCCINIC-2      0.285  0.278  0.278
## genotypeKO:activityExercise:metabolitevaline          0.293  0.285  0.285
##                                                      mtbltg mtblts mLACTI
## genotypeKO                                                               
## activityExercise                                                         
## metabolitearginine                                                       
## metaboliteCITRIC                                                         
## metaboliteFUMARIC                                                        
## metaboliteglutamine                                                      
## metaboliteisoleucine                                  0.526              
## metaboliteLACTIC                                      0.513  0.513       
## metaboliteLCAC total                                  0.526  0.526  0.513
## metaboliteleucine                                     0.526  0.526  0.513
## metaboliteMALIC                                       0.513  0.513  0.500
## metaboliteMCAC Total                                  0.526  0.526  0.513
## metaboliteMETHYLSUCCINIC                              0.513  0.513  0.500
## metabolitePYRUVIC_P2P                                 0.513  0.513  0.500
## metaboliteSUCCINIC-2                                  0.513  0.513  0.500
## metabolitevaline                                      0.526  0.526  0.513
## genotypeKO:activityExercise                          -0.411 -0.411 -0.401
## genotypeKO:metabolitearginine                        -0.377 -0.377 -0.368
## genotypeKO:metaboliteCITRIC                          -0.372 -0.372 -0.363
## genotypeKO:metaboliteFUMARIC                         -0.372 -0.372 -0.363
## genotypeKO:metaboliteglutamine                       -0.717 -0.377 -0.368
## genotypeKO:metaboliteisoleucine                      -0.377 -0.717 -0.368
## genotypeKO:metaboliteLACTIC                          -0.372 -0.372 -0.725
## genotypeKO:metaboliteLCAC total                      -0.377 -0.377 -0.368
## genotypeKO:metaboliteleucine                         -0.377 -0.377 -0.368
## genotypeKO:metaboliteMALIC                           -0.372 -0.372 -0.363
## genotypeKO:metaboliteMCAC Total                      -0.377 -0.377 -0.368
## genotypeKO:metaboliteMETHYLSUCCINIC                  -0.372 -0.372 -0.363
## genotypeKO:metabolitePYRUVIC_P2P                     -0.372 -0.372 -0.363
## genotypeKO:metaboliteSUCCINIC-2                      -0.372 -0.372 -0.363
## genotypeKO:metabolitevaline                          -0.377 -0.377 -0.368
## activityExercise:metabolitearginine                  -0.386 -0.386 -0.376
## activityExercise:metaboliteCITRIC                    -0.380 -0.380 -0.371
## activityExercise:metaboliteFUMARIC                   -0.380 -0.380 -0.371
## activityExercise:metaboliteglutamine                 -0.733 -0.386 -0.376
## activityExercise:metaboliteisoleucine                -0.386 -0.733 -0.376
## activityExercise:metaboliteLACTIC                    -0.380 -0.380 -0.742
## activityExercise:metaboliteLCAC total                -0.386 -0.386 -0.376
## activityExercise:metaboliteleucine                   -0.386 -0.386 -0.376
## activityExercise:metaboliteMALIC                     -0.380 -0.380 -0.371
## activityExercise:metaboliteMCAC Total                -0.386 -0.386 -0.376
## activityExercise:metaboliteMETHYLSUCCINIC            -0.380 -0.380 -0.371
## activityExercise:metabolitePYRUVIC_P2P               -0.380 -0.380 -0.371
## activityExercise:metaboliteSUCCINIC-2                -0.380 -0.380 -0.371
## activityExercise:metabolitevaline                    -0.386 -0.386 -0.376
## genotypeKO:activityExercise:metabolitearginine        0.306  0.306  0.298
## genotypeKO:activityExercise:metaboliteCITRIC          0.290  0.290  0.283
## genotypeKO:activityExercise:metaboliteFUMARIC         0.292  0.292  0.285
## genotypeKO:activityExercise:metaboliteglutamine       0.561  0.295  0.288
## genotypeKO:activityExercise:metaboliteisoleucine      0.294  0.559  0.287
## genotypeKO:activityExercise:metaboliteLACTIC          0.293  0.293  0.571
## genotypeKO:activityExercise:metaboliteLCAC total      0.290  0.290  0.283
## genotypeKO:activityExercise:metaboliteleucine         0.292  0.292  0.285
## genotypeKO:activityExercise:metaboliteMALIC           0.290  0.290  0.283
## genotypeKO:activityExercise:metaboliteMCAC Total      0.290  0.290  0.283
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC  0.292  0.292  0.285
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P     0.286  0.286  0.279
## genotypeKO:activityExercise:metaboliteSUCCINIC-2      0.285  0.285  0.278
## genotypeKO:activityExercise:metabolitevaline          0.293  0.293  0.285
##                                                      mLCACt mtbltl mMALIC
## genotypeKO                                                               
## activityExercise                                                         
## metabolitearginine                                                       
## metaboliteCITRIC                                                         
## metaboliteFUMARIC                                                        
## metaboliteglutamine                                                      
## metaboliteisoleucine                                                     
## metaboliteLACTIC                                                         
## metaboliteLCAC total                                                     
## metaboliteleucine                                     0.526              
## metaboliteMALIC                                       0.513  0.513       
## metaboliteMCAC Total                                  0.526  0.526  0.513
## metaboliteMETHYLSUCCINIC                              0.513  0.513  0.500
## metabolitePYRUVIC_P2P                                 0.513  0.513  0.500
## metaboliteSUCCINIC-2                                  0.513  0.513  0.500
## metabolitevaline                                      0.526  0.526  0.513
## genotypeKO:activityExercise                          -0.411 -0.411 -0.401
## genotypeKO:metabolitearginine                        -0.377 -0.377 -0.368
## genotypeKO:metaboliteCITRIC                          -0.372 -0.372 -0.363
## genotypeKO:metaboliteFUMARIC                         -0.372 -0.372 -0.363
## genotypeKO:metaboliteglutamine                       -0.377 -0.377 -0.368
## genotypeKO:metaboliteisoleucine                      -0.377 -0.377 -0.368
## genotypeKO:metaboliteLACTIC                          -0.372 -0.372 -0.363
## genotypeKO:metaboliteLCAC total                      -0.717 -0.377 -0.368
## genotypeKO:metaboliteleucine                         -0.377 -0.717 -0.368
## genotypeKO:metaboliteMALIC                           -0.372 -0.372 -0.725
## genotypeKO:metaboliteMCAC Total                      -0.377 -0.377 -0.368
## genotypeKO:metaboliteMETHYLSUCCINIC                  -0.372 -0.372 -0.363
## genotypeKO:metabolitePYRUVIC_P2P                     -0.372 -0.372 -0.363
## genotypeKO:metaboliteSUCCINIC-2                      -0.372 -0.372 -0.363
## genotypeKO:metabolitevaline                          -0.377 -0.377 -0.368
## activityExercise:metabolitearginine                  -0.386 -0.386 -0.376
## activityExercise:metaboliteCITRIC                    -0.380 -0.380 -0.371
## activityExercise:metaboliteFUMARIC                   -0.380 -0.380 -0.371
## activityExercise:metaboliteglutamine                 -0.386 -0.386 -0.376
## activityExercise:metaboliteisoleucine                -0.386 -0.386 -0.376
## activityExercise:metaboliteLACTIC                    -0.380 -0.380 -0.371
## activityExercise:metaboliteLCAC total                -0.733 -0.386 -0.376
## activityExercise:metaboliteleucine                   -0.386 -0.733 -0.376
## activityExercise:metaboliteMALIC                     -0.380 -0.380 -0.742
## activityExercise:metaboliteMCAC Total                -0.386 -0.386 -0.376
## activityExercise:metaboliteMETHYLSUCCINIC            -0.380 -0.380 -0.371
## activityExercise:metabolitePYRUVIC_P2P               -0.380 -0.380 -0.371
## activityExercise:metaboliteSUCCINIC-2                -0.380 -0.380 -0.371
## activityExercise:metabolitevaline                    -0.386 -0.386 -0.376
## genotypeKO:activityExercise:metabolitearginine        0.306  0.306  0.298
## genotypeKO:activityExercise:metaboliteCITRIC          0.290  0.290  0.283
## genotypeKO:activityExercise:metaboliteFUMARIC         0.292  0.292  0.285
## genotypeKO:activityExercise:metaboliteglutamine       0.295  0.295  0.288
## genotypeKO:activityExercise:metaboliteisoleucine      0.294  0.294  0.287
## genotypeKO:activityExercise:metaboliteLACTIC          0.293  0.293  0.286
## genotypeKO:activityExercise:metaboliteLCAC total      0.552  0.290  0.283
## genotypeKO:activityExercise:metaboliteleucine         0.292  0.555  0.285
## genotypeKO:activityExercise:metaboliteMALIC           0.290  0.290  0.566
## genotypeKO:activityExercise:metaboliteMCAC Total      0.290  0.290  0.283
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC  0.292  0.292  0.285
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P     0.286  0.286  0.279
## genotypeKO:activityExercise:metaboliteSUCCINIC-2      0.285  0.285  0.278
## genotypeKO:activityExercise:metabolitevaline          0.293  0.293  0.285
##                                                      mMCACT mMETHY mPYRUV
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
## metaboliteMETHYLSUCCINIC                              0.513              
## metabolitePYRUVIC_P2P                                 0.513  0.500       
## metaboliteSUCCINIC-2                                  0.513  0.500  0.500
## metabolitevaline                                      0.526  0.513  0.513
## genotypeKO:activityExercise                          -0.411 -0.401 -0.401
## genotypeKO:metabolitearginine                        -0.377 -0.368 -0.368
## genotypeKO:metaboliteCITRIC                          -0.372 -0.363 -0.363
## genotypeKO:metaboliteFUMARIC                         -0.372 -0.363 -0.363
## genotypeKO:metaboliteglutamine                       -0.377 -0.368 -0.368
## genotypeKO:metaboliteisoleucine                      -0.377 -0.368 -0.368
## genotypeKO:metaboliteLACTIC                          -0.372 -0.363 -0.363
## genotypeKO:metaboliteLCAC total                      -0.377 -0.368 -0.368
## genotypeKO:metaboliteleucine                         -0.377 -0.368 -0.368
## genotypeKO:metaboliteMALIC                           -0.372 -0.363 -0.363
## genotypeKO:metaboliteMCAC Total                      -0.717 -0.368 -0.368
## genotypeKO:metaboliteMETHYLSUCCINIC                  -0.372 -0.725 -0.363
## genotypeKO:metabolitePYRUVIC_P2P                     -0.372 -0.363 -0.725
## genotypeKO:metaboliteSUCCINIC-2                      -0.372 -0.363 -0.363
## genotypeKO:metabolitevaline                          -0.377 -0.368 -0.368
## activityExercise:metabolitearginine                  -0.386 -0.376 -0.376
## activityExercise:metaboliteCITRIC                    -0.380 -0.371 -0.371
## activityExercise:metaboliteFUMARIC                   -0.380 -0.371 -0.371
## activityExercise:metaboliteglutamine                 -0.386 -0.376 -0.376
## activityExercise:metaboliteisoleucine                -0.386 -0.376 -0.376
## activityExercise:metaboliteLACTIC                    -0.380 -0.371 -0.371
## activityExercise:metaboliteLCAC total                -0.386 -0.376 -0.376
## activityExercise:metaboliteleucine                   -0.386 -0.376 -0.376
## activityExercise:metaboliteMALIC                     -0.380 -0.371 -0.371
## activityExercise:metaboliteMCAC Total                -0.733 -0.376 -0.376
## activityExercise:metaboliteMETHYLSUCCINIC            -0.380 -0.742 -0.371
## activityExercise:metabolitePYRUVIC_P2P               -0.380 -0.371 -0.742
## activityExercise:metaboliteSUCCINIC-2                -0.380 -0.371 -0.371
## activityExercise:metabolitevaline                    -0.386 -0.376 -0.376
## genotypeKO:activityExercise:metabolitearginine        0.306  0.298  0.298
## genotypeKO:activityExercise:metaboliteCITRIC          0.290  0.283  0.283
## genotypeKO:activityExercise:metaboliteFUMARIC         0.292  0.285  0.285
## genotypeKO:activityExercise:metaboliteglutamine       0.295  0.288  0.288
## genotypeKO:activityExercise:metaboliteisoleucine      0.294  0.287  0.287
## genotypeKO:activityExercise:metaboliteLACTIC          0.293  0.286  0.286
## genotypeKO:activityExercise:metaboliteLCAC total      0.290  0.283  0.283
## genotypeKO:activityExercise:metaboliteleucine         0.292  0.285  0.285
## genotypeKO:activityExercise:metaboliteMALIC           0.290  0.283  0.283
## genotypeKO:activityExercise:metaboliteMCAC Total      0.551  0.283  0.283
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC  0.292  0.570  0.285
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P     0.286  0.279  0.558
## genotypeKO:activityExercise:metaboliteSUCCINIC-2      0.285  0.278  0.278
## genotypeKO:activityExercise:metabolitevaline          0.293  0.285  0.285
##                                                      mSUCCI mtbltv gnKO:E
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
## metabolitevaline                                      0.513              
## genotypeKO:activityExercise                          -0.401 -0.411       
## genotypeKO:metabolitearginine                        -0.368 -0.377  0.560
## genotypeKO:metaboliteCITRIC                          -0.363 -0.372  0.552
## genotypeKO:metaboliteFUMARIC                         -0.363 -0.372  0.552
## genotypeKO:metaboliteglutamine                       -0.368 -0.377  0.560
## genotypeKO:metaboliteisoleucine                      -0.368 -0.377  0.560
## genotypeKO:metaboliteLACTIC                          -0.363 -0.372  0.552
## genotypeKO:metaboliteLCAC total                      -0.368 -0.377  0.560
## genotypeKO:metaboliteleucine                         -0.368 -0.377  0.560
## genotypeKO:metaboliteMALIC                           -0.363 -0.372  0.552
## genotypeKO:metaboliteMCAC Total                      -0.368 -0.377  0.560
## genotypeKO:metaboliteMETHYLSUCCINIC                  -0.363 -0.372  0.552
## genotypeKO:metabolitePYRUVIC_P2P                     -0.363 -0.372  0.552
## genotypeKO:metaboliteSUCCINIC-2                      -0.725 -0.372  0.552
## genotypeKO:metabolitevaline                          -0.368 -0.717  0.560
## activityExercise:metabolitearginine                  -0.376 -0.386  0.548
## activityExercise:metaboliteCITRIC                    -0.371 -0.380  0.540
## activityExercise:metaboliteFUMARIC                   -0.371 -0.380  0.540
## activityExercise:metaboliteglutamine                 -0.376 -0.386  0.548
## activityExercise:metaboliteisoleucine                -0.376 -0.386  0.548
## activityExercise:metaboliteLACTIC                    -0.371 -0.380  0.540
## activityExercise:metaboliteLCAC total                -0.376 -0.386  0.548
## activityExercise:metaboliteleucine                   -0.376 -0.386  0.548
## activityExercise:metaboliteMALIC                     -0.371 -0.380  0.540
## activityExercise:metaboliteMCAC Total                -0.376 -0.386  0.548
## activityExercise:metaboliteMETHYLSUCCINIC            -0.371 -0.380  0.540
## activityExercise:metabolitePYRUVIC_P2P               -0.371 -0.380  0.540
## activityExercise:metaboliteSUCCINIC-2                -0.742 -0.380  0.540
## activityExercise:metabolitevaline                    -0.376 -0.733  0.548
## genotypeKO:activityExercise:metabolitearginine        0.298  0.306 -0.693
## genotypeKO:activityExercise:metaboliteCITRIC          0.283  0.290 -0.700
## genotypeKO:activityExercise:metaboliteFUMARIC         0.285  0.292 -0.695
## genotypeKO:activityExercise:metaboliteglutamine       0.288  0.295 -0.694
## genotypeKO:activityExercise:metaboliteisoleucine      0.287  0.294 -0.728
## genotypeKO:activityExercise:metaboliteLACTIC          0.286  0.293 -0.693
## genotypeKO:activityExercise:metaboliteLCAC total      0.283  0.290 -0.716
## genotypeKO:activityExercise:metaboliteleucine         0.285  0.292 -0.721
## genotypeKO:activityExercise:metaboliteMALIC           0.283  0.290 -0.697
## genotypeKO:activityExercise:metaboliteMCAC Total      0.283  0.290 -0.709
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC  0.285  0.292 -0.694
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P     0.279  0.286 -0.702
## genotypeKO:activityExercise:metaboliteSUCCINIC-2      0.555  0.285 -0.713
## genotypeKO:activityExercise:metabolitevaline          0.285  0.556 -0.690
##                                                      gntypKO:mtbltr gKO:CI
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
## genotypeKO:metaboliteCITRIC                           0.507               
## genotypeKO:metaboliteFUMARIC                          0.507          0.500
## genotypeKO:metaboliteglutamine                        0.514          0.507
## genotypeKO:metaboliteisoleucine                       0.514          0.507
## genotypeKO:metaboliteLACTIC                           0.507          0.500
## genotypeKO:metaboliteLCAC total                       0.514          0.507
## genotypeKO:metaboliteleucine                          0.514          0.507
## genotypeKO:metaboliteMALIC                            0.507          0.500
## genotypeKO:metaboliteMCAC Total                       0.514          0.507
## genotypeKO:metaboliteMETHYLSUCCINIC                   0.507          0.500
## genotypeKO:metabolitePYRUVIC_P2P                      0.507          0.500
## genotypeKO:metaboliteSUCCINIC-2                       0.507          0.500
## genotypeKO:metabolitevaline                           0.514          0.507
## activityExercise:metabolitearginine                   0.525          0.273
## activityExercise:metaboliteCITRIC                     0.273          0.538
## activityExercise:metaboliteFUMARIC                    0.273          0.269
## activityExercise:metaboliteglutamine                  0.277          0.273
## activityExercise:metaboliteisoleucine                 0.277          0.273
## activityExercise:metaboliteLACTIC                     0.273          0.269
## activityExercise:metaboliteLCAC total                 0.277          0.273
## activityExercise:metaboliteleucine                    0.277          0.273
## activityExercise:metaboliteMALIC                      0.273          0.269
## activityExercise:metaboliteMCAC Total                 0.277          0.273
## activityExercise:metaboliteMETHYLSUCCINIC             0.273          0.269
## activityExercise:metabolitePYRUVIC_P2P                0.273          0.269
## activityExercise:metaboliteSUCCINIC-2                 0.273          0.269
## activityExercise:metabolitevaline                     0.277          0.273
## genotypeKO:activityExercise:metabolitearginine       -0.811         -0.411
## genotypeKO:activityExercise:metaboliteCITRIC         -0.395         -0.780
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.397         -0.392
## genotypeKO:activityExercise:metaboliteglutamine      -0.402         -0.397
## genotypeKO:activityExercise:metaboliteisoleucine     -0.401         -0.395
## genotypeKO:activityExercise:metaboliteLACTIC         -0.399         -0.394
## genotypeKO:activityExercise:metaboliteLCAC total     -0.395         -0.390
## genotypeKO:activityExercise:metaboliteleucine        -0.398         -0.392
## genotypeKO:activityExercise:metaboliteMALIC          -0.395         -0.390
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.395         -0.390
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.398         -0.393
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.389         -0.384
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.388         -0.383
## genotypeKO:activityExercise:metabolitevaline         -0.399         -0.393
##                                                      gKO:FU gntypKO:mtbltg
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
## genotypeKO:metaboliteglutamine                        0.507               
## genotypeKO:metaboliteisoleucine                       0.507  0.514        
## genotypeKO:metaboliteLACTIC                           0.500  0.507        
## genotypeKO:metaboliteLCAC total                       0.507  0.514        
## genotypeKO:metaboliteleucine                          0.507  0.514        
## genotypeKO:metaboliteMALIC                            0.500  0.507        
## genotypeKO:metaboliteMCAC Total                       0.507  0.514        
## genotypeKO:metaboliteMETHYLSUCCINIC                   0.500  0.507        
## genotypeKO:metabolitePYRUVIC_P2P                      0.500  0.507        
## genotypeKO:metaboliteSUCCINIC-2                       0.500  0.507        
## genotypeKO:metabolitevaline                           0.507  0.514        
## activityExercise:metabolitearginine                   0.273  0.277        
## activityExercise:metaboliteCITRIC                     0.269  0.273        
## activityExercise:metaboliteFUMARIC                    0.538  0.273        
## activityExercise:metaboliteglutamine                  0.273  0.525        
## activityExercise:metaboliteisoleucine                 0.273  0.277        
## activityExercise:metaboliteLACTIC                     0.269  0.273        
## activityExercise:metaboliteLCAC total                 0.273  0.277        
## activityExercise:metaboliteleucine                    0.273  0.277        
## activityExercise:metaboliteMALIC                      0.269  0.273        
## activityExercise:metaboliteMCAC Total                 0.273  0.277        
## activityExercise:metaboliteMETHYLSUCCINIC             0.269  0.273        
## activityExercise:metabolitePYRUVIC_P2P                0.269  0.273        
## activityExercise:metaboliteSUCCINIC-2                 0.269  0.273        
## activityExercise:metabolitevaline                     0.273  0.277        
## genotypeKO:activityExercise:metabolitearginine       -0.411 -0.416        
## genotypeKO:activityExercise:metaboliteCITRIC         -0.390 -0.395        
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.784 -0.397        
## genotypeKO:activityExercise:metaboliteglutamine      -0.397 -0.783        
## genotypeKO:activityExercise:metaboliteisoleucine     -0.395 -0.401        
## genotypeKO:activityExercise:metaboliteLACTIC         -0.394 -0.399        
## genotypeKO:activityExercise:metaboliteLCAC total     -0.390 -0.395        
## genotypeKO:activityExercise:metaboliteleucine        -0.392 -0.398        
## genotypeKO:activityExercise:metaboliteMALIC          -0.390 -0.395        
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.390 -0.395        
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.393 -0.398        
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.384 -0.389        
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.383 -0.388        
## genotypeKO:activityExercise:metabolitevaline         -0.393 -0.399        
##                                                      gntypKO:mtblts gKO:LA
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
## genotypeKO:metaboliteLACTIC                           0.507               
## genotypeKO:metaboliteLCAC total                       0.514          0.507
## genotypeKO:metaboliteleucine                          0.514          0.507
## genotypeKO:metaboliteMALIC                            0.507          0.500
## genotypeKO:metaboliteMCAC Total                       0.514          0.507
## genotypeKO:metaboliteMETHYLSUCCINIC                   0.507          0.500
## genotypeKO:metabolitePYRUVIC_P2P                      0.507          0.500
## genotypeKO:metaboliteSUCCINIC-2                       0.507          0.500
## genotypeKO:metabolitevaline                           0.514          0.507
## activityExercise:metabolitearginine                   0.277          0.273
## activityExercise:metaboliteCITRIC                     0.273          0.269
## activityExercise:metaboliteFUMARIC                    0.273          0.269
## activityExercise:metaboliteglutamine                  0.277          0.273
## activityExercise:metaboliteisoleucine                 0.525          0.273
## activityExercise:metaboliteLACTIC                     0.273          0.538
## activityExercise:metaboliteLCAC total                 0.277          0.273
## activityExercise:metaboliteleucine                    0.277          0.273
## activityExercise:metaboliteMALIC                      0.273          0.269
## activityExercise:metaboliteMCAC Total                 0.277          0.273
## activityExercise:metaboliteMETHYLSUCCINIC             0.273          0.269
## activityExercise:metabolitePYRUVIC_P2P                0.273          0.269
## activityExercise:metaboliteSUCCINIC-2                 0.273          0.269
## activityExercise:metabolitevaline                     0.277          0.273
## genotypeKO:activityExercise:metabolitearginine       -0.416         -0.411
## genotypeKO:activityExercise:metaboliteCITRIC         -0.395         -0.390
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.397         -0.392
## genotypeKO:activityExercise:metaboliteglutamine      -0.402         -0.397
## genotypeKO:activityExercise:metaboliteisoleucine     -0.780         -0.395
## genotypeKO:activityExercise:metaboliteLACTIC         -0.399         -0.788
## genotypeKO:activityExercise:metaboliteLCAC total     -0.395         -0.390
## genotypeKO:activityExercise:metaboliteleucine        -0.398         -0.392
## genotypeKO:activityExercise:metaboliteMALIC          -0.395         -0.390
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.395         -0.390
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.398         -0.393
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.389         -0.384
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.388         -0.383
## genotypeKO:activityExercise:metabolitevaline         -0.399         -0.393
##                                                      gKO:Lt gntypKO:mtbltl
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
## genotypeKO:metaboliteleucine                          0.514               
## genotypeKO:metaboliteMALIC                            0.507  0.507        
## genotypeKO:metaboliteMCAC Total                       0.514  0.514        
## genotypeKO:metaboliteMETHYLSUCCINIC                   0.507  0.507        
## genotypeKO:metabolitePYRUVIC_P2P                      0.507  0.507        
## genotypeKO:metaboliteSUCCINIC-2                       0.507  0.507        
## genotypeKO:metabolitevaline                           0.514  0.514        
## activityExercise:metabolitearginine                   0.277  0.277        
## activityExercise:metaboliteCITRIC                     0.273  0.273        
## activityExercise:metaboliteFUMARIC                    0.273  0.273        
## activityExercise:metaboliteglutamine                  0.277  0.277        
## activityExercise:metaboliteisoleucine                 0.277  0.277        
## activityExercise:metaboliteLACTIC                     0.273  0.273        
## activityExercise:metaboliteLCAC total                 0.525  0.277        
## activityExercise:metaboliteleucine                    0.277  0.525        
## activityExercise:metaboliteMALIC                      0.273  0.273        
## activityExercise:metaboliteMCAC Total                 0.277  0.277        
## activityExercise:metaboliteMETHYLSUCCINIC             0.273  0.273        
## activityExercise:metabolitePYRUVIC_P2P                0.273  0.273        
## activityExercise:metaboliteSUCCINIC-2                 0.273  0.273        
## activityExercise:metabolitevaline                     0.277  0.277        
## genotypeKO:activityExercise:metabolitearginine       -0.416 -0.416        
## genotypeKO:activityExercise:metaboliteCITRIC         -0.395 -0.395        
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.397 -0.397        
## genotypeKO:activityExercise:metaboliteglutamine      -0.402 -0.402        
## genotypeKO:activityExercise:metaboliteisoleucine     -0.401 -0.401        
## genotypeKO:activityExercise:metaboliteLACTIC         -0.399 -0.399        
## genotypeKO:activityExercise:metaboliteLCAC total     -0.770 -0.395        
## genotypeKO:activityExercise:metaboliteleucine        -0.398 -0.775        
## genotypeKO:activityExercise:metaboliteMALIC          -0.395 -0.395        
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.395 -0.395        
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.398 -0.398        
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.389 -0.389        
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.388 -0.388        
## genotypeKO:activityExercise:metabolitevaline         -0.399 -0.399        
##                                                      gKO:MA gKO:MT gKO:ME
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
## genotypeKO:metaboliteMCAC Total                       0.507              
## genotypeKO:metaboliteMETHYLSUCCINIC                   0.500  0.507       
## genotypeKO:metabolitePYRUVIC_P2P                      0.500  0.507  0.500
## genotypeKO:metaboliteSUCCINIC-2                       0.500  0.507  0.500
## genotypeKO:metabolitevaline                           0.507  0.514  0.507
## activityExercise:metabolitearginine                   0.273  0.277  0.273
## activityExercise:metaboliteCITRIC                     0.269  0.273  0.269
## activityExercise:metaboliteFUMARIC                    0.269  0.273  0.269
## activityExercise:metaboliteglutamine                  0.273  0.277  0.273
## activityExercise:metaboliteisoleucine                 0.273  0.277  0.273
## activityExercise:metaboliteLACTIC                     0.269  0.273  0.269
## activityExercise:metaboliteLCAC total                 0.273  0.277  0.273
## activityExercise:metaboliteleucine                    0.273  0.277  0.273
## activityExercise:metaboliteMALIC                      0.538  0.273  0.269
## activityExercise:metaboliteMCAC Total                 0.273  0.525  0.273
## activityExercise:metaboliteMETHYLSUCCINIC             0.269  0.273  0.538
## activityExercise:metabolitePYRUVIC_P2P                0.269  0.273  0.269
## activityExercise:metaboliteSUCCINIC-2                 0.269  0.273  0.269
## activityExercise:metabolitevaline                     0.273  0.277  0.273
## genotypeKO:activityExercise:metabolitearginine       -0.411 -0.416 -0.411
## genotypeKO:activityExercise:metaboliteCITRIC         -0.390 -0.395 -0.390
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.392 -0.397 -0.392
## genotypeKO:activityExercise:metaboliteglutamine      -0.397 -0.402 -0.397
## genotypeKO:activityExercise:metaboliteisoleucine     -0.395 -0.401 -0.395
## genotypeKO:activityExercise:metaboliteLACTIC         -0.394 -0.399 -0.394
## genotypeKO:activityExercise:metaboliteLCAC total     -0.390 -0.395 -0.390
## genotypeKO:activityExercise:metaboliteleucine        -0.392 -0.398 -0.392
## genotypeKO:activityExercise:metaboliteMALIC          -0.780 -0.395 -0.390
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.390 -0.769 -0.390
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.393 -0.398 -0.786
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.384 -0.389 -0.384
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.383 -0.388 -0.383
## genotypeKO:activityExercise:metabolitevaline         -0.393 -0.399 -0.393
##                                                      gKO:PY gKO:SU
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
## genotypeKO:metaboliteSUCCINIC-2                       0.500       
## genotypeKO:metabolitevaline                           0.507  0.507
## activityExercise:metabolitearginine                   0.273  0.273
## activityExercise:metaboliteCITRIC                     0.269  0.269
## activityExercise:metaboliteFUMARIC                    0.269  0.269
## activityExercise:metaboliteglutamine                  0.273  0.273
## activityExercise:metaboliteisoleucine                 0.273  0.273
## activityExercise:metaboliteLACTIC                     0.269  0.269
## activityExercise:metaboliteLCAC total                 0.273  0.273
## activityExercise:metaboliteleucine                    0.273  0.273
## activityExercise:metaboliteMALIC                      0.269  0.269
## activityExercise:metaboliteMCAC Total                 0.273  0.273
## activityExercise:metaboliteMETHYLSUCCINIC             0.269  0.269
## activityExercise:metabolitePYRUVIC_P2P                0.538  0.269
## activityExercise:metaboliteSUCCINIC-2                 0.269  0.538
## activityExercise:metabolitevaline                     0.273  0.273
## genotypeKO:activityExercise:metabolitearginine       -0.411 -0.411
## genotypeKO:activityExercise:metaboliteCITRIC         -0.390 -0.390
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.392 -0.392
## genotypeKO:activityExercise:metaboliteglutamine      -0.397 -0.397
## genotypeKO:activityExercise:metaboliteisoleucine     -0.395 -0.395
## genotypeKO:activityExercise:metaboliteLACTIC         -0.394 -0.394
## genotypeKO:activityExercise:metaboliteLCAC total     -0.390 -0.390
## genotypeKO:activityExercise:metaboliteleucine        -0.392 -0.392
## genotypeKO:activityExercise:metaboliteMALIC          -0.390 -0.390
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.390 -0.390
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.393 -0.393
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.769 -0.384
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.383 -0.765
## genotypeKO:activityExercise:metabolitevaline         -0.393 -0.393
##                                                      gntypKO:mtbltv
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
## activityExercise:metabolitearginine                   0.277        
## activityExercise:metaboliteCITRIC                     0.273        
## activityExercise:metaboliteFUMARIC                    0.273        
## activityExercise:metaboliteglutamine                  0.277        
## activityExercise:metaboliteisoleucine                 0.277        
## activityExercise:metaboliteLACTIC                     0.273        
## activityExercise:metaboliteLCAC total                 0.277        
## activityExercise:metaboliteleucine                    0.277        
## activityExercise:metaboliteMALIC                      0.273        
## activityExercise:metaboliteMCAC Total                 0.277        
## activityExercise:metaboliteMETHYLSUCCINIC             0.273        
## activityExercise:metabolitePYRUVIC_P2P                0.273        
## activityExercise:metaboliteSUCCINIC-2                 0.273        
## activityExercise:metabolitevaline                     0.525        
## genotypeKO:activityExercise:metabolitearginine       -0.416        
## genotypeKO:activityExercise:metaboliteCITRIC         -0.395        
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.397        
## genotypeKO:activityExercise:metaboliteglutamine      -0.402        
## genotypeKO:activityExercise:metaboliteisoleucine     -0.401        
## genotypeKO:activityExercise:metaboliteLACTIC         -0.399        
## genotypeKO:activityExercise:metaboliteLCAC total     -0.395        
## genotypeKO:activityExercise:metaboliteleucine        -0.398        
## genotypeKO:activityExercise:metaboliteMALIC          -0.395        
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.395        
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.398        
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.389        
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.388        
## genotypeKO:activityExercise:metabolitevaline         -0.776        
##                                                      actvtyExrcs:mtbltr
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
## activityExercise:metaboliteCITRIC                     0.507            
## activityExercise:metaboliteFUMARIC                    0.507            
## activityExercise:metaboliteglutamine                  0.514            
## activityExercise:metaboliteisoleucine                 0.514            
## activityExercise:metaboliteLACTIC                     0.507            
## activityExercise:metaboliteLCAC total                 0.514            
## activityExercise:metaboliteleucine                    0.514            
## activityExercise:metaboliteMALIC                      0.507            
## activityExercise:metaboliteMCAC Total                 0.514            
## activityExercise:metaboliteMETHYLSUCCINIC             0.507            
## activityExercise:metabolitePYRUVIC_P2P                0.507            
## activityExercise:metaboliteSUCCINIC-2                 0.507            
## activityExercise:metabolitevaline                     0.514            
## genotypeKO:activityExercise:metabolitearginine       -0.793            
## genotypeKO:activityExercise:metaboliteCITRIC         -0.387            
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.389            
## genotypeKO:activityExercise:metaboliteglutamine      -0.394            
## genotypeKO:activityExercise:metaboliteisoleucine     -0.392            
## genotypeKO:activityExercise:metaboliteLACTIC         -0.391            
## genotypeKO:activityExercise:metaboliteLCAC total     -0.387            
## genotypeKO:activityExercise:metaboliteleucine        -0.389            
## genotypeKO:activityExercise:metaboliteMALIC          -0.387            
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.387            
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.390            
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.381            
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.379            
## genotypeKO:activityExercise:metabolitevaline         -0.390            
##                                                      aE:CIT aE:FUM
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
## activityExercise:metaboliteFUMARIC                    0.500       
## activityExercise:metaboliteglutamine                  0.507  0.507
## activityExercise:metaboliteisoleucine                 0.507  0.507
## activityExercise:metaboliteLACTIC                     0.500  0.500
## activityExercise:metaboliteLCAC total                 0.507  0.507
## activityExercise:metaboliteleucine                    0.507  0.507
## activityExercise:metaboliteMALIC                      0.500  0.500
## activityExercise:metaboliteMCAC Total                 0.507  0.507
## activityExercise:metaboliteMETHYLSUCCINIC             0.500  0.500
## activityExercise:metabolitePYRUVIC_P2P                0.500  0.500
## activityExercise:metaboliteSUCCINIC-2                 0.500  0.500
## activityExercise:metabolitevaline                     0.507  0.507
## genotypeKO:activityExercise:metabolitearginine       -0.402 -0.402
## genotypeKO:activityExercise:metaboliteCITRIC         -0.763 -0.381
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.384 -0.767
## genotypeKO:activityExercise:metaboliteglutamine      -0.388 -0.388
## genotypeKO:activityExercise:metaboliteisoleucine     -0.387 -0.387
## genotypeKO:activityExercise:metaboliteLACTIC         -0.385 -0.385
## genotypeKO:activityExercise:metaboliteLCAC total     -0.382 -0.382
## genotypeKO:activityExercise:metaboliteleucine        -0.384 -0.384
## genotypeKO:activityExercise:metaboliteMALIC          -0.382 -0.382
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.381 -0.381
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.384 -0.384
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.376 -0.376
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.374 -0.374
## genotypeKO:activityExercise:metabolitevaline         -0.385 -0.385
##                                                      actvtyExrcs:mtbltg
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
## activityExercise:metaboliteisoleucine                 0.514            
## activityExercise:metaboliteLACTIC                     0.507            
## activityExercise:metaboliteLCAC total                 0.514            
## activityExercise:metaboliteleucine                    0.514            
## activityExercise:metaboliteMALIC                      0.507            
## activityExercise:metaboliteMCAC Total                 0.514            
## activityExercise:metaboliteMETHYLSUCCINIC             0.507            
## activityExercise:metabolitePYRUVIC_P2P                0.507            
## activityExercise:metaboliteSUCCINIC-2                 0.507            
## activityExercise:metabolitevaline                     0.514            
## genotypeKO:activityExercise:metabolitearginine       -0.408            
## genotypeKO:activityExercise:metaboliteCITRIC         -0.387            
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.389            
## genotypeKO:activityExercise:metaboliteglutamine      -0.765            
## genotypeKO:activityExercise:metaboliteisoleucine     -0.392            
## genotypeKO:activityExercise:metaboliteLACTIC         -0.391            
## genotypeKO:activityExercise:metaboliteLCAC total     -0.387            
## genotypeKO:activityExercise:metaboliteleucine        -0.389            
## genotypeKO:activityExercise:metaboliteMALIC          -0.387            
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.387            
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.390            
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.381            
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.379            
## genotypeKO:activityExercise:metabolitevaline         -0.390            
##                                                      actvtyExrcs:mtblts
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
## activityExercise:metaboliteLACTIC                     0.507            
## activityExercise:metaboliteLCAC total                 0.514            
## activityExercise:metaboliteleucine                    0.514            
## activityExercise:metaboliteMALIC                      0.507            
## activityExercise:metaboliteMCAC Total                 0.514            
## activityExercise:metaboliteMETHYLSUCCINIC             0.507            
## activityExercise:metabolitePYRUVIC_P2P                0.507            
## activityExercise:metaboliteSUCCINIC-2                 0.507            
## activityExercise:metabolitevaline                     0.514            
## genotypeKO:activityExercise:metabolitearginine       -0.408            
## genotypeKO:activityExercise:metaboliteCITRIC         -0.387            
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.389            
## genotypeKO:activityExercise:metaboliteglutamine      -0.394            
## genotypeKO:activityExercise:metaboliteisoleucine     -0.762            
## genotypeKO:activityExercise:metaboliteLACTIC         -0.391            
## genotypeKO:activityExercise:metaboliteLCAC total     -0.387            
## genotypeKO:activityExercise:metaboliteleucine        -0.389            
## genotypeKO:activityExercise:metaboliteMALIC          -0.387            
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.387            
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.390            
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.381            
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.379            
## genotypeKO:activityExercise:metabolitevaline         -0.390            
##                                                      aE:LAC aE:LCt
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
## activityExercise:metaboliteLCAC total                 0.507       
## activityExercise:metaboliteleucine                    0.507  0.514
## activityExercise:metaboliteMALIC                      0.500  0.507
## activityExercise:metaboliteMCAC Total                 0.507  0.514
## activityExercise:metaboliteMETHYLSUCCINIC             0.500  0.507
## activityExercise:metabolitePYRUVIC_P2P                0.500  0.507
## activityExercise:metaboliteSUCCINIC-2                 0.500  0.507
## activityExercise:metabolitevaline                     0.507  0.514
## genotypeKO:activityExercise:metabolitearginine       -0.402 -0.408
## genotypeKO:activityExercise:metaboliteCITRIC         -0.381 -0.387
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.384 -0.389
## genotypeKO:activityExercise:metaboliteglutamine      -0.388 -0.394
## genotypeKO:activityExercise:metaboliteisoleucine     -0.387 -0.392
## genotypeKO:activityExercise:metaboliteLACTIC         -0.771 -0.391
## genotypeKO:activityExercise:metaboliteLCAC total     -0.382 -0.753
## genotypeKO:activityExercise:metaboliteleucine        -0.384 -0.389
## genotypeKO:activityExercise:metaboliteMALIC          -0.382 -0.387
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.381 -0.387
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.384 -0.390
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.376 -0.381
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.374 -0.379
## genotypeKO:activityExercise:metabolitevaline         -0.385 -0.390
##                                                      actvtyExrcs:mtbltl
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
## activityExercise:metaboliteMALIC                      0.507            
## activityExercise:metaboliteMCAC Total                 0.514            
## activityExercise:metaboliteMETHYLSUCCINIC             0.507            
## activityExercise:metabolitePYRUVIC_P2P                0.507            
## activityExercise:metaboliteSUCCINIC-2                 0.507            
## activityExercise:metabolitevaline                     0.514            
## genotypeKO:activityExercise:metabolitearginine       -0.408            
## genotypeKO:activityExercise:metaboliteCITRIC         -0.387            
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.389            
## genotypeKO:activityExercise:metaboliteglutamine      -0.394            
## genotypeKO:activityExercise:metaboliteisoleucine     -0.392            
## genotypeKO:activityExercise:metaboliteLACTIC         -0.391            
## genotypeKO:activityExercise:metaboliteLCAC total     -0.387            
## genotypeKO:activityExercise:metaboliteleucine        -0.757            
## genotypeKO:activityExercise:metaboliteMALIC          -0.387            
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.387            
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.390            
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.381            
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.379            
## genotypeKO:activityExercise:metabolitevaline         -0.390            
##                                                      aE:MAL aE:MCT aE:MET
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
## activityExercise:metaboliteMCAC Total                 0.507              
## activityExercise:metaboliteMETHYLSUCCINIC             0.500  0.507       
## activityExercise:metabolitePYRUVIC_P2P                0.500  0.507  0.500
## activityExercise:metaboliteSUCCINIC-2                 0.500  0.507  0.500
## activityExercise:metabolitevaline                     0.507  0.514  0.507
## genotypeKO:activityExercise:metabolitearginine       -0.402 -0.408 -0.402
## genotypeKO:activityExercise:metaboliteCITRIC         -0.381 -0.387 -0.381
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.384 -0.389 -0.384
## genotypeKO:activityExercise:metaboliteglutamine      -0.388 -0.394 -0.388
## genotypeKO:activityExercise:metaboliteisoleucine     -0.387 -0.392 -0.387
## genotypeKO:activityExercise:metaboliteLACTIC         -0.385 -0.391 -0.385
## genotypeKO:activityExercise:metaboliteLCAC total     -0.382 -0.387 -0.382
## genotypeKO:activityExercise:metaboliteleucine        -0.384 -0.389 -0.384
## genotypeKO:activityExercise:metaboliteMALIC          -0.763 -0.387 -0.382
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.381 -0.752 -0.381
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.384 -0.390 -0.769
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.376 -0.381 -0.376
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.374 -0.379 -0.374
## genotypeKO:activityExercise:metabolitevaline         -0.385 -0.390 -0.385
##                                                      aE:PYR aE:SUC
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
## activityExercise:metabolitePYRUVIC_P2P                            
## activityExercise:metaboliteSUCCINIC-2                 0.500       
## activityExercise:metabolitevaline                     0.507  0.507
## genotypeKO:activityExercise:metabolitearginine       -0.402 -0.402
## genotypeKO:activityExercise:metaboliteCITRIC         -0.381 -0.381
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.384 -0.384
## genotypeKO:activityExercise:metaboliteglutamine      -0.388 -0.388
## genotypeKO:activityExercise:metaboliteisoleucine     -0.387 -0.387
## genotypeKO:activityExercise:metaboliteLACTIC         -0.385 -0.385
## genotypeKO:activityExercise:metaboliteLCAC total     -0.382 -0.382
## genotypeKO:activityExercise:metaboliteleucine        -0.384 -0.384
## genotypeKO:activityExercise:metaboliteMALIC          -0.382 -0.382
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.381 -0.381
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.384 -0.384
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.752 -0.376
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.374 -0.748
## genotypeKO:activityExercise:metabolitevaline         -0.385 -0.385
##                                                      actvtyExrcs:mtbltv
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
## activityExercise:metabolitePYRUVIC_P2P                                 
## activityExercise:metaboliteSUCCINIC-2                                  
## activityExercise:metabolitevaline                                      
## genotypeKO:activityExercise:metabolitearginine       -0.408            
## genotypeKO:activityExercise:metaboliteCITRIC         -0.387            
## genotypeKO:activityExercise:metaboliteFUMARIC        -0.389            
## genotypeKO:activityExercise:metaboliteglutamine      -0.394            
## genotypeKO:activityExercise:metaboliteisoleucine     -0.392            
## genotypeKO:activityExercise:metaboliteLACTIC         -0.391            
## genotypeKO:activityExercise:metaboliteLCAC total     -0.387            
## genotypeKO:activityExercise:metaboliteleucine        -0.389            
## genotypeKO:activityExercise:metaboliteMALIC          -0.387            
## genotypeKO:activityExercise:metaboliteMCAC Total     -0.387            
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC -0.390            
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P    -0.381            
## genotypeKO:activityExercise:metaboliteSUCCINIC-2     -0.379            
## genotypeKO:activityExercise:metabolitevaline         -0.759            
##                                                      gntypKO:ctvtyExrcs:mtbltr
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
## activityExercise:metabolitePYRUVIC_P2P                                        
## activityExercise:metaboliteSUCCINIC-2                                         
## activityExercise:metabolitevaline                                             
## genotypeKO:activityExercise:metabolitearginine                                
## genotypeKO:activityExercise:metaboliteCITRIC          0.499                   
## genotypeKO:activityExercise:metaboliteFUMARIC         0.492                   
## genotypeKO:activityExercise:metaboliteglutamine       0.489                   
## genotypeKO:activityExercise:metaboliteisoleucine      0.508                   
## genotypeKO:activityExercise:metaboliteLACTIC          0.461                   
## genotypeKO:activityExercise:metaboliteLCAC total      0.482                   
## genotypeKO:activityExercise:metaboliteleucine         0.509                   
## genotypeKO:activityExercise:metaboliteMALIC           0.492                   
## genotypeKO:activityExercise:metaboliteMCAC Total      0.486                   
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC  0.474                   
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P     0.515                   
## genotypeKO:activityExercise:metaboliteSUCCINIC-2      0.486                   
## genotypeKO:activityExercise:metabolitevaline          0.500                   
##                                                      gKO:E:C gKO:E:F
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
## activityExercise:metabolitePYRUVIC_P2P                              
## activityExercise:metaboliteSUCCINIC-2                               
## activityExercise:metabolitevaline                                   
## genotypeKO:activityExercise:metabolitearginine                      
## genotypeKO:activityExercise:metaboliteCITRIC                        
## genotypeKO:activityExercise:metaboliteFUMARIC         0.475         
## genotypeKO:activityExercise:metaboliteglutamine       0.480   0.478 
## genotypeKO:activityExercise:metaboliteisoleucine      0.527   0.506 
## genotypeKO:activityExercise:metaboliteLACTIC          0.501   0.478 
## genotypeKO:activityExercise:metaboliteLCAC total      0.514   0.505 
## genotypeKO:activityExercise:metaboliteleucine         0.540   0.485 
## genotypeKO:activityExercise:metaboliteMALIC           0.511   0.504 
## genotypeKO:activityExercise:metaboliteMCAC Total      0.474   0.494 
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC  0.487   0.506 
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P     0.490   0.487 
## genotypeKO:activityExercise:metaboliteSUCCINIC-2      0.480   0.493 
## genotypeKO:activityExercise:metabolitevaline          0.497   0.465 
##                                                      gntypKO:ctvtyExrcs:mtbltg
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
## activityExercise:metabolitePYRUVIC_P2P                                        
## activityExercise:metaboliteSUCCINIC-2                                         
## activityExercise:metabolitevaline                                             
## genotypeKO:activityExercise:metabolitearginine                                
## genotypeKO:activityExercise:metaboliteCITRIC                                  
## genotypeKO:activityExercise:metaboliteFUMARIC                                 
## genotypeKO:activityExercise:metaboliteglutamine                               
## genotypeKO:activityExercise:metaboliteisoleucine      0.504                   
## genotypeKO:activityExercise:metaboliteLACTIC          0.510                   
## genotypeKO:activityExercise:metaboliteLCAC total      0.480                   
## genotypeKO:activityExercise:metaboliteleucine         0.500                   
## genotypeKO:activityExercise:metaboliteMALIC           0.490                   
## genotypeKO:activityExercise:metaboliteMCAC Total      0.496                   
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC  0.474                   
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P     0.521                   
## genotypeKO:activityExercise:metaboliteSUCCINIC-2      0.491                   
## genotypeKO:activityExercise:metabolitevaline          0.496                   
##                                                      gntypKO:ctvtyExrcs:mtblts
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
## activityExercise:metabolitePYRUVIC_P2P                                        
## activityExercise:metaboliteSUCCINIC-2                                         
## activityExercise:metabolitevaline                                             
## genotypeKO:activityExercise:metabolitearginine                                
## genotypeKO:activityExercise:metaboliteCITRIC                                  
## genotypeKO:activityExercise:metaboliteFUMARIC                                 
## genotypeKO:activityExercise:metaboliteglutamine                               
## genotypeKO:activityExercise:metaboliteisoleucine                              
## genotypeKO:activityExercise:metaboliteLACTIC          0.499                   
## genotypeKO:activityExercise:metaboliteLCAC total      0.537                   
## genotypeKO:activityExercise:metaboliteleucine         0.541                   
## genotypeKO:activityExercise:metaboliteMALIC           0.526                   
## genotypeKO:activityExercise:metaboliteMCAC Total      0.517                   
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC  0.512                   
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P     0.503                   
## genotypeKO:activityExercise:metaboliteSUCCINIC-2      0.539                   
## genotypeKO:activityExercise:metabolitevaline          0.511                   
##                                                      gKO:E:L gKO:Et
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
## activityExercise:metabolitePYRUVIC_P2P                             
## activityExercise:metaboliteSUCCINIC-2                              
## activityExercise:metabolitevaline                                  
## genotypeKO:activityExercise:metabolitearginine                     
## genotypeKO:activityExercise:metaboliteCITRIC                       
## genotypeKO:activityExercise:metaboliteFUMARIC                      
## genotypeKO:activityExercise:metaboliteglutamine                    
## genotypeKO:activityExercise:metaboliteisoleucine                   
## genotypeKO:activityExercise:metaboliteLACTIC                       
## genotypeKO:activityExercise:metaboliteLCAC total      0.502        
## genotypeKO:activityExercise:metaboliteleucine         0.508   0.526
## genotypeKO:activityExercise:metaboliteMALIC           0.498   0.484
## genotypeKO:activityExercise:metaboliteMCAC Total      0.524   0.510
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC  0.512   0.527
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P     0.499   0.496
## genotypeKO:activityExercise:metaboliteSUCCINIC-2      0.482   0.518
## genotypeKO:activityExercise:metabolitevaline          0.475   0.504
##                                                      gntypKO:ctvtyExrcs:mtbltl
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
## activityExercise:metabolitePYRUVIC_P2P                                        
## activityExercise:metaboliteSUCCINIC-2                                         
## activityExercise:metabolitevaline                                             
## genotypeKO:activityExercise:metabolitearginine                                
## genotypeKO:activityExercise:metaboliteCITRIC                                  
## genotypeKO:activityExercise:metaboliteFUMARIC                                 
## genotypeKO:activityExercise:metaboliteglutamine                               
## genotypeKO:activityExercise:metaboliteisoleucine                              
## genotypeKO:activityExercise:metaboliteLACTIC                                  
## genotypeKO:activityExercise:metaboliteLCAC total                              
## genotypeKO:activityExercise:metaboliteleucine                                 
## genotypeKO:activityExercise:metaboliteMALIC           0.499                   
## genotypeKO:activityExercise:metaboliteMCAC Total      0.527                   
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC  0.491                   
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P     0.506                   
## genotypeKO:activityExercise:metaboliteSUCCINIC-2      0.503                   
## genotypeKO:activityExercise:metabolitevaline          0.516                   
##                                                      gKO:E:MA gKO:ET
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
## activityExercise:metabolitePYRUVIC_P2P                              
## activityExercise:metaboliteSUCCINIC-2                               
## activityExercise:metabolitevaline                                   
## genotypeKO:activityExercise:metabolitearginine                      
## genotypeKO:activityExercise:metaboliteCITRIC                        
## genotypeKO:activityExercise:metaboliteFUMARIC                       
## genotypeKO:activityExercise:metaboliteglutamine                     
## genotypeKO:activityExercise:metaboliteisoleucine                    
## genotypeKO:activityExercise:metaboliteLACTIC                        
## genotypeKO:activityExercise:metaboliteLCAC total                    
## genotypeKO:activityExercise:metaboliteleucine                       
## genotypeKO:activityExercise:metaboliteMALIC                         
## genotypeKO:activityExercise:metaboliteMCAC Total      0.487         
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC  0.511    0.488
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P     0.476    0.519
## genotypeKO:activityExercise:metaboliteSUCCINIC-2      0.504    0.517
## genotypeKO:activityExercise:metabolitevaline          0.507    0.506
##                                                      gKO:E:ME gKO:E:P
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
## activityExercise:metabolitePYRUVIC_P2P                               
## activityExercise:metaboliteSUCCINIC-2                                
## activityExercise:metabolitevaline                                    
## genotypeKO:activityExercise:metabolitearginine                       
## genotypeKO:activityExercise:metaboliteCITRIC                         
## genotypeKO:activityExercise:metaboliteFUMARIC                        
## genotypeKO:activityExercise:metaboliteglutamine                      
## genotypeKO:activityExercise:metaboliteisoleucine                     
## genotypeKO:activityExercise:metaboliteLACTIC                         
## genotypeKO:activityExercise:metaboliteLCAC total                     
## genotypeKO:activityExercise:metaboliteleucine                        
## genotypeKO:activityExercise:metaboliteMALIC                          
## genotypeKO:activityExercise:metaboliteMCAC Total                     
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC                 
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P     0.480          
## genotypeKO:activityExercise:metaboliteSUCCINIC-2      0.506    0.506 
## genotypeKO:activityExercise:metabolitevaline          0.492    0.490 
##                                                      gKO:E:S
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
## activityExercise:metabolitePYRUVIC_P2P                      
## activityExercise:metaboliteSUCCINIC-2                       
## activityExercise:metabolitevaline                           
## genotypeKO:activityExercise:metabolitearginine              
## genotypeKO:activityExercise:metaboliteCITRIC                
## genotypeKO:activityExercise:metaboliteFUMARIC               
## genotypeKO:activityExercise:metaboliteglutamine             
## genotypeKO:activityExercise:metaboliteisoleucine            
## genotypeKO:activityExercise:metaboliteLACTIC                
## genotypeKO:activityExercise:metaboliteLCAC total            
## genotypeKO:activityExercise:metaboliteleucine               
## genotypeKO:activityExercise:metaboliteMALIC                 
## genotypeKO:activityExercise:metaboliteMCAC Total            
## genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC        
## genotypeKO:activityExercise:metabolitePYRUVIC_P2P           
## genotypeKO:activityExercise:metaboliteSUCCINIC-2            
## genotypeKO:activityExercise:metabolitevaline          0.501 
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -8.94348308 -0.08755239  0.01523379  0.12984547  3.96700365 
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
                   chow * metabolite +
                   chow * metabolite * genotype)
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
##                          numDF denDF      F-value p-value
## (Intercept)                  1   546     11108818  <.0001
## genotype                     1    39 858664163314  <.0001
## chow                         1    39  40260488735  <.0001
## metabolite                  14   546 289717015158  <.0001
## genotype:chow                1    39            2  0.2069
## genotype:metabolite         14   546            4  <.0001
## chow:metabolite             14   546            3  <.0001
## genotype:chow:metabolite    14   546            1  0.1949
```

```r
summary(M)
```

```
## Linear mixed-effects model fit by REML
##  Data: D2 
##        AIC      BIC    logLik
##   1409.398 2139.458 -537.6992
## 
## Random effects:
##  Formula: ~1 | id
##         (Intercept)  Residual
## StdDev:  0.07226382 0.5209843
## 
## Correlation Structure: General
##  Formula: ~1 | id 
##  Parameter estimate(s):
##  Correlation: 
##    1      2      3      4      5      6      7      8      9      10    
## 2   0.026                                                               
## 3   0.010  0.009                                                        
## 4  -0.004  0.032 -0.007                                                 
## 5   0.009  0.062  0.014  0.023                                          
## 6  -0.124 -0.011  0.025 -0.112 -0.003                                   
## 7  -0.039 -0.002  0.009  0.022  0.026 -0.001                            
## 8   0.039 -0.001  0.023 -0.039  0.019 -0.069  0.007                     
## 9   0.012  0.043 -0.002  0.011  0.033  0.000  0.004 -0.003              
## 10 -0.028 -0.005  0.021 -0.172 -0.012 -0.246  0.006 -0.045  0.009       
## 11 -0.053  0.029 -0.003 -0.003 -0.007  0.034  0.040 -0.006  0.009  0.043
## 12  0.062 -0.021 -0.001 -0.042  0.012 -0.133  0.025  0.079  0.018  0.011
## 13 -0.023 -0.001 -0.008 -0.257  0.000 -0.225 -0.002 -0.107 -0.004 -0.083
## 14 -0.358 -0.031  0.000  0.025 -0.096  0.226 -0.028 -0.021  0.015  0.033
## 15  0.091  0.022  0.092 -0.027  0.024  0.018  0.084  0.132  0.002 -0.072
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
## 12  0.052                     
## 13 -0.008 -0.090              
## 14  0.160 -0.016 -0.047       
## 15 -0.022  0.034 -0.005 -0.037
## Fixed effects: list(fixed) 
##                                                            Value Std.Error
## (Intercept)                                            0.2557123 0.1585866
## genotypeKO                                             0.6626366 0.2242753
## chowYellow (C8)                                        0.7549455 0.2009892
## metabolitearginine                                     0.1236579 0.2221485
## metaboliteCITRIC                                       0.0395968 0.2221485
## metaboliteFUMARIC                                     -1.0066586 0.2221485
## metaboliteglutamine                                    1.7844750 0.2221485
## metaboliteisoleucine                                  -0.0276526 0.2221485
## metaboliteLACTIC                                       2.1516434 0.2221485
## metaboliteLC even AC total                             0.4067190 0.2221485
## metaboliteLC odd AC total                             -0.5647188 0.2221485
## metaboliteleucine                                      0.0958975 0.2221485
## metaboliteMALIC                                        0.7454015 0.2221485
## metaboliteMCAC total                                  -0.7176324 0.2221485
## metaboliteMETHYLSUCCINIC                              -1.4975941 0.2221485
## metaboliteSUCCINIC-2                                   0.9674087 0.2221485
## metabolitevaline                                       0.1464322 0.2221485
## genotypeKO:chowYellow (C8)                            -0.6978098 0.3056066
## genotypeKO:metabolitearginine                         -0.6736076 0.3141654
## genotypeKO:metaboliteCITRIC                            0.0088516 0.3141654
## genotypeKO:metaboliteFUMARIC                          -0.9692841 0.3141654
## genotypeKO:metaboliteglutamine                        -0.7667286 0.3141654
## genotypeKO:metaboliteisoleucine                       -0.7247245 0.3141654
## genotypeKO:metaboliteLACTIC                           -0.8237407 0.3141654
## genotypeKO:metaboliteLC even AC total                 -0.1083540 0.3141654
## genotypeKO:metaboliteLC odd AC total                   0.0839471 0.3141654
## genotypeKO:metaboliteleucine                          -0.6996743 0.3141654
## genotypeKO:metaboliteMALIC                            -0.7985721 0.3141654
## genotypeKO:metaboliteMCAC total                       -0.5409875 0.3141654
## genotypeKO:metaboliteMETHYLSUCCINIC                   -0.9184289 0.3141654
## genotypeKO:metaboliteSUCCINIC-2                       -0.7749538 0.3141654
## genotypeKO:metabolitevaline                           -0.7162648 0.3141654
## chowYellow (C8):metabolitearginine                    -1.2626957 0.2673925
## chowYellow (C8):metaboliteCITRIC                      -0.5439097 0.2949119
## chowYellow (C8):metaboliteFUMARIC                     -0.9714318 0.2891552
## chowYellow (C8):metaboliteglutamine                   -0.9631569 0.2771722
## chowYellow (C8):metaboliteisoleucine                  -1.1705246 0.2876369
## chowYellow (C8):metaboliteLACTIC                      -0.8315275 0.2854145
## chowYellow (C8):metaboliteLC even AC total            -0.7875069 0.2840107
## chowYellow (C8):metaboliteLC odd AC total             -0.8928525 0.2937578
## chowYellow (C8):metaboliteleucine                     -1.0879909 0.2979676
## chowYellow (C8):metaboliteMALIC                       -0.8293624 0.2825788
## chowYellow (C8):metaboliteMCAC total                  -0.9477404 0.2855944
## chowYellow (C8):metaboliteMETHYLSUCCINIC              -0.3178135 0.3020975
## chowYellow (C8):metaboliteSUCCINIC-2                  -0.7365276 0.2795860
## chowYellow (C8):metabolitevaline                      -0.5978008 0.2843059
## genotypeKO:chowYellow (C8):metabolitearginine          1.1183598 0.4199741
## genotypeKO:chowYellow (C8):metaboliteCITRIC            0.5398332 0.4351838
## genotypeKO:chowYellow (C8):metaboliteFUMARIC           0.9596491 0.4323726
## genotypeKO:chowYellow (C8):metaboliteglutamine         0.9859728 0.4269021
## genotypeKO:chowYellow (C8):metaboliteisoleucine        1.0857875 0.4296033
## genotypeKO:chowYellow (C8):metaboliteLACTIC            0.6497277 0.4302308
## genotypeKO:chowYellow (C8):metaboliteLC even AC total  0.7968450 0.4309968
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total   0.1361594 0.4360711
## genotypeKO:chowYellow (C8):metaboliteleucine           1.0353980 0.4380953
## genotypeKO:chowYellow (C8):metaboliteMALIC             0.6371135 0.4288386
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.8934399 0.4326781
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.9976601 0.4415882
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.6302176 0.4266922
## genotypeKO:chowYellow (C8):metabolitevaline            0.7417809 0.4293746
##                                                        DF   t-value
## (Intercept)                                           546  1.612446
## genotypeKO                                             39  2.954568
## chowYellow (C8)                                        39  3.756150
## metabolitearginine                                    546  0.556645
## metaboliteCITRIC                                      546  0.178245
## metaboliteFUMARIC                                     546 -4.531468
## metaboliteglutamine                                   546  8.032804
## metaboliteisoleucine                                  546 -0.124478
## metaboliteLACTIC                                      546  9.685610
## metaboliteLC even AC total                            546  1.830843
## metaboliteLC odd AC total                             546 -2.542079
## metaboliteleucine                                     546  0.431682
## metaboliteMALIC                                       546  3.355421
## metaboliteMCAC total                                  546 -3.230418
## metaboliteMETHYLSUCCINIC                              546 -6.741411
## metaboliteSUCCINIC-2                                  546  4.354784
## metabolitevaline                                      546  0.659164
## genotypeKO:chowYellow (C8)                             39 -2.283360
## genotypeKO:metabolitearginine                         546 -2.144118
## genotypeKO:metaboliteCITRIC                           546  0.028175
## genotypeKO:metaboliteFUMARIC                          546 -3.085267
## genotypeKO:metaboliteglutamine                        546 -2.440526
## genotypeKO:metaboliteisoleucine                       546 -2.306825
## genotypeKO:metaboliteLACTIC                           546 -2.621997
## genotypeKO:metaboliteLC even AC total                 546 -0.344895
## genotypeKO:metaboliteLC odd AC total                  546  0.267207
## genotypeKO:metaboliteleucine                          546 -2.227089
## genotypeKO:metaboliteMALIC                            546 -2.541885
## genotypeKO:metaboliteMCAC total                       546 -1.721983
## genotypeKO:metaboliteMETHYLSUCCINIC                   546 -2.923393
## genotypeKO:metaboliteSUCCINIC-2                       546 -2.466707
## genotypeKO:metabolitevaline                           546 -2.279897
## chowYellow (C8):metabolitearginine                    546 -4.722255
## chowYellow (C8):metaboliteCITRIC                      546 -1.844312
## chowYellow (C8):metaboliteFUMARIC                     546 -3.359551
## chowYellow (C8):metaboliteglutamine                   546 -3.474941
## chowYellow (C8):metaboliteisoleucine                  546 -4.069452
## chowYellow (C8):metaboliteLACTIC                      546 -2.913403
## chowYellow (C8):metaboliteLC even AC total            546 -2.772807
## chowYellow (C8):metaboliteLC odd AC total             546 -3.039417
## chowYellow (C8):metaboliteleucine                     546 -3.651374
## chowYellow (C8):metaboliteMALIC                       546 -2.934977
## chowYellow (C8):metaboliteMCAC total                  546 -3.318483
## chowYellow (C8):metaboliteMETHYLSUCCINIC              546 -1.052023
## chowYellow (C8):metaboliteSUCCINIC-2                  546 -2.634351
## chowYellow (C8):metabolitevaline                      546 -2.102668
## genotypeKO:chowYellow (C8):metabolitearginine         546  2.662926
## genotypeKO:chowYellow (C8):metaboliteCITRIC           546  1.240472
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          546  2.219496
## genotypeKO:chowYellow (C8):metaboliteglutamine        546  2.309600
## genotypeKO:chowYellow (C8):metaboliteisoleucine       546  2.527419
## genotypeKO:chowYellow (C8):metaboliteLACTIC           546  1.510184
## genotypeKO:chowYellow (C8):metaboliteLC even AC total 546  1.848842
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  546  0.312241
## genotypeKO:chowYellow (C8):metaboliteleucine          546  2.363408
## genotypeKO:chowYellow (C8):metaboliteMALIC            546  1.485672
## genotypeKO:chowYellow (C8):metaboliteMCAC total       546  2.064907
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   546  2.259255
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       546  1.476984
## genotypeKO:chowYellow (C8):metabolitevaline           546  1.727584
##                                                       p-value
## (Intercept)                                            0.1074
## genotypeKO                                             0.0053
## chowYellow (C8)                                        0.0006
## metabolitearginine                                     0.5780
## metaboliteCITRIC                                       0.8586
## metaboliteFUMARIC                                      0.0000
## metaboliteglutamine                                    0.0000
## metaboliteisoleucine                                   0.9010
## metaboliteLACTIC                                       0.0000
## metaboliteLC even AC total                             0.0677
## metaboliteLC odd AC total                              0.0113
## metaboliteleucine                                      0.6661
## metaboliteMALIC                                        0.0008
## metaboliteMCAC total                                   0.0013
## metaboliteMETHYLSUCCINIC                               0.0000
## metaboliteSUCCINIC-2                                   0.0000
## metabolitevaline                                       0.5101
## genotypeKO:chowYellow (C8)                             0.0279
## genotypeKO:metabolitearginine                          0.0325
## genotypeKO:metaboliteCITRIC                            0.9775
## genotypeKO:metaboliteFUMARIC                           0.0021
## genotypeKO:metaboliteglutamine                         0.0150
## genotypeKO:metaboliteisoleucine                        0.0214
## genotypeKO:metaboliteLACTIC                            0.0090
## genotypeKO:metaboliteLC even AC total                  0.7303
## genotypeKO:metaboliteLC odd AC total                   0.7894
## genotypeKO:metaboliteleucine                           0.0263
## genotypeKO:metaboliteMALIC                             0.0113
## genotypeKO:metaboliteMCAC total                        0.0856
## genotypeKO:metaboliteMETHYLSUCCINIC                    0.0036
## genotypeKO:metaboliteSUCCINIC-2                        0.0139
## genotypeKO:metabolitevaline                            0.0230
## chowYellow (C8):metabolitearginine                     0.0000
## chowYellow (C8):metaboliteCITRIC                       0.0657
## chowYellow (C8):metaboliteFUMARIC                      0.0008
## chowYellow (C8):metaboliteglutamine                    0.0006
## chowYellow (C8):metaboliteisoleucine                   0.0001
## chowYellow (C8):metaboliteLACTIC                       0.0037
## chowYellow (C8):metaboliteLC even AC total             0.0057
## chowYellow (C8):metaboliteLC odd AC total              0.0025
## chowYellow (C8):metaboliteleucine                      0.0003
## chowYellow (C8):metaboliteMALIC                        0.0035
## chowYellow (C8):metaboliteMCAC total                   0.0010
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.2933
## chowYellow (C8):metaboliteSUCCINIC-2                   0.0087
## chowYellow (C8):metabolitevaline                       0.0360
## genotypeKO:chowYellow (C8):metabolitearginine          0.0080
## genotypeKO:chowYellow (C8):metaboliteCITRIC            0.2153
## genotypeKO:chowYellow (C8):metaboliteFUMARIC           0.0269
## genotypeKO:chowYellow (C8):metaboliteglutamine         0.0213
## genotypeKO:chowYellow (C8):metaboliteisoleucine        0.0118
## genotypeKO:chowYellow (C8):metaboliteLACTIC            0.1316
## genotypeKO:chowYellow (C8):metaboliteLC even AC total  0.0650
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total   0.7550
## genotypeKO:chowYellow (C8):metaboliteleucine           0.0185
## genotypeKO:chowYellow (C8):metaboliteMALIC             0.1379
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.0394
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.0243
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.1403
## genotypeKO:chowYellow (C8):metabolitevaline            0.0846
##  Correlation: 
##                                                       (Intr) gntyKO
## genotypeKO                                            -0.707       
## chowYellow (C8)                                       -0.789  0.558
## metabolitearginine                                    -0.700  0.495
## metaboliteCITRIC                                      -0.700  0.495
## metaboliteFUMARIC                                     -0.700  0.495
## metaboliteglutamine                                   -0.700  0.495
## metaboliteisoleucine                                  -0.700  0.495
## metaboliteLACTIC                                      -0.700  0.495
## metaboliteLC even AC total                            -0.700  0.495
## metaboliteLC odd AC total                             -0.700  0.495
## metaboliteleucine                                     -0.700  0.495
## metaboliteMALIC                                       -0.700  0.495
## metaboliteMCAC total                                  -0.700  0.495
## metaboliteMETHYLSUCCINIC                              -0.700  0.495
## metaboliteSUCCINIC-2                                  -0.700  0.495
## metabolitevaline                                      -0.700  0.495
## genotypeKO:chowYellow (C8)                             0.519 -0.734
## genotypeKO:metabolitearginine                          0.495 -0.700
## genotypeKO:metaboliteCITRIC                            0.495 -0.700
## genotypeKO:metaboliteFUMARIC                           0.495 -0.700
## genotypeKO:metaboliteglutamine                         0.495 -0.700
## genotypeKO:metaboliteisoleucine                        0.495 -0.700
## genotypeKO:metaboliteLACTIC                            0.495 -0.700
## genotypeKO:metaboliteLC even AC total                  0.495 -0.700
## genotypeKO:metaboliteLC odd AC total                   0.495 -0.700
## genotypeKO:metaboliteleucine                           0.495 -0.700
## genotypeKO:metaboliteMALIC                             0.495 -0.700
## genotypeKO:metaboliteMCAC total                        0.495 -0.700
## genotypeKO:metaboliteMETHYLSUCCINIC                    0.495 -0.700
## genotypeKO:metaboliteSUCCINIC-2                        0.495 -0.700
## genotypeKO:metabolitevaline                            0.495 -0.700
## chowYellow (C8):metabolitearginine                     0.582 -0.411
## chowYellow (C8):metaboliteCITRIC                       0.528 -0.373
## chowYellow (C8):metaboliteFUMARIC                      0.538 -0.380
## chowYellow (C8):metaboliteglutamine                    0.561 -0.397
## chowYellow (C8):metaboliteisoleucine                   0.541 -0.382
## chowYellow (C8):metaboliteLACTIC                       0.545 -0.385
## chowYellow (C8):metaboliteLC even AC total             0.548 -0.387
## chowYellow (C8):metaboliteLC odd AC total              0.530 -0.375
## chowYellow (C8):metaboliteleucine                      0.522 -0.369
## chowYellow (C8):metaboliteMALIC                        0.551 -0.389
## chowYellow (C8):metaboliteMCAC total                   0.545 -0.385
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.515 -0.364
## chowYellow (C8):metaboliteSUCCINIC-2                   0.557 -0.394
## chowYellow (C8):metabolitevaline                       0.547 -0.387
## genotypeKO:chowYellow (C8):metabolitearginine         -0.370  0.524
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.358  0.506
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.360  0.509
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.364  0.515
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.362  0.512
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.362  0.511
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.361  0.511
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.357  0.505
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.355  0.502
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.363  0.513
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.360  0.509
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.352  0.498
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.365  0.516
## genotypeKO:chowYellow (C8):metabolitevaline           -0.362  0.512
##                                                       chY(C8) mtbltr
## genotypeKO                                                          
## chowYellow (C8)                                                     
## metabolitearginine                                     0.553        
## metaboliteCITRIC                                       0.553   0.500
## metaboliteFUMARIC                                      0.553   0.500
## metaboliteglutamine                                    0.553   0.500
## metaboliteisoleucine                                   0.553   0.500
## metaboliteLACTIC                                       0.553   0.500
## metaboliteLC even AC total                             0.553   0.500
## metaboliteLC odd AC total                              0.553   0.500
## metaboliteleucine                                      0.553   0.500
## metaboliteMALIC                                        0.553   0.500
## metaboliteMCAC total                                   0.553   0.500
## metaboliteMETHYLSUCCINIC                               0.553   0.500
## metaboliteSUCCINIC-2                                   0.553   0.500
## metabolitevaline                                       0.553   0.500
## genotypeKO:chowYellow (C8)                            -0.659  -0.363
## genotypeKO:metabolitearginine                         -0.391  -0.707
## genotypeKO:metaboliteCITRIC                           -0.391  -0.354
## genotypeKO:metaboliteFUMARIC                          -0.391  -0.354
## genotypeKO:metaboliteglutamine                        -0.391  -0.354
## genotypeKO:metaboliteisoleucine                       -0.391  -0.354
## genotypeKO:metaboliteLACTIC                           -0.391  -0.354
## genotypeKO:metaboliteLC even AC total                 -0.391  -0.354
## genotypeKO:metaboliteLC odd AC total                  -0.391  -0.354
## genotypeKO:metaboliteleucine                          -0.391  -0.354
## genotypeKO:metaboliteMALIC                            -0.391  -0.354
## genotypeKO:metaboliteMCAC total                       -0.391  -0.354
## genotypeKO:metaboliteMETHYLSUCCINIC                   -0.391  -0.354
## genotypeKO:metaboliteSUCCINIC-2                       -0.391  -0.354
## genotypeKO:metabolitevaline                           -0.391  -0.354
## chowYellow (C8):metabolitearginine                    -0.711  -0.831
## chowYellow (C8):metaboliteCITRIC                      -0.710  -0.377
## chowYellow (C8):metaboliteFUMARIC                     -0.668  -0.384
## chowYellow (C8):metaboliteglutamine                   -0.643  -0.401
## chowYellow (C8):metaboliteisoleucine                  -0.724  -0.386
## chowYellow (C8):metaboliteLACTIC                      -0.684  -0.389
## chowYellow (C8):metaboliteLC even AC total            -0.657  -0.391
## chowYellow (C8):metaboliteLC odd AC total             -0.714  -0.378
## chowYellow (C8):metaboliteleucine                     -0.683  -0.373
## chowYellow (C8):metaboliteMALIC                       -0.668  -0.393
## chowYellow (C8):metaboliteMCAC total                  -0.695  -0.389
## chowYellow (C8):metaboliteMETHYLSUCCINIC              -0.691  -0.368
## chowYellow (C8):metaboliteSUCCINIC-2                  -0.654  -0.397
## chowYellow (C8):metabolitevaline                      -0.632  -0.391
## genotypeKO:chowYellow (C8):metabolitearginine          0.455   0.529
## genotypeKO:chowYellow (C8):metaboliteCITRIC            0.483   0.255
## genotypeKO:chowYellow (C8):metaboliteFUMARIC           0.451   0.257
## genotypeKO:chowYellow (C8):metaboliteglutamine         0.418   0.260
## genotypeKO:chowYellow (C8):metaboliteisoleucine        0.482   0.259
## genotypeKO:chowYellow (C8):metaboliteLACTIC            0.456   0.258
## genotypeKO:chowYellow (C8):metaboliteLC even AC total  0.436   0.258
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total   0.488   0.255
## genotypeKO:chowYellow (C8):metaboliteleucine           0.463   0.254
## genotypeKO:chowYellow (C8):metaboliteMALIC             0.442   0.259
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.463   0.257
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.475   0.252
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.431   0.260
## genotypeKO:chowYellow (C8):metabolitevaline            0.420   0.259
##                                                       mCITRI mFUMAR mtbltg
## genotypeKO                                                                
## chowYellow (C8)                                                           
## metabolitearginine                                                        
## metaboliteCITRIC                                                          
## metaboliteFUMARIC                                      0.500              
## metaboliteglutamine                                    0.500  0.500       
## metaboliteisoleucine                                   0.500  0.500  0.500
## metaboliteLACTIC                                       0.500  0.500  0.500
## metaboliteLC even AC total                             0.500  0.500  0.500
## metaboliteLC odd AC total                              0.500  0.500  0.500
## metaboliteleucine                                      0.500  0.500  0.500
## metaboliteMALIC                                        0.500  0.500  0.500
## metaboliteMCAC total                                   0.500  0.500  0.500
## metaboliteMETHYLSUCCINIC                               0.500  0.500  0.500
## metaboliteSUCCINIC-2                                   0.500  0.500  0.500
## metabolitevaline                                       0.500  0.500  0.500
## genotypeKO:chowYellow (C8)                            -0.363 -0.363 -0.363
## genotypeKO:metabolitearginine                         -0.354 -0.354 -0.354
## genotypeKO:metaboliteCITRIC                           -0.707 -0.354 -0.354
## genotypeKO:metaboliteFUMARIC                          -0.354 -0.707 -0.354
## genotypeKO:metaboliteglutamine                        -0.354 -0.354 -0.707
## genotypeKO:metaboliteisoleucine                       -0.354 -0.354 -0.354
## genotypeKO:metaboliteLACTIC                           -0.354 -0.354 -0.354
## genotypeKO:metaboliteLC even AC total                 -0.354 -0.354 -0.354
## genotypeKO:metaboliteLC odd AC total                  -0.354 -0.354 -0.354
## genotypeKO:metaboliteleucine                          -0.354 -0.354 -0.354
## genotypeKO:metaboliteMALIC                            -0.354 -0.354 -0.354
## genotypeKO:metaboliteMCAC total                       -0.354 -0.354 -0.354
## genotypeKO:metaboliteMETHYLSUCCINIC                   -0.354 -0.354 -0.354
## genotypeKO:metaboliteSUCCINIC-2                       -0.354 -0.354 -0.354
## genotypeKO:metabolitevaline                           -0.354 -0.354 -0.354
## chowYellow (C8):metabolitearginine                    -0.415 -0.415 -0.415
## chowYellow (C8):metaboliteCITRIC                      -0.753 -0.377 -0.377
## chowYellow (C8):metaboliteFUMARIC                     -0.384 -0.768 -0.384
## chowYellow (C8):metaboliteglutamine                   -0.401 -0.401 -0.801
## chowYellow (C8):metaboliteisoleucine                  -0.386 -0.386 -0.386
## chowYellow (C8):metaboliteLACTIC                      -0.389 -0.389 -0.389
## chowYellow (C8):metaboliteLC even AC total            -0.391 -0.391 -0.391
## chowYellow (C8):metaboliteLC odd AC total             -0.378 -0.378 -0.378
## chowYellow (C8):metaboliteleucine                     -0.373 -0.373 -0.373
## chowYellow (C8):metaboliteMALIC                       -0.393 -0.393 -0.393
## chowYellow (C8):metaboliteMCAC total                  -0.389 -0.389 -0.389
## chowYellow (C8):metaboliteMETHYLSUCCINIC              -0.368 -0.368 -0.368
## chowYellow (C8):metaboliteSUCCINIC-2                  -0.397 -0.397 -0.397
## chowYellow (C8):metabolitevaline                      -0.391 -0.391 -0.391
## genotypeKO:chowYellow (C8):metabolitearginine          0.264  0.264  0.264
## genotypeKO:chowYellow (C8):metaboliteCITRIC            0.510  0.255  0.255
## genotypeKO:chowYellow (C8):metaboliteFUMARIC           0.257  0.514  0.257
## genotypeKO:chowYellow (C8):metaboliteglutamine         0.260  0.260  0.520
## genotypeKO:chowYellow (C8):metaboliteisoleucine        0.259  0.259  0.259
## genotypeKO:chowYellow (C8):metaboliteLACTIC            0.258  0.258  0.258
## genotypeKO:chowYellow (C8):metaboliteLC even AC total  0.258  0.258  0.258
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total   0.255  0.255  0.255
## genotypeKO:chowYellow (C8):metaboliteleucine           0.254  0.254  0.254
## genotypeKO:chowYellow (C8):metaboliteMALIC             0.259  0.259  0.259
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.257  0.257  0.257
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.252  0.252  0.252
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.260  0.260  0.260
## genotypeKO:chowYellow (C8):metabolitevaline            0.259  0.259  0.259
##                                                       mtblts mLACTI mLCeAt
## genotypeKO                                                                
## chowYellow (C8)                                                           
## metabolitearginine                                                        
## metaboliteCITRIC                                                          
## metaboliteFUMARIC                                                         
## metaboliteglutamine                                                       
## metaboliteisoleucine                                                      
## metaboliteLACTIC                                       0.500              
## metaboliteLC even AC total                             0.500  0.500       
## metaboliteLC odd AC total                              0.500  0.500  0.500
## metaboliteleucine                                      0.500  0.500  0.500
## metaboliteMALIC                                        0.500  0.500  0.500
## metaboliteMCAC total                                   0.500  0.500  0.500
## metaboliteMETHYLSUCCINIC                               0.500  0.500  0.500
## metaboliteSUCCINIC-2                                   0.500  0.500  0.500
## metabolitevaline                                       0.500  0.500  0.500
## genotypeKO:chowYellow (C8)                            -0.363 -0.363 -0.363
## genotypeKO:metabolitearginine                         -0.354 -0.354 -0.354
## genotypeKO:metaboliteCITRIC                           -0.354 -0.354 -0.354
## genotypeKO:metaboliteFUMARIC                          -0.354 -0.354 -0.354
## genotypeKO:metaboliteglutamine                        -0.354 -0.354 -0.354
## genotypeKO:metaboliteisoleucine                       -0.707 -0.354 -0.354
## genotypeKO:metaboliteLACTIC                           -0.354 -0.707 -0.354
## genotypeKO:metaboliteLC even AC total                 -0.354 -0.354 -0.707
## genotypeKO:metaboliteLC odd AC total                  -0.354 -0.354 -0.354
## genotypeKO:metaboliteleucine                          -0.354 -0.354 -0.354
## genotypeKO:metaboliteMALIC                            -0.354 -0.354 -0.354
## genotypeKO:metaboliteMCAC total                       -0.354 -0.354 -0.354
## genotypeKO:metaboliteMETHYLSUCCINIC                   -0.354 -0.354 -0.354
## genotypeKO:metaboliteSUCCINIC-2                       -0.354 -0.354 -0.354
## genotypeKO:metabolitevaline                           -0.354 -0.354 -0.354
## chowYellow (C8):metabolitearginine                    -0.415 -0.415 -0.415
## chowYellow (C8):metaboliteCITRIC                      -0.377 -0.377 -0.377
## chowYellow (C8):metaboliteFUMARIC                     -0.384 -0.384 -0.384
## chowYellow (C8):metaboliteglutamine                   -0.401 -0.401 -0.401
## chowYellow (C8):metaboliteisoleucine                  -0.772 -0.386 -0.386
## chowYellow (C8):metaboliteLACTIC                      -0.389 -0.778 -0.389
## chowYellow (C8):metaboliteLC even AC total            -0.391 -0.391 -0.782
## chowYellow (C8):metaboliteLC odd AC total             -0.378 -0.378 -0.378
## chowYellow (C8):metaboliteleucine                     -0.373 -0.373 -0.373
## chowYellow (C8):metaboliteMALIC                       -0.393 -0.393 -0.393
## chowYellow (C8):metaboliteMCAC total                  -0.389 -0.389 -0.389
## chowYellow (C8):metaboliteMETHYLSUCCINIC              -0.368 -0.368 -0.368
## chowYellow (C8):metaboliteSUCCINIC-2                  -0.397 -0.397 -0.397
## chowYellow (C8):metabolitevaline                      -0.391 -0.391 -0.391
## genotypeKO:chowYellow (C8):metabolitearginine          0.264  0.264  0.264
## genotypeKO:chowYellow (C8):metaboliteCITRIC            0.255  0.255  0.255
## genotypeKO:chowYellow (C8):metaboliteFUMARIC           0.257  0.257  0.257
## genotypeKO:chowYellow (C8):metaboliteglutamine         0.260  0.260  0.260
## genotypeKO:chowYellow (C8):metaboliteisoleucine        0.517  0.259  0.259
## genotypeKO:chowYellow (C8):metaboliteLACTIC            0.258  0.516  0.258
## genotypeKO:chowYellow (C8):metaboliteLC even AC total  0.258  0.258  0.515
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total   0.255  0.255  0.255
## genotypeKO:chowYellow (C8):metaboliteleucine           0.254  0.254  0.254
## genotypeKO:chowYellow (C8):metaboliteMALIC             0.259  0.259  0.259
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.257  0.257  0.257
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.252  0.252  0.252
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.260  0.260  0.260
## genotypeKO:chowYellow (C8):metabolitevaline            0.259  0.259  0.259
##                                                       mLCoAt mtbltl mMALIC
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
## metaboliteleucine                                      0.500              
## metaboliteMALIC                                        0.500  0.500       
## metaboliteMCAC total                                   0.500  0.500  0.500
## metaboliteMETHYLSUCCINIC                               0.500  0.500  0.500
## metaboliteSUCCINIC-2                                   0.500  0.500  0.500
## metabolitevaline                                       0.500  0.500  0.500
## genotypeKO:chowYellow (C8)                            -0.363 -0.363 -0.363
## genotypeKO:metabolitearginine                         -0.354 -0.354 -0.354
## genotypeKO:metaboliteCITRIC                           -0.354 -0.354 -0.354
## genotypeKO:metaboliteFUMARIC                          -0.354 -0.354 -0.354
## genotypeKO:metaboliteglutamine                        -0.354 -0.354 -0.354
## genotypeKO:metaboliteisoleucine                       -0.354 -0.354 -0.354
## genotypeKO:metaboliteLACTIC                           -0.354 -0.354 -0.354
## genotypeKO:metaboliteLC even AC total                 -0.354 -0.354 -0.354
## genotypeKO:metaboliteLC odd AC total                  -0.707 -0.354 -0.354
## genotypeKO:metaboliteleucine                          -0.354 -0.707 -0.354
## genotypeKO:metaboliteMALIC                            -0.354 -0.354 -0.707
## genotypeKO:metaboliteMCAC total                       -0.354 -0.354 -0.354
## genotypeKO:metaboliteMETHYLSUCCINIC                   -0.354 -0.354 -0.354
## genotypeKO:metaboliteSUCCINIC-2                       -0.354 -0.354 -0.354
## genotypeKO:metabolitevaline                           -0.354 -0.354 -0.354
## chowYellow (C8):metabolitearginine                    -0.415 -0.415 -0.415
## chowYellow (C8):metaboliteCITRIC                      -0.377 -0.377 -0.377
## chowYellow (C8):metaboliteFUMARIC                     -0.384 -0.384 -0.384
## chowYellow (C8):metaboliteglutamine                   -0.401 -0.401 -0.401
## chowYellow (C8):metaboliteisoleucine                  -0.386 -0.386 -0.386
## chowYellow (C8):metaboliteLACTIC                      -0.389 -0.389 -0.389
## chowYellow (C8):metaboliteLC even AC total            -0.391 -0.391 -0.391
## chowYellow (C8):metaboliteLC odd AC total             -0.756 -0.378 -0.378
## chowYellow (C8):metaboliteleucine                     -0.373 -0.746 -0.373
## chowYellow (C8):metaboliteMALIC                       -0.393 -0.393 -0.786
## chowYellow (C8):metaboliteMCAC total                  -0.389 -0.389 -0.389
## chowYellow (C8):metaboliteMETHYLSUCCINIC              -0.368 -0.368 -0.368
## chowYellow (C8):metaboliteSUCCINIC-2                  -0.397 -0.397 -0.397
## chowYellow (C8):metabolitevaline                      -0.391 -0.391 -0.391
## genotypeKO:chowYellow (C8):metabolitearginine          0.264  0.264  0.264
## genotypeKO:chowYellow (C8):metaboliteCITRIC            0.255  0.255  0.255
## genotypeKO:chowYellow (C8):metaboliteFUMARIC           0.257  0.257  0.257
## genotypeKO:chowYellow (C8):metaboliteglutamine         0.260  0.260  0.260
## genotypeKO:chowYellow (C8):metaboliteisoleucine        0.259  0.259  0.259
## genotypeKO:chowYellow (C8):metaboliteLACTIC            0.258  0.258  0.258
## genotypeKO:chowYellow (C8):metaboliteLC even AC total  0.258  0.258  0.258
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total   0.509  0.255  0.255
## genotypeKO:chowYellow (C8):metaboliteleucine           0.254  0.507  0.254
## genotypeKO:chowYellow (C8):metaboliteMALIC             0.259  0.259  0.518
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.257  0.257  0.257
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.252  0.252  0.252
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.260  0.260  0.260
## genotypeKO:chowYellow (C8):metabolitevaline            0.259  0.259  0.259
##                                                       mMCACt mMETHY mSUCCI
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
## metaboliteMETHYLSUCCINIC                               0.500              
## metaboliteSUCCINIC-2                                   0.500  0.500       
## metabolitevaline                                       0.500  0.500  0.500
## genotypeKO:chowYellow (C8)                            -0.363 -0.363 -0.363
## genotypeKO:metabolitearginine                         -0.354 -0.354 -0.354
## genotypeKO:metaboliteCITRIC                           -0.354 -0.354 -0.354
## genotypeKO:metaboliteFUMARIC                          -0.354 -0.354 -0.354
## genotypeKO:metaboliteglutamine                        -0.354 -0.354 -0.354
## genotypeKO:metaboliteisoleucine                       -0.354 -0.354 -0.354
## genotypeKO:metaboliteLACTIC                           -0.354 -0.354 -0.354
## genotypeKO:metaboliteLC even AC total                 -0.354 -0.354 -0.354
## genotypeKO:metaboliteLC odd AC total                  -0.354 -0.354 -0.354
## genotypeKO:metaboliteleucine                          -0.354 -0.354 -0.354
## genotypeKO:metaboliteMALIC                            -0.354 -0.354 -0.354
## genotypeKO:metaboliteMCAC total                       -0.707 -0.354 -0.354
## genotypeKO:metaboliteMETHYLSUCCINIC                   -0.354 -0.707 -0.354
## genotypeKO:metaboliteSUCCINIC-2                       -0.354 -0.354 -0.707
## genotypeKO:metabolitevaline                           -0.354 -0.354 -0.354
## chowYellow (C8):metabolitearginine                    -0.415 -0.415 -0.415
## chowYellow (C8):metaboliteCITRIC                      -0.377 -0.377 -0.377
## chowYellow (C8):metaboliteFUMARIC                     -0.384 -0.384 -0.384
## chowYellow (C8):metaboliteglutamine                   -0.401 -0.401 -0.401
## chowYellow (C8):metaboliteisoleucine                  -0.386 -0.386 -0.386
## chowYellow (C8):metaboliteLACTIC                      -0.389 -0.389 -0.389
## chowYellow (C8):metaboliteLC even AC total            -0.391 -0.391 -0.391
## chowYellow (C8):metaboliteLC odd AC total             -0.378 -0.378 -0.378
## chowYellow (C8):metaboliteleucine                     -0.373 -0.373 -0.373
## chowYellow (C8):metaboliteMALIC                       -0.393 -0.393 -0.393
## chowYellow (C8):metaboliteMCAC total                  -0.778 -0.389 -0.389
## chowYellow (C8):metaboliteMETHYLSUCCINIC              -0.368 -0.735 -0.368
## chowYellow (C8):metaboliteSUCCINIC-2                  -0.397 -0.397 -0.795
## chowYellow (C8):metabolitevaline                      -0.391 -0.391 -0.391
## genotypeKO:chowYellow (C8):metabolitearginine          0.264  0.264  0.264
## genotypeKO:chowYellow (C8):metaboliteCITRIC            0.255  0.255  0.255
## genotypeKO:chowYellow (C8):metaboliteFUMARIC           0.257  0.257  0.257
## genotypeKO:chowYellow (C8):metaboliteglutamine         0.260  0.260  0.260
## genotypeKO:chowYellow (C8):metaboliteisoleucine        0.259  0.259  0.259
## genotypeKO:chowYellow (C8):metaboliteLACTIC            0.258  0.258  0.258
## genotypeKO:chowYellow (C8):metaboliteLC even AC total  0.258  0.258  0.258
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total   0.255  0.255  0.255
## genotypeKO:chowYellow (C8):metaboliteleucine           0.254  0.254  0.254
## genotypeKO:chowYellow (C8):metaboliteMALIC             0.259  0.259  0.259
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.513  0.257  0.257
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.252  0.503  0.252
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.260  0.260  0.521
## genotypeKO:chowYellow (C8):metabolitevaline            0.259  0.259  0.259
##                                                       mtbltv gnKO:Y(C8)
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
## genotypeKO:chowYellow (C8)                            -0.363           
## genotypeKO:metabolitearginine                         -0.354  0.514    
## genotypeKO:metaboliteCITRIC                           -0.354  0.514    
## genotypeKO:metaboliteFUMARIC                          -0.354  0.514    
## genotypeKO:metaboliteglutamine                        -0.354  0.514    
## genotypeKO:metaboliteisoleucine                       -0.354  0.514    
## genotypeKO:metaboliteLACTIC                           -0.354  0.514    
## genotypeKO:metaboliteLC even AC total                 -0.354  0.514    
## genotypeKO:metaboliteLC odd AC total                  -0.354  0.514    
## genotypeKO:metaboliteleucine                          -0.354  0.514    
## genotypeKO:metaboliteMALIC                            -0.354  0.514    
## genotypeKO:metaboliteMCAC total                       -0.354  0.514    
## genotypeKO:metaboliteMETHYLSUCCINIC                   -0.354  0.514    
## genotypeKO:metaboliteSUCCINIC-2                       -0.354  0.514    
## genotypeKO:metabolitevaline                           -0.707  0.514    
## chowYellow (C8):metabolitearginine                    -0.415  0.473    
## chowYellow (C8):metaboliteCITRIC                      -0.377  0.472    
## chowYellow (C8):metaboliteFUMARIC                     -0.384  0.440    
## chowYellow (C8):metaboliteglutamine                   -0.401  0.425    
## chowYellow (C8):metaboliteisoleucine                  -0.386  0.474    
## chowYellow (C8):metaboliteLACTIC                      -0.389  0.451    
## chowYellow (C8):metaboliteLC even AC total            -0.391  0.434    
## chowYellow (C8):metaboliteLC odd AC total             -0.378  0.470    
## chowYellow (C8):metaboliteleucine                     -0.373  0.451    
## chowYellow (C8):metaboliteMALIC                       -0.393  0.442    
## chowYellow (C8):metaboliteMCAC total                  -0.389  0.461    
## chowYellow (C8):metaboliteMETHYLSUCCINIC              -0.368  0.455    
## chowYellow (C8):metaboliteSUCCINIC-2                  -0.397  0.432    
## chowYellow (C8):metabolitevaline                      -0.781  0.418    
## genotypeKO:chowYellow (C8):metabolitearginine          0.264 -0.707    
## genotypeKO:chowYellow (C8):metaboliteCITRIC            0.255 -0.709    
## genotypeKO:chowYellow (C8):metaboliteFUMARIC           0.257 -0.689    
## genotypeKO:chowYellow (C8):metaboliteglutamine         0.260 -0.675    
## genotypeKO:chowYellow (C8):metaboliteisoleucine        0.259 -0.711    
## genotypeKO:chowYellow (C8):metaboliteLACTIC            0.258 -0.695    
## genotypeKO:chowYellow (C8):metaboliteLC even AC total  0.258 -0.681    
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total   0.255 -0.710    
## genotypeKO:chowYellow (C8):metaboliteleucine           0.254 -0.691    
## genotypeKO:chowYellow (C8):metaboliteMALIC             0.259 -0.688    
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.257 -0.694    
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.252 -0.697    
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.260 -0.682    
## genotypeKO:chowYellow (C8):metabolitevaline            0.517 -0.671    
##                                                       gntypKO:mtbltr
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
## genotypeKO:metaboliteCITRIC                            0.500        
## genotypeKO:metaboliteFUMARIC                           0.500        
## genotypeKO:metaboliteglutamine                         0.500        
## genotypeKO:metaboliteisoleucine                        0.500        
## genotypeKO:metaboliteLACTIC                            0.500        
## genotypeKO:metaboliteLC even AC total                  0.500        
## genotypeKO:metaboliteLC odd AC total                   0.500        
## genotypeKO:metaboliteleucine                           0.500        
## genotypeKO:metaboliteMALIC                             0.500        
## genotypeKO:metaboliteMCAC total                        0.500        
## genotypeKO:metaboliteMETHYLSUCCINIC                    0.500        
## genotypeKO:metaboliteSUCCINIC-2                        0.500        
## genotypeKO:metabolitevaline                            0.500        
## chowYellow (C8):metabolitearginine                     0.587        
## chowYellow (C8):metaboliteCITRIC                       0.266        
## chowYellow (C8):metaboliteFUMARIC                      0.272        
## chowYellow (C8):metaboliteglutamine                    0.283        
## chowYellow (C8):metaboliteisoleucine                   0.273        
## chowYellow (C8):metaboliteLACTIC                       0.275        
## chowYellow (C8):metaboliteLC even AC total             0.277        
## chowYellow (C8):metaboliteLC odd AC total              0.267        
## chowYellow (C8):metaboliteleucine                      0.264        
## chowYellow (C8):metaboliteMALIC                        0.278        
## chowYellow (C8):metaboliteMCAC total                   0.275        
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.260        
## chowYellow (C8):metaboliteSUCCINIC-2                   0.281        
## chowYellow (C8):metabolitevaline                       0.276        
## genotypeKO:chowYellow (C8):metabolitearginine         -0.748        
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.361        
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.363        
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.368        
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.366        
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.365        
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.364        
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.360        
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.359        
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.366        
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.363        
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.356        
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.368        
## genotypeKO:chowYellow (C8):metabolitevaline           -0.366        
##                                                       gKO:CI gKO:FU
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
## genotypeKO:metaboliteFUMARIC                           0.500       
## genotypeKO:metaboliteglutamine                         0.500  0.500
## genotypeKO:metaboliteisoleucine                        0.500  0.500
## genotypeKO:metaboliteLACTIC                            0.500  0.500
## genotypeKO:metaboliteLC even AC total                  0.500  0.500
## genotypeKO:metaboliteLC odd AC total                   0.500  0.500
## genotypeKO:metaboliteleucine                           0.500  0.500
## genotypeKO:metaboliteMALIC                             0.500  0.500
## genotypeKO:metaboliteMCAC total                        0.500  0.500
## genotypeKO:metaboliteMETHYLSUCCINIC                    0.500  0.500
## genotypeKO:metaboliteSUCCINIC-2                        0.500  0.500
## genotypeKO:metabolitevaline                            0.500  0.500
## chowYellow (C8):metabolitearginine                     0.294  0.294
## chowYellow (C8):metaboliteCITRIC                       0.533  0.266
## chowYellow (C8):metaboliteFUMARIC                      0.272  0.543
## chowYellow (C8):metaboliteglutamine                    0.283  0.283
## chowYellow (C8):metaboliteisoleucine                   0.273  0.273
## chowYellow (C8):metaboliteLACTIC                       0.275  0.275
## chowYellow (C8):metaboliteLC even AC total             0.277  0.277
## chowYellow (C8):metaboliteLC odd AC total              0.267  0.267
## chowYellow (C8):metaboliteleucine                      0.264  0.264
## chowYellow (C8):metaboliteMALIC                        0.278  0.278
## chowYellow (C8):metaboliteMCAC total                   0.275  0.275
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.260  0.260
## chowYellow (C8):metaboliteSUCCINIC-2                   0.281  0.281
## chowYellow (C8):metabolitevaline                       0.276  0.276
## genotypeKO:chowYellow (C8):metabolitearginine         -0.374 -0.374
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.722 -0.361
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.363 -0.727
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.368 -0.368
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.366 -0.366
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.365 -0.365
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.364 -0.364
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.360 -0.360
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.359 -0.359
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.366 -0.366
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.363 -0.363
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.356 -0.356
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.368 -0.368
## genotypeKO:chowYellow (C8):metabolitevaline           -0.366 -0.366
##                                                       gntypKO:mtbltg
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
## genotypeKO:metaboliteisoleucine                        0.500        
## genotypeKO:metaboliteLACTIC                            0.500        
## genotypeKO:metaboliteLC even AC total                  0.500        
## genotypeKO:metaboliteLC odd AC total                   0.500        
## genotypeKO:metaboliteleucine                           0.500        
## genotypeKO:metaboliteMALIC                             0.500        
## genotypeKO:metaboliteMCAC total                        0.500        
## genotypeKO:metaboliteMETHYLSUCCINIC                    0.500        
## genotypeKO:metaboliteSUCCINIC-2                        0.500        
## genotypeKO:metabolitevaline                            0.500        
## chowYellow (C8):metabolitearginine                     0.294        
## chowYellow (C8):metaboliteCITRIC                       0.266        
## chowYellow (C8):metaboliteFUMARIC                      0.272        
## chowYellow (C8):metaboliteglutamine                    0.567        
## chowYellow (C8):metaboliteisoleucine                   0.273        
## chowYellow (C8):metaboliteLACTIC                       0.275        
## chowYellow (C8):metaboliteLC even AC total             0.277        
## chowYellow (C8):metaboliteLC odd AC total              0.267        
## chowYellow (C8):metaboliteleucine                      0.264        
## chowYellow (C8):metaboliteMALIC                        0.278        
## chowYellow (C8):metaboliteMCAC total                   0.275        
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.260        
## chowYellow (C8):metaboliteSUCCINIC-2                   0.281        
## chowYellow (C8):metabolitevaline                       0.276        
## genotypeKO:chowYellow (C8):metabolitearginine         -0.374        
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.361        
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.363        
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.736        
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.366        
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.365        
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.364        
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.360        
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.359        
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.366        
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.363        
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.356        
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.368        
## genotypeKO:chowYellow (C8):metabolitevaline           -0.366        
##                                                       gntypKO:mtblts
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
## genotypeKO:metaboliteLACTIC                            0.500        
## genotypeKO:metaboliteLC even AC total                  0.500        
## genotypeKO:metaboliteLC odd AC total                   0.500        
## genotypeKO:metaboliteleucine                           0.500        
## genotypeKO:metaboliteMALIC                             0.500        
## genotypeKO:metaboliteMCAC total                        0.500        
## genotypeKO:metaboliteMETHYLSUCCINIC                    0.500        
## genotypeKO:metaboliteSUCCINIC-2                        0.500        
## genotypeKO:metabolitevaline                            0.500        
## chowYellow (C8):metabolitearginine                     0.294        
## chowYellow (C8):metaboliteCITRIC                       0.266        
## chowYellow (C8):metaboliteFUMARIC                      0.272        
## chowYellow (C8):metaboliteglutamine                    0.283        
## chowYellow (C8):metaboliteisoleucine                   0.546        
## chowYellow (C8):metaboliteLACTIC                       0.275        
## chowYellow (C8):metaboliteLC even AC total             0.277        
## chowYellow (C8):metaboliteLC odd AC total              0.267        
## chowYellow (C8):metaboliteleucine                      0.264        
## chowYellow (C8):metaboliteMALIC                        0.278        
## chowYellow (C8):metaboliteMCAC total                   0.275        
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.260        
## chowYellow (C8):metaboliteSUCCINIC-2                   0.281        
## chowYellow (C8):metabolitevaline                       0.276        
## genotypeKO:chowYellow (C8):metabolitearginine         -0.374        
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.361        
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.363        
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.368        
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.731        
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.365        
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.364        
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.360        
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.359        
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.366        
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.363        
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.356        
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.368        
## genotypeKO:chowYellow (C8):metabolitevaline           -0.366        
##                                                       gKO:LA gKOeAt gKOoAt
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
## genotypeKO:metaboliteLC even AC total                  0.500              
## genotypeKO:metaboliteLC odd AC total                   0.500  0.500       
## genotypeKO:metaboliteleucine                           0.500  0.500  0.500
## genotypeKO:metaboliteMALIC                             0.500  0.500  0.500
## genotypeKO:metaboliteMCAC total                        0.500  0.500  0.500
## genotypeKO:metaboliteMETHYLSUCCINIC                    0.500  0.500  0.500
## genotypeKO:metaboliteSUCCINIC-2                        0.500  0.500  0.500
## genotypeKO:metabolitevaline                            0.500  0.500  0.500
## chowYellow (C8):metabolitearginine                     0.294  0.294  0.294
## chowYellow (C8):metaboliteCITRIC                       0.266  0.266  0.266
## chowYellow (C8):metaboliteFUMARIC                      0.272  0.272  0.272
## chowYellow (C8):metaboliteglutamine                    0.283  0.283  0.283
## chowYellow (C8):metaboliteisoleucine                   0.273  0.273  0.273
## chowYellow (C8):metaboliteLACTIC                       0.550  0.275  0.275
## chowYellow (C8):metaboliteLC even AC total             0.277  0.553  0.277
## chowYellow (C8):metaboliteLC odd AC total              0.267  0.267  0.535
## chowYellow (C8):metaboliteleucine                      0.264  0.264  0.264
## chowYellow (C8):metaboliteMALIC                        0.278  0.278  0.278
## chowYellow (C8):metaboliteMCAC total                   0.275  0.275  0.275
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.260  0.260  0.260
## chowYellow (C8):metaboliteSUCCINIC-2                   0.281  0.281  0.281
## chowYellow (C8):metabolitevaline                       0.276  0.276  0.276
## genotypeKO:chowYellow (C8):metabolitearginine         -0.374 -0.374 -0.374
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.361 -0.361 -0.361
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.363 -0.363 -0.363
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.368 -0.368 -0.368
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.366 -0.366 -0.366
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.730 -0.365 -0.365
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.364 -0.729 -0.364
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.360 -0.360 -0.720
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.359 -0.359 -0.359
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.366 -0.366 -0.366
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.363 -0.363 -0.363
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.356 -0.356 -0.356
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.368 -0.368 -0.368
## genotypeKO:chowYellow (C8):metabolitevaline           -0.366 -0.366 -0.366
##                                                       gntypKO:mtbltl
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
## genotypeKO:metaboliteMALIC                             0.500        
## genotypeKO:metaboliteMCAC total                        0.500        
## genotypeKO:metaboliteMETHYLSUCCINIC                    0.500        
## genotypeKO:metaboliteSUCCINIC-2                        0.500        
## genotypeKO:metabolitevaline                            0.500        
## chowYellow (C8):metabolitearginine                     0.294        
## chowYellow (C8):metaboliteCITRIC                       0.266        
## chowYellow (C8):metaboliteFUMARIC                      0.272        
## chowYellow (C8):metaboliteglutamine                    0.283        
## chowYellow (C8):metaboliteisoleucine                   0.273        
## chowYellow (C8):metaboliteLACTIC                       0.275        
## chowYellow (C8):metaboliteLC even AC total             0.277        
## chowYellow (C8):metaboliteLC odd AC total              0.267        
## chowYellow (C8):metaboliteleucine                      0.527        
## chowYellow (C8):metaboliteMALIC                        0.278        
## chowYellow (C8):metaboliteMCAC total                   0.275        
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.260        
## chowYellow (C8):metaboliteSUCCINIC-2                   0.281        
## chowYellow (C8):metabolitevaline                       0.276        
## genotypeKO:chowYellow (C8):metabolitearginine         -0.374        
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.361        
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.363        
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.368        
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.366        
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.365        
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.364        
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.360        
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.717        
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.366        
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.363        
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.356        
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.368        
## genotypeKO:chowYellow (C8):metabolitevaline           -0.366        
##                                                       gKO:MA gKO:Mt gKO:ME
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
## genotypeKO:metaboliteMCAC total                        0.500              
## genotypeKO:metaboliteMETHYLSUCCINIC                    0.500  0.500       
## genotypeKO:metaboliteSUCCINIC-2                        0.500  0.500  0.500
## genotypeKO:metabolitevaline                            0.500  0.500  0.500
## chowYellow (C8):metabolitearginine                     0.294  0.294  0.294
## chowYellow (C8):metaboliteCITRIC                       0.266  0.266  0.266
## chowYellow (C8):metaboliteFUMARIC                      0.272  0.272  0.272
## chowYellow (C8):metaboliteglutamine                    0.283  0.283  0.283
## chowYellow (C8):metaboliteisoleucine                   0.273  0.273  0.273
## chowYellow (C8):metaboliteLACTIC                       0.275  0.275  0.275
## chowYellow (C8):metaboliteLC even AC total             0.277  0.277  0.277
## chowYellow (C8):metaboliteLC odd AC total              0.267  0.267  0.267
## chowYellow (C8):metaboliteleucine                      0.264  0.264  0.264
## chowYellow (C8):metaboliteMALIC                        0.556  0.278  0.278
## chowYellow (C8):metaboliteMCAC total                   0.275  0.550  0.275
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.260  0.260  0.520
## chowYellow (C8):metaboliteSUCCINIC-2                   0.281  0.281  0.281
## chowYellow (C8):metabolitevaline                       0.276  0.276  0.276
## genotypeKO:chowYellow (C8):metabolitearginine         -0.374 -0.374 -0.374
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.361 -0.361 -0.361
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.363 -0.363 -0.363
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.368 -0.368 -0.368
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.366 -0.366 -0.366
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.365 -0.365 -0.365
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.364 -0.364 -0.364
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.360 -0.360 -0.360
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.359 -0.359 -0.359
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.733 -0.366 -0.366
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.363 -0.726 -0.363
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.356 -0.356 -0.711
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.368 -0.368 -0.368
## genotypeKO:chowYellow (C8):metabolitevaline           -0.366 -0.366 -0.366
##                                                       gKO:SU
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
## genotypeKO:metabolitevaline                            0.500
## chowYellow (C8):metabolitearginine                     0.294
## chowYellow (C8):metaboliteCITRIC                       0.266
## chowYellow (C8):metaboliteFUMARIC                      0.272
## chowYellow (C8):metaboliteglutamine                    0.283
## chowYellow (C8):metaboliteisoleucine                   0.273
## chowYellow (C8):metaboliteLACTIC                       0.275
## chowYellow (C8):metaboliteLC even AC total             0.277
## chowYellow (C8):metaboliteLC odd AC total              0.267
## chowYellow (C8):metaboliteleucine                      0.264
## chowYellow (C8):metaboliteMALIC                        0.278
## chowYellow (C8):metaboliteMCAC total                   0.275
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.260
## chowYellow (C8):metaboliteSUCCINIC-2                   0.562
## chowYellow (C8):metabolitevaline                       0.276
## genotypeKO:chowYellow (C8):metabolitearginine         -0.374
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.361
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.363
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.368
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.366
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.365
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.364
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.360
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.359
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.366
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.363
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.356
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.736
## genotypeKO:chowYellow (C8):metabolitevaline           -0.366
##                                                       gntypKO:mtbltv
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
## chowYellow (C8):metabolitearginine                     0.294        
## chowYellow (C8):metaboliteCITRIC                       0.266        
## chowYellow (C8):metaboliteFUMARIC                      0.272        
## chowYellow (C8):metaboliteglutamine                    0.283        
## chowYellow (C8):metaboliteisoleucine                   0.273        
## chowYellow (C8):metaboliteLACTIC                       0.275        
## chowYellow (C8):metaboliteLC even AC total             0.277        
## chowYellow (C8):metaboliteLC odd AC total              0.267        
## chowYellow (C8):metaboliteleucine                      0.264        
## chowYellow (C8):metaboliteMALIC                        0.278        
## chowYellow (C8):metaboliteMCAC total                   0.275        
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.260        
## chowYellow (C8):metaboliteSUCCINIC-2                   0.281        
## chowYellow (C8):metabolitevaline                       0.553        
## genotypeKO:chowYellow (C8):metabolitearginine         -0.374        
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.361        
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.363        
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.368        
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.366        
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.365        
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.364        
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.360        
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.359        
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.366        
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.363        
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.356        
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.368        
## genotypeKO:chowYellow (C8):metabolitevaline           -0.732        
##                                                       chwYllw(C8):mtbltr
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
## chowYellow (C8):metaboliteCITRIC                       0.534            
## chowYellow (C8):metaboliteFUMARIC                      0.416            
## chowYellow (C8):metaboliteglutamine                    0.460            
## chowYellow (C8):metaboliteisoleucine                   0.497            
## chowYellow (C8):metaboliteLACTIC                       0.490            
## chowYellow (C8):metaboliteLC even AC total             0.494            
## chowYellow (C8):metaboliteLC odd AC total              0.555            
## chowYellow (C8):metaboliteleucine                      0.410            
## chowYellow (C8):metaboliteMALIC                        0.530            
## chowYellow (C8):metaboliteMCAC total                   0.471            
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.524            
## chowYellow (C8):metaboliteSUCCINIC-2                   0.514            
## chowYellow (C8):metabolitevaline                       0.454            
## genotypeKO:chowYellow (C8):metabolitearginine         -0.638            
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.364            
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.284            
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.297            
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.332            
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.328            
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.329            
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.380            
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.278            
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.352            
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.320            
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.361            
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.340            
## genotypeKO:chowYellow (C8):metabolitevaline           -0.304            
##                                                       cY(C8):C cY(C8):F
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
## chowYellow (C8):metaboliteFUMARIC                      0.524           
## chowYellow (C8):metaboliteglutamine                    0.444    0.464  
## chowYellow (C8):metaboliteisoleucine                   0.547    0.526  
## chowYellow (C8):metaboliteLACTIC                       0.513    0.494  
## chowYellow (C8):metaboliteLC even AC total             0.465    0.459  
## chowYellow (C8):metaboliteLC odd AC total              0.516    0.517  
## chowYellow (C8):metaboliteleucine                      0.472    0.434  
## chowYellow (C8):metaboliteMALIC                        0.440    0.532  
## chowYellow (C8):metaboliteMCAC total                   0.445    0.511  
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.504    0.470  
## chowYellow (C8):metaboliteSUCCINIC-2                   0.457    0.454  
## chowYellow (C8):metabolitevaline                       0.448    0.429  
## genotypeKO:chowYellow (C8):metabolitearginine         -0.345   -0.267  
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.679   -0.355  
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.357   -0.671  
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.287   -0.302  
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.372   -0.348  
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.345   -0.329  
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.313   -0.308  
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.355   -0.351  
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.321   -0.295  
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.294   -0.352  
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.304   -0.341  
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.348   -0.323  
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.303   -0.299  
## genotypeKO:chowYellow (C8):metabolitevaline           -0.300   -0.285  
##                                                       chwYllw(C8):mtbltg
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
## chowYellow (C8):metaboliteisoleucine                   0.422            
## chowYellow (C8):metaboliteLACTIC                       0.478            
## chowYellow (C8):metaboliteLC even AC total             0.426            
## chowYellow (C8):metaboliteLC odd AC total              0.535            
## chowYellow (C8):metaboliteleucine                      0.440            
## chowYellow (C8):metaboliteMALIC                        0.390            
## chowYellow (C8):metaboliteMCAC total                   0.460            
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.448            
## chowYellow (C8):metaboliteSUCCINIC-2                   0.374            
## chowYellow (C8):metabolitevaline                       0.368            
## genotypeKO:chowYellow (C8):metabolitearginine         -0.293            
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.301            
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.312            
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.649            
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.286            
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.318            
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.281            
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.366            
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.299            
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.258            
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.307            
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.308            
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.246            
## genotypeKO:chowYellow (C8):metabolitevaline           -0.244            
##                                                       chwYllw(C8):mtblts
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
## chowYellow (C8):metaboliteLACTIC                       0.473            
## chowYellow (C8):metaboliteLC even AC total             0.470            
## chowYellow (C8):metaboliteLC odd AC total              0.555            
## chowYellow (C8):metaboliteleucine                      0.470            
## chowYellow (C8):metaboliteMALIC                        0.481            
## chowYellow (C8):metaboliteMCAC total                   0.451            
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.489            
## chowYellow (C8):metaboliteSUCCINIC-2                   0.478            
## chowYellow (C8):metabolitevaline                       0.491            
## genotypeKO:chowYellow (C8):metabolitearginine         -0.315            
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.370            
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.348            
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.276            
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.663            
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.312            
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.307            
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.380            
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.319            
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.315            
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.292            
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.333            
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.311            
## genotypeKO:chowYellow (C8):metabolitevaline           -0.323            
##                                                       cY(C8):L cY(eAt
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
## chowYellow (C8):metaboliteLC even AC total             0.480         
## chowYellow (C8):metaboliteLC odd AC total              0.471    0.483
## chowYellow (C8):metaboliteleucine                      0.470    0.484
## chowYellow (C8):metaboliteMALIC                        0.458    0.435
## chowYellow (C8):metaboliteMCAC total                   0.455    0.445
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.494    0.480
## chowYellow (C8):metaboliteSUCCINIC-2                   0.476    0.427
## chowYellow (C8):metabolitevaline                       0.476    0.432
## genotypeKO:chowYellow (C8):metabolitearginine         -0.313   -0.317
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.350   -0.316
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.332   -0.310
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.314   -0.276
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.315   -0.316
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.664   -0.321
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.318   -0.661
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.315   -0.333
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.316   -0.328
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.302   -0.289
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.299   -0.298
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.338   -0.330
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.313   -0.282
## genotypeKO:chowYellow (C8):metabolitevaline           -0.318   -0.286
##                                                       cY(oAt
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
## chowYellow (C8):metaboliteleucine                      0.519
## chowYellow (C8):metaboliteMALIC                        0.515
## chowYellow (C8):metaboliteMCAC total                   0.440
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.525
## chowYellow (C8):metaboliteSUCCINIC-2                   0.556
## chowYellow (C8):metabolitevaline                       0.495
## genotypeKO:chowYellow (C8):metabolitearginine         -0.356
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.352
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.349
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.348
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.366
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.314
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.319
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.677
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.352
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.341
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.292
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.360
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.365
## genotypeKO:chowYellow (C8):metabolitevaline           -0.330
##                                                       chwYllw(C8):mtbltl
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
## chowYellow (C8):metaboliteMALIC                        0.482            
## chowYellow (C8):metaboliteMCAC total                   0.489            
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.464            
## chowYellow (C8):metaboliteSUCCINIC-2                   0.459            
## chowYellow (C8):metabolitevaline                       0.447            
## genotypeKO:chowYellow (C8):metabolitearginine         -0.264            
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.321            
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.293            
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.287            
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.311            
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.313            
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.324            
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.354            
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.679            
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.319            
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.326            
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.319            
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.302            
## genotypeKO:chowYellow (C8):metabolitevaline           -0.297            
##                                                       cY(C8):MA cY(C8t
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
## chowYellow (C8):metaboliteMCAC total                   0.451          
## chowYellow (C8):metaboliteMETHYLSUCCINIC               0.461     0.506
## chowYellow (C8):metaboliteSUCCINIC-2                   0.420     0.489
## chowYellow (C8):metabolitevaline                       0.429     0.467
## genotypeKO:chowYellow (C8):metabolitearginine         -0.339    -0.302
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.298    -0.302
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.358    -0.346
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.251    -0.297
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.308    -0.304
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.306    -0.305
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.289    -0.294
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.352    -0.305
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.328    -0.333
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.661    -0.299
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.303    -0.667
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.317    -0.348
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.277    -0.323
## genotypeKO:chowYellow (C8):metabolitevaline           -0.283    -0.311
##                                                       cY(C8):ME cY(C8):S
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
## chowYellow (C8):metaboliteMETHYLSUCCINIC                                
## chowYellow (C8):metaboliteSUCCINIC-2                   0.442            
## chowYellow (C8):metabolitevaline                       0.461     0.401  
## genotypeKO:chowYellow (C8):metabolitearginine         -0.333    -0.330  
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.343    -0.313  
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.314    -0.309  
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.293    -0.246  
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.326    -0.318  
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.328    -0.317  
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.315    -0.286  
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.362    -0.368  
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.313    -0.310  
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.304    -0.278  
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.334    -0.324  
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.685    -0.304  
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.290    -0.656  
## genotypeKO:chowYellow (C8):metabolitevaline           -0.303    -0.273  
##                                                       chwYllw(C8):mtbltv
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
## chowYellow (C8):metaboliteMETHYLSUCCINIC                                
## chowYellow (C8):metaboliteSUCCINIC-2                                    
## chowYellow (C8):metabolitevaline                                        
## genotypeKO:chowYellow (C8):metabolitearginine         -0.291            
## genotypeKO:chowYellow (C8):metaboliteCITRIC           -0.305            
## genotypeKO:chowYellow (C8):metaboliteFUMARIC          -0.290            
## genotypeKO:chowYellow (C8):metaboliteglutamine        -0.240            
## genotypeKO:chowYellow (C8):metaboliteisoleucine       -0.329            
## genotypeKO:chowYellow (C8):metaboliteLACTIC           -0.318            
## genotypeKO:chowYellow (C8):metaboliteLC even AC total -0.288            
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total  -0.338            
## genotypeKO:chowYellow (C8):metaboliteleucine          -0.303            
## genotypeKO:chowYellow (C8):metaboliteMALIC            -0.284            
## genotypeKO:chowYellow (C8):metaboliteMCAC total       -0.313            
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   -0.317            
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       -0.265            
## genotypeKO:chowYellow (C8):metabolitevaline           -0.664            
##                                                       gntypKO:chwYllw(C8):mtbltr
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
## chowYellow (C8):metaboliteMETHYLSUCCINIC                                        
## chowYellow (C8):metaboliteSUCCINIC-2                                            
## chowYellow (C8):metabolitevaline                                                
## genotypeKO:chowYellow (C8):metabolitearginine                                   
## genotypeKO:chowYellow (C8):metaboliteCITRIC            0.516                    
## genotypeKO:chowYellow (C8):metaboliteFUMARIC           0.466                    
## genotypeKO:chowYellow (C8):metaboliteglutamine         0.479                    
## genotypeKO:chowYellow (C8):metaboliteisoleucine        0.498                    
## genotypeKO:chowYellow (C8):metaboliteLACTIC            0.496                    
## genotypeKO:chowYellow (C8):metaboliteLC even AC total  0.488                    
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total   0.525                    
## genotypeKO:chowYellow (C8):metaboliteleucine           0.459                    
## genotypeKO:chowYellow (C8):metaboliteMALIC             0.513                    
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.488                    
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.509                    
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.507                    
## genotypeKO:chowYellow (C8):metabolitevaline            0.480                    
##                                                       gKO:Y(C8):C
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
## chowYellow (C8):metaboliteMETHYLSUCCINIC                         
## chowYellow (C8):metaboliteSUCCINIC-2                             
## chowYellow (C8):metabolitevaline                                 
## genotypeKO:chowYellow (C8):metabolitearginine                    
## genotypeKO:chowYellow (C8):metaboliteCITRIC                      
## genotypeKO:chowYellow (C8):metaboliteFUMARIC           0.505     
## genotypeKO:chowYellow (C8):metaboliteglutamine         0.471     
## genotypeKO:chowYellow (C8):metaboliteisoleucine        0.529     
## genotypeKO:chowYellow (C8):metaboliteLACTIC            0.512     
## genotypeKO:chowYellow (C8):metaboliteLC even AC total  0.487     
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total   0.515     
## genotypeKO:chowYellow (C8):metaboliteleucine           0.486     
## genotypeKO:chowYellow (C8):metaboliteMALIC             0.476     
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.477     
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.505     
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.485     
## genotypeKO:chowYellow (C8):metabolitevaline            0.478     
##                                                       gKO:Y(C8):F
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
## chowYellow (C8):metaboliteMETHYLSUCCINIC                         
## chowYellow (C8):metaboliteSUCCINIC-2                             
## chowYellow (C8):metabolitevaline                                 
## genotypeKO:chowYellow (C8):metabolitearginine                    
## genotypeKO:chowYellow (C8):metaboliteCITRIC                      
## genotypeKO:chowYellow (C8):metaboliteFUMARIC                     
## genotypeKO:chowYellow (C8):metaboliteglutamine         0.483     
## genotypeKO:chowYellow (C8):metaboliteisoleucine        0.509     
## genotypeKO:chowYellow (C8):metaboliteLACTIC            0.499     
## genotypeKO:chowYellow (C8):metaboliteLC even AC total  0.485     
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total   0.512     
## genotypeKO:chowYellow (C8):metaboliteleucine           0.471     
## genotypeKO:chowYellow (C8):metaboliteMALIC             0.515     
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.500     
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.486     
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.483     
## genotypeKO:chowYellow (C8):metabolitevaline            0.466     
##                                                       gntypKO:chwYllw(C8):mtbltg
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
## chowYellow (C8):metaboliteMETHYLSUCCINIC                                        
## chowYellow (C8):metaboliteSUCCINIC-2                                            
## chowYellow (C8):metabolitevaline                                                
## genotypeKO:chowYellow (C8):metabolitearginine                                   
## genotypeKO:chowYellow (C8):metaboliteCITRIC                                     
## genotypeKO:chowYellow (C8):metaboliteFUMARIC                                    
## genotypeKO:chowYellow (C8):metaboliteglutamine                                  
## genotypeKO:chowYellow (C8):metaboliteisoleucine        0.472                    
## genotypeKO:chowYellow (C8):metaboliteLACTIC            0.493                    
## genotypeKO:chowYellow (C8):metaboliteLC even AC total  0.464                    
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total   0.518                    
## genotypeKO:chowYellow (C8):metaboliteleucine           0.472                    
## genotypeKO:chowYellow (C8):metaboliteMALIC             0.451                    
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.484                    
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.477                    
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.448                    
## genotypeKO:chowYellow (C8):metabolitevaline            0.442                    
##                                                       gntypKO:chwYllw(C8):mtblts
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
## chowYellow (C8):metaboliteMETHYLSUCCINIC                                        
## chowYellow (C8):metaboliteSUCCINIC-2                                            
## chowYellow (C8):metabolitevaline                                                
## genotypeKO:chowYellow (C8):metabolitearginine                                   
## genotypeKO:chowYellow (C8):metaboliteCITRIC                                     
## genotypeKO:chowYellow (C8):metaboliteFUMARIC                                    
## genotypeKO:chowYellow (C8):metaboliteglutamine                                  
## genotypeKO:chowYellow (C8):metaboliteisoleucine                                 
## genotypeKO:chowYellow (C8):metaboliteLACTIC            0.489                    
## genotypeKO:chowYellow (C8):metaboliteLC even AC total  0.486                    
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total   0.531                    
## genotypeKO:chowYellow (C8):metaboliteleucine           0.485                    
## genotypeKO:chowYellow (C8):metaboliteMALIC             0.483                    
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.473                    
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.495                    
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.490                    
## genotypeKO:chowYellow (C8):metabolitevaline            0.499                    
##                                                       gKO:Y(C8):L gK(eAt
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
## chowYellow (C8):metaboliteMETHYLSUCCINIC                                
## chowYellow (C8):metaboliteSUCCINIC-2                                    
## chowYellow (C8):metabolitevaline                                        
## genotypeKO:chowYellow (C8):metabolitearginine                           
## genotypeKO:chowYellow (C8):metaboliteCITRIC                             
## genotypeKO:chowYellow (C8):metaboliteFUMARIC                            
## genotypeKO:chowYellow (C8):metaboliteglutamine                          
## genotypeKO:chowYellow (C8):metaboliteisoleucine                         
## genotypeKO:chowYellow (C8):metaboliteLACTIC                             
## genotypeKO:chowYellow (C8):metaboliteLC even AC total  0.491            
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total   0.485       0.497
## genotypeKO:chowYellow (C8):metaboliteleucine           0.485       0.493
## genotypeKO:chowYellow (C8):metaboliteMALIC             0.483       0.472
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.475       0.470
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.497       0.488
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.490       0.470
## genotypeKO:chowYellow (C8):metabolitevaline            0.492       0.469
##                                                       gK(oAt
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
## chowYellow (C8):metaboliteMETHYLSUCCINIC                    
## chowYellow (C8):metaboliteSUCCINIC-2                        
## chowYellow (C8):metabolitevaline                            
## genotypeKO:chowYellow (C8):metabolitearginine               
## genotypeKO:chowYellow (C8):metaboliteCITRIC                 
## genotypeKO:chowYellow (C8):metaboliteFUMARIC                
## genotypeKO:chowYellow (C8):metaboliteglutamine              
## genotypeKO:chowYellow (C8):metaboliteisoleucine             
## genotypeKO:chowYellow (C8):metaboliteLACTIC                 
## genotypeKO:chowYellow (C8):metaboliteLC even AC total       
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total        
## genotypeKO:chowYellow (C8):metaboliteleucine           0.510
## genotypeKO:chowYellow (C8):metaboliteMALIC             0.510
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.473
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.517
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.520
## genotypeKO:chowYellow (C8):metabolitevaline            0.505
##                                                       gntypKO:chwYllw(C8):mtbltl
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
## chowYellow (C8):metaboliteMETHYLSUCCINIC                                        
## chowYellow (C8):metaboliteSUCCINIC-2                                            
## chowYellow (C8):metabolitevaline                                                
## genotypeKO:chowYellow (C8):metabolitearginine                                   
## genotypeKO:chowYellow (C8):metaboliteCITRIC                                     
## genotypeKO:chowYellow (C8):metaboliteFUMARIC                                    
## genotypeKO:chowYellow (C8):metaboliteglutamine                                  
## genotypeKO:chowYellow (C8):metaboliteisoleucine                                 
## genotypeKO:chowYellow (C8):metaboliteLACTIC                                     
## genotypeKO:chowYellow (C8):metaboliteLC even AC total                           
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total                            
## genotypeKO:chowYellow (C8):metaboliteleucine                                    
## genotypeKO:chowYellow (C8):metaboliteMALIC             0.492                    
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.493                    
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.481                    
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.480                    
## genotypeKO:chowYellow (C8):metabolitevaline            0.477                    
##                                                       gKO:Y(C8):MA gKO:(t
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
## chowYellow (C8):metaboliteMETHYLSUCCINIC                                 
## chowYellow (C8):metaboliteSUCCINIC-2                                     
## chowYellow (C8):metabolitevaline                                         
## genotypeKO:chowYellow (C8):metabolitearginine                            
## genotypeKO:chowYellow (C8):metaboliteCITRIC                              
## genotypeKO:chowYellow (C8):metaboliteFUMARIC                             
## genotypeKO:chowYellow (C8):metaboliteglutamine                           
## genotypeKO:chowYellow (C8):metaboliteisoleucine                          
## genotypeKO:chowYellow (C8):metaboliteLACTIC                              
## genotypeKO:chowYellow (C8):metaboliteLC even AC total                    
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total                     
## genotypeKO:chowYellow (C8):metaboliteleucine                             
## genotypeKO:chowYellow (C8):metaboliteMALIC                               
## genotypeKO:chowYellow (C8):metaboliteMCAC total        0.477             
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC    0.483        0.498
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.466        0.491
## genotypeKO:chowYellow (C8):metabolitevaline            0.468        0.485
##                                                       gKO:Y(C8):ME
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
## chowYellow (C8):metaboliteMETHYLSUCCINIC                          
## chowYellow (C8):metaboliteSUCCINIC-2                              
## chowYellow (C8):metabolitevaline                                  
## genotypeKO:chowYellow (C8):metabolitearginine                     
## genotypeKO:chowYellow (C8):metaboliteCITRIC                       
## genotypeKO:chowYellow (C8):metaboliteFUMARIC                      
## genotypeKO:chowYellow (C8):metaboliteglutamine                    
## genotypeKO:chowYellow (C8):metaboliteisoleucine                   
## genotypeKO:chowYellow (C8):metaboliteLACTIC                       
## genotypeKO:chowYellow (C8):metaboliteLC even AC total             
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total              
## genotypeKO:chowYellow (C8):metaboliteleucine                      
## genotypeKO:chowYellow (C8):metaboliteMALIC                        
## genotypeKO:chowYellow (C8):metaboliteMCAC total                   
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC               
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2        0.474      
## genotypeKO:chowYellow (C8):metabolitevaline            0.481      
##                                                       gKO:Y(C8):S
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
## chowYellow (C8):metaboliteMETHYLSUCCINIC                         
## chowYellow (C8):metaboliteSUCCINIC-2                             
## chowYellow (C8):metabolitevaline                                 
## genotypeKO:chowYellow (C8):metabolitearginine                    
## genotypeKO:chowYellow (C8):metaboliteCITRIC                      
## genotypeKO:chowYellow (C8):metaboliteFUMARIC                     
## genotypeKO:chowYellow (C8):metaboliteglutamine                   
## genotypeKO:chowYellow (C8):metaboliteisoleucine                  
## genotypeKO:chowYellow (C8):metaboliteLACTIC                      
## genotypeKO:chowYellow (C8):metaboliteLC even AC total            
## genotypeKO:chowYellow (C8):metaboliteLC odd AC total             
## genotypeKO:chowYellow (C8):metaboliteleucine                     
## genotypeKO:chowYellow (C8):metaboliteMALIC                       
## genotypeKO:chowYellow (C8):metaboliteMCAC total                  
## genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC              
## genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2                  
## genotypeKO:chowYellow (C8):metabolitevaline            0.462     
## 
## Standardized Within-Group Residuals:
##          Min           Q1          Med           Q3          Max 
## -13.73899714  -0.10407005   0.02163056   0.18144089   4.43841691 
## 
## Number of Observations: 645
## Number of Groups: 43
```
