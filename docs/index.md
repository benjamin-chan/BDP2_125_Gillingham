---
title: "Metabolomics of very long-chain aclCoA dehydrogenase knockout mice"
date: "2016-12-12 15:28:01"
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
D2 <- L2[["data"]]
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



|              |       Min.| X1st.Qu.| Median|   Mean| X3rd.Qu.|    Max.|
|:-------------|----------:|--------:|------:|------:|--------:|-------:|
|nominal       |  0.0000001|  1.18900|  3.690| 23.980|   12.850| 408.800|
|log-transform | -6.9550000|  0.07533|  0.567|  0.526|    1.109|   2.612|

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
table(D1$genotype, D1$metabolite_type)
```

```
##     
##      Acylcarnitines Amino acids Organic acids
##   WT             42         105           160
##   KO             42         105           168
```

```r
table(D1$activity, D1$metabolite_type)
```

```
##           
##            Acylcarnitines Amino acids Organic acids
##   Rest                 40         100           152
##   Exercise             44         110           176
```

```r
table(D1$chow, D1$metabolite_type)
```

```
##          
##           Acylcarnitines Amino acids Organic acids
##   Regular             84         210           328
```

```r
table(D1$metabolite, D1$metabolite_type)
```

```
##                   
##                    Acylcarnitines Amino acids Organic acids
##   3-HYDROXYBUTYRIC              0           0            41
##   arginine                      0          42             0
##   CITRIC                        0           0            41
##   FUMARIC                       0           0            41
##   glutamine                     0          42             0
##   isoleucine                    0          42             0
##   LACTIC                        0           0            41
##   LCAC total                   42           0             0
##   leucine                       0          42             0
##   MALIC                         0           0            41
##   MCAC Total                   42           0             0
##   METHYLSUCCINIC                0           0            41
##   PYRUVIC_P2P                   0           0            41
##   SUCCINIC-2                    0           0            41
##   valine                        0          42             0
```


```r
table(D2$genotype, D2$metabolite_type)
```

```
##     
##      Acylcarnitines Amino acids Organic acids
##   WT             88         110           154
##   KO             84         105           147
```

```r
table(D2$activity, D2$metabolite_type)
```

```
##           
##            Acylcarnitines Amino acids Organic acids
##   Exercise            172         215           301
```

```r
table(D2$chow, D2$metabolite_type)
```

```
##              
##               Acylcarnitines Amino acids Organic acids
##   White (C7)              88         110           154
##   Yellow (C8)             84         105           147
```

```r
table(D2$metabolite)
```

```
## 
## 3-HYDROXYBUTYRIC         arginine           CITRIC          FUMARIC 
##               43               43               43               43 
##        glutamine       isoleucine           LACTIC LC even AC total 
##               43               43               43               43 
##  LC odd AC total       LCAC total          leucine            MALIC 
##               43               43               43               43 
##       MCAC total   METHYLSUCCINIC       SUCCINIC-2           valine 
##               43               43               43               43
```

```r
table(D2$metabolite, D2$metabolite_type)
```

```
##                   
##                    Acylcarnitines Amino acids Organic acids
##   3-HYDROXYBUTYRIC              0           0            43
##   arginine                      0          43             0
##   CITRIC                        0           0            43
##   FUMARIC                       0           0            43
##   glutamine                     0          43             0
##   isoleucine                    0          43             0
##   LACTIC                        0           0            43
##   LC even AC total             43           0             0
##   LC odd AC total              43           0             0
##   LCAC total                   43           0             0
##   leucine                       0          43             0
##   MALIC                         0           0            43
##   MCAC total                   43           0             0
##   METHYLSUCCINIC                0           0            43
##   SUCCINIC-2                    0           0            43
##   valine                        0          43             0
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
##   2158.771 2202.971 -1069.385
## 
## Random effects:
##  Formula: ~1 | id
##          (Intercept) Residual
## StdDev: 4.239754e-05 1.341349
## 
## Fixed effects: list(fixed) 
##                                              Value Std.Error  DF
## (Intercept)                              0.1890360 0.2220616 576
## genotypeKO                               0.0975869 0.3137401  38
## activityExercise                        -0.0827994 0.1535940  38
## metabolite_typeAmino acids               0.5201771 0.2448957 576
## metabolite_typeOrganic acids             0.3002810 0.2325936 576
## genotypeKO:activityExercise              0.0500204 0.2156160  38
## genotypeKO:metabolite_typeAmino acids   -0.1732246 0.3463348 576
## genotypeKO:metabolite_typeOrganic acids -0.4103840 0.3280974 576
##                                            t-value p-value
## (Intercept)                              0.8512776  0.3950
## genotypeKO                               0.3110439  0.7575
## activityExercise                        -0.5390798  0.5930
## metabolite_typeAmino acids               2.1240761  0.0341
## metabolite_typeOrganic acids             1.2910116  0.1972
## genotypeKO:activityExercise              0.2319883  0.8178
## genotypeKO:metabolite_typeAmino acids   -0.5001652  0.6171
## genotypeKO:metabolite_typeOrganic acids -1.2507994  0.2115
##  Correlation: 
##                                         (Intr) gntyKO actvtE mtb_Aa mtb_Oa
## genotypeKO                              -0.708                            
## activityExercise                        -0.362  0.256                     
## metabolite_typeAmino acids              -0.788  0.558  0.000              
## metabolite_typeOrganic acids            -0.823  0.583 -0.017  0.752       
## genotypeKO:activityExercise              0.258 -0.360 -0.712  0.000  0.012
## genotypeKO:metabolite_typeAmino acids    0.557 -0.788  0.000 -0.707 -0.532
## genotypeKO:metabolite_typeOrganic acids  0.584 -0.829  0.012 -0.533 -0.709
##                                         gnKO:E gKO:_Aa
## genotypeKO                                            
## activityExercise                                      
## metabolite_typeAmino acids                            
## metabolite_typeOrganic acids                          
## genotypeKO:activityExercise                           
## genotypeKO:metabolite_typeAmino acids    0.000        
## genotypeKO:metabolite_typeOrganic acids -0.009  0.754 
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -5.63274292 -0.36589546 -0.09890158  0.58652147  1.76166847 
## 
## Number of Observations: 622
## Number of Groups: 42
```

Estimate model.
Specify the correlation structure using `cs`.
Use `corSymm`, *general correlation matrix, with no additional structure*.


```r
rm(M)
ctrl <- lmeControl(opt = "optim")
cs <-
  corSymm(form = random, fixed = TRUE) %>%
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
summary(M)
```

```
## Linear mixed-effects model fit by REML
##  Data: D1 
##        AIC      BIC    logLik
##   2158.775 2202.975 -1069.388
## 
## Random effects:
##  Formula: ~1 | id
##         (Intercept) Residual
## StdDev: 0.004510822 1.341347
## 
## Correlation Structure: General
##  Formula: ~1 | id 
##  Parameter estimate(s):
##  Correlation: 
##    1 2 3 4 5 6 7 8 9 10 11 12 13 14
## 2  0                               
## 3  0 0                             
## 4  0 0 0                           
## 5  0 0 0 0                         
## 6  0 0 0 0 0                       
## 7  0 0 0 0 0 0                     
## 8  0 0 0 0 0 0 0                   
## 9  0 0 0 0 0 0 0 0                 
## 10 0 0 0 0 0 0 0 0 0               
## 11 0 0 0 0 0 0 0 0 0 0             
## 12 0 0 0 0 0 0 0 0 0 0  0          
## 13 0 0 0 0 0 0 0 0 0 0  0  0       
## 14 0 0 0 0 0 0 0 0 0 0  0  0  0    
## 15 0 0 0 0 0 0 0 0 0 0  0  0  0  0 
## Fixed effects: list(fixed) 
##                                              Value Std.Error  DF
## (Intercept)                              0.1890362 0.2220659 576
## genotypeKO                               0.0975868 0.3137461  38
## activityExercise                        -0.0827997 0.1536066  38
## metabolite_typeAmino acids               0.5201771 0.2448953 576
## metabolite_typeOrganic acids             0.3002812 0.2325933 576
## genotypeKO:activityExercise              0.0500206 0.2156338  38
## genotypeKO:metabolite_typeAmino acids   -0.1732246 0.3463343 576
## genotypeKO:metabolite_typeOrganic acids -0.4103842 0.3280970 576
##                                            t-value p-value
## (Intercept)                              0.8512618  0.3950
## genotypeKO                               0.3110376  0.7575
## activityExercise                        -0.5390370  0.5930
## metabolite_typeAmino acids               2.1240794  0.0341
## metabolite_typeOrganic acids             1.2910139  0.1972
## genotypeKO:activityExercise              0.2319702  0.8178
## genotypeKO:metabolite_typeAmino acids   -0.5001660  0.6171
## genotypeKO:metabolite_typeOrganic acids -1.2508017  0.2115
##  Correlation: 
##                                         (Intr) gntyKO actvtE mtb_Aa mtb_Oa
## genotypeKO                              -0.708                            
## activityExercise                        -0.362  0.256                     
## metabolite_typeAmino acids              -0.788  0.558  0.000              
## metabolite_typeOrganic acids            -0.823  0.583 -0.017  0.752       
## genotypeKO:activityExercise              0.258 -0.360 -0.712  0.000  0.012
## genotypeKO:metabolite_typeAmino acids    0.557 -0.788  0.000 -0.707 -0.532
## genotypeKO:metabolite_typeOrganic acids  0.584 -0.829  0.012 -0.533 -0.709
##                                         gnKO:E gKO:_Aa
## genotypeKO                                            
## activityExercise                                      
## metabolite_typeAmino acids                            
## metabolite_typeOrganic acids                          
## genotypeKO:activityExercise                           
## genotypeKO:metabolite_typeAmino acids    0.000        
## genotypeKO:metabolite_typeOrganic acids -0.009  0.754 
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -5.63270748 -0.36590175 -0.09892581  0.58649943  1.76164765 
## 
## Number of Observations: 622
## Number of Groups: 42
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
##   2039.157 2084.378 -1009.578
## 
## Random effects:
##  Formula: ~1 | id
##          (Intercept) Residual
## StdDev: 3.091411e-05 1.039132
## 
## Fixed effects: list(fixed) 
##                                              Value Std.Error  DF   t-value
## (Intercept)                              0.0918984 0.1238467 641  0.742033
## genotypeKO                               0.4765655 0.1763923  39  2.701737
## chowYellow (C8)                          0.0445021 0.1107718  39  0.401745
## metabolite_typeAmino acids               0.5639168 0.1486160 641  3.794455
## metabolite_typeOrganic acids             0.4143852 0.1388599 641  2.984198
## genotypeKO:chowYellow (C8)              -0.0389779 0.1586011  39 -0.245760
## genotypeKO:metabolite_typeAmino acids   -0.4936399 0.2126621 641 -2.321240
## genotypeKO:metabolite_typeOrganic acids -0.4376502 0.1987016 641 -2.202550
##                                         p-value
## (Intercept)                              0.4583
## genotypeKO                               0.0102
## chowYellow (C8)                          0.6901
## metabolite_typeAmino acids               0.0002
## metabolite_typeOrganic acids             0.0030
## genotypeKO:chowYellow (C8)               0.8072
## genotypeKO:metabolite_typeAmino acids    0.0206
## genotypeKO:metabolite_typeOrganic acids  0.0280
##  Correlation: 
##                                         (Intr) gntyKO cY(C8) mtb_Aa mtb_Oa
## genotypeKO                              -0.702                            
## chowYellow (C8)                         -0.447  0.314                     
## metabolite_typeAmino acids              -0.667  0.468  0.000              
## metabolite_typeOrganic acids            -0.714  0.501  0.000  0.595       
## genotypeKO:chowYellow (C8)               0.312 -0.439 -0.698  0.000  0.000
## genotypeKO:metabolite_typeAmino acids    0.466 -0.670  0.000 -0.699 -0.416
## genotypeKO:metabolite_typeOrganic acids  0.499 -0.717  0.000 -0.416 -0.699
##                                         gKO:Y( gKO:_Aa
## genotypeKO                                            
## chowYellow (C8)                                       
## metabolite_typeAmino acids                            
## metabolite_typeOrganic acids                          
## genotypeKO:chowYellow (C8)                            
## genotypeKO:metabolite_typeAmino acids    0.000        
## genotypeKO:metabolite_typeOrganic acids  0.000  0.595 
## 
## Standardized Within-Group Residuals:
##        Min         Q1        Med         Q3        Max 
## -7.2228197 -0.4887542  0.0200908  0.5899636  2.0259443 
## 
## Number of Observations: 688
## Number of Groups: 43
```

Estimate model.
Specify the correlation structure using `cs`.
Use `corSymm`, *general correlation matrix, with no additional structure*.


```r
rm(M)
cs <-
  corSymm(form = random, fixed = TRUE) %>%
  Initialize(data = D2)
Dim(cs)
```

```
## $N
## [1] 688
## 
## $M
## [1] 43
## 
## $maxLen
## [1] 16
## 
## $sumLenSq
## [1] 11008
## 
## $len
## groups
##       1101       1102       1103       1184       1185       1186 
##         16         16         16         16         16         16 
##       1195       1196       1197       1203 1192/1198B       1176 
##         16         16         16         16         16         16 
##       1177       1179       1180       1190       1191       1193 
##         16         16         16         16         16         16 
##       1194       1199       1200 1192A/1198       1107       1113 
##         16         16         16         16         16         16 
##       1114       1115       1117       1118       1119       1204 
##         16         16         16         16         16         16 
##       1205       1206       1120       1126       1127       1128 
##         16         16         16         16         16         16 
##       1142       1143       1144       1145       1146       1158 
##         16         16         16         16         16         16 
##       1159 
##         16 
## 
## $start
##  [1]  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
## [24] 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42
```

```r
M <- lme(fixed, data = D2, random = random, correlation = cs, control = ctrl)
summary(M)
```

```
## Linear mixed-effects model fit by REML
##  Data: D2 
##        AIC      BIC    logLik
##   2039.161 2084.382 -1009.581
## 
## Random effects:
##  Formula: ~1 | id
##         (Intercept) Residual
## StdDev: 0.003164857 1.039131
## 
## Correlation Structure: General
##  Formula: ~1 | id 
##  Parameter estimate(s):
##  Correlation: 
##    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
## 2  0                                  
## 3  0 0                                
## 4  0 0 0                              
## 5  0 0 0 0                            
## 6  0 0 0 0 0                          
## 7  0 0 0 0 0 0                        
## 8  0 0 0 0 0 0 0                      
## 9  0 0 0 0 0 0 0 0                    
## 10 0 0 0 0 0 0 0 0 0                  
## 11 0 0 0 0 0 0 0 0 0 0                
## 12 0 0 0 0 0 0 0 0 0 0  0             
## 13 0 0 0 0 0 0 0 0 0 0  0  0          
## 14 0 0 0 0 0 0 0 0 0 0  0  0  0       
## 15 0 0 0 0 0 0 0 0 0 0  0  0  0  0    
## 16 0 0 0 0 0 0 0 0 0 0  0  0  0  0  0 
## Fixed effects: list(fixed) 
##                                              Value Std.Error  DF   t-value
## (Intercept)                              0.0918984 0.1238502 641  0.742012
## genotypeKO                               0.4765655 0.1763972  39  2.701661
## chowYellow (C8)                          0.0445021 0.1107799  39  0.401716
## metabolite_typeAmino acids               0.5639168 0.1486158 641  3.794459
## metabolite_typeOrganic acids             0.4143852 0.1388597 641  2.984201
## genotypeKO:chowYellow (C8)              -0.0389779 0.1586127  39 -0.245742
## genotypeKO:metabolite_typeAmino acids   -0.4936399 0.2126619 641 -2.321243
## genotypeKO:metabolite_typeOrganic acids -0.4376502 0.1987013 641 -2.202553
##                                         p-value
## (Intercept)                              0.4584
## genotypeKO                               0.0102
## chowYellow (C8)                          0.6901
## metabolite_typeAmino acids               0.0002
## metabolite_typeOrganic acids             0.0030
## genotypeKO:chowYellow (C8)               0.8072
## genotypeKO:metabolite_typeAmino acids    0.0206
## genotypeKO:metabolite_typeOrganic acids  0.0280
##  Correlation: 
##                                         (Intr) gntyKO cY(C8) mtb_Aa mtb_Oa
## genotypeKO                              -0.702                            
## chowYellow (C8)                         -0.447  0.314                     
## metabolite_typeAmino acids              -0.667  0.468  0.000              
## metabolite_typeOrganic acids            -0.713  0.501  0.000  0.595       
## genotypeKO:chowYellow (C8)               0.312 -0.439 -0.698  0.000  0.000
## genotypeKO:metabolite_typeAmino acids    0.466 -0.670  0.000 -0.699 -0.416
## genotypeKO:metabolite_typeOrganic acids  0.499 -0.717  0.000 -0.416 -0.699
##                                         gKO:Y( gKO:_Aa
## genotypeKO                                            
## chowYellow (C8)                                       
## metabolite_typeAmino acids                            
## metabolite_typeOrganic acids                          
## genotypeKO:chowYellow (C8)                            
## genotypeKO:metabolite_typeAmino acids    0.000        
## genotypeKO:metabolite_typeOrganic acids  0.000  0.595 
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -7.22276211 -0.48875493  0.02009388  0.58995641  2.02592409 
## 
## Number of Observations: 688
## Number of Groups: 43
```
