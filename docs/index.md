---
title: "Metabolomics of very long-chain aclCoA dehydrogenase knockout mice"
date: "2016-11-23 14:10:51"
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
library(nlme)
library(ggplot2)
library(knitr)
```

Reproducibility steps.


```r
sessionInfo()
```

```
## R version 3.3.1 (2016-06-21)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 7 x64 (build 7601) Service Pack 1
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] rmarkdown_1.0       checkpoint_0.3.16   knitr_1.13         
## [4] ggplot2_2.1.0       nlme_3.1-128        dplyr_0.5.0        
## [7] magrittr_1.5        readxl_0.1.1        RevoUtilsMath_8.0.3
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.6      munsell_0.4.3    colorspace_1.2-6 lattice_0.20-33 
##  [5] R6_2.1.2         stringr_1.0.0    highr_0.6        plyr_1.8.4      
##  [9] tools_3.3.1      grid_3.3.1       gtable_0.2.0     DBI_0.5         
## [13] htmltools_0.3.5  lazyeval_0.2.0   assertthat_0.1   digest_0.6.10   
## [17] tibble_1.0       reshape2_1.4.1   formatR_1.4      evaluate_0.9    
## [21] labeling_0.3     stringi_1.1.1    RevoUtils_10.0.1 scales_0.4.0
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
## Source: local data frame [6 x 8]
## 
##       id genotype activity    chow metabolite_type
##    <chr>   <fctr>   <fctr>  <fctr>           <chr>
## 1 1208B   Pos/Pos Exercise Regular  acylcarnitines
## 2 1208B   Pos/Pos Exercise Regular  acylcarnitines
## 3 1208B   Pos/Pos Exercise Regular  acylcarnitines
## 4 1208B   Pos/Pos Exercise Regular  acylcarnitines
## 5 1208B   Pos/Pos Exercise Regular  acylcarnitines
## 6 1208B   Pos/Pos Exercise Regular  acylcarnitines
## Variables not shown: metabolite <chr>, value <dbl>, logValue <dbl>.
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
## Source: local data frame [6 x 8]
## 
##      id genotype activity       chow metabolite_type
##   <chr>   <fctr>   <fctr>     <fctr>           <chr>
## 1  1120  Neg/Neg Exercise White (C7)  acylcarnitines
## 2  1120  Neg/Neg Exercise White (C7)  acylcarnitines
## 3  1120  Neg/Neg Exercise White (C7)  acylcarnitines
## 4  1120  Neg/Neg Exercise White (C7)  acylcarnitines
## 5  1120  Neg/Neg Exercise White (C7)  acylcarnitines
## 6  1120  Neg/Neg Exercise White (C7)  acylcarnitines
## Variables not shown: metabolite <chr>, value <dbl>, logValue <dbl>.
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
##          
##           acylcarnitines amino acids organic acids
##   Neg/Neg           1659        1071           714
##   Pos/Pos           1659        1071           680
##           
##            acylcarnitines amino acids organic acids
##   Rest                  0        1020           646
##   Exercise           3318        1122           748
##          
##           acylcarnitines amino acids organic acids
##   Regular           3318        2142          1394
```

```r
tableFixed(D2)
```

```
##          
##           acylcarnitines Amino Acids organic acids
##   Neg/Neg           1848        1100           770
##   Pos/Pos           1848        1100           770
##           
##            acylcarnitines Amino Acids organic acids
##   Exercise           3696        2200          1540
##              
##               acylcarnitines Amino Acids organic acids
##   White (C7)            1848        1100           770
##   Yellow (C8)           1848        1100           770
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
## StdDev: 5.139997e-05 1.373658
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
Use `corCompSymm`, *compound symmetry structure corresponding to a constant correlation*.


```r
cs <-
  corCompSymm(form = random, fixed = FALSE) %>%
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
##   23835.58 23910.73 -11906.79
## 
## Random effects:
##  Formula: ~1 | id
##          (Intercept) Residual
## StdDev: 0.0001492959 1.373392
## 
## Correlation Structure: Compound symmetry
##  Formula: ~1 | id 
##  Parameter estimate(s):
##          Rho 
## -0.004067683 
## Fixed effects: list(fixed) 
##                                                   Value  Std.Error   DF
## (Intercept)                                  -1.4723218 0.05269083 6801
## genotypePos/Pos                               0.0609111 0.07532329   45
## activityExercise                             -0.0026047 0.04766522 6801
## metabolite_typeamino acids                    1.5001335 0.05594654 6801
## metabolite_typeorganic acids                  1.1048377 0.06286770 6801
## genotypePos/Pos:activityExercise              0.0184061 0.06669116 6801
## genotypePos/Pos:metabolite_typeamino acids   -0.0638759 0.07996824 6801
## genotypePos/Pos:metabolite_typeorganic acids  0.0224827 0.08973140 6801
##                                                 t-value p-value
## (Intercept)                                  -27.942657  0.0000
## genotypePos/Pos                                0.808662  0.4230
## activityExercise                              -0.054646  0.9564
## metabolite_typeamino acids                    26.813696  0.0000
## metabolite_typeorganic acids                  17.574013  0.0000
## genotypePos/Pos:activityExercise               0.275989  0.7826
## genotypePos/Pos:metabolite_typeamino acids    -0.798766  0.4245
## genotypePos/Pos:metabolite_typeorganic acids   0.250555  0.8022
##  Correlation: 
##                                              (Intr) gntP/P actvtE
## genotypePos/Pos                              -0.701              
## activityExercise                             -0.854  0.597       
## metabolite_typeamino acids                   -0.668  0.466  0.373
## metabolite_typeorganic acids                 -0.587  0.410  0.328
## genotypePos/Pos:activityExercise              0.611 -0.855 -0.715
## genotypePos/Pos:metabolite_typeamino acids    0.466 -0.676 -0.261
## genotypePos/Pos:metabolite_typeorganic acids  0.410 -0.582 -0.230
##                                              mtblt_typma mtblt_typra
## genotypePos/Pos                                                     
## activityExercise                                                    
## metabolite_typeamino acids                                          
## metabolite_typeorganic acids                  0.379                 
## genotypePos/Pos:activityExercise             -0.267      -0.235     
## genotypePos/Pos:metabolite_typeamino acids   -0.693      -0.259     
## genotypePos/Pos:metabolite_typeorganic acids -0.259      -0.695     
##                                              gP/P:E gntypPs/Ps:mtblt_typma
## genotypePos/Pos                                                           
## activityExercise                                                          
## metabolite_typeamino acids                                                
## metabolite_typeorganic acids                                              
## genotypePos/Pos:activityExercise                                          
## genotypePos/Pos:metabolite_typeamino acids    0.376                       
## genotypePos/Pos:metabolite_typeorganic acids  0.318  0.385                
## 
## Standardized Within-Group Residuals:
##        Min         Q1        Med         Q3        Max 
## -5.0934659 -0.5157460  0.0635438  0.4954288  3.8418980 
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
## StdDev: 9.091496e-05 1.407474
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
Use `corCompSymm`, *compound symmetry structure corresponding to a constant correlation*.


```r
cs <-
  corCompSymm(form = random, fixed = FALSE) %>%
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
##   26233.52 26309.57 -13105.76
## 
## Random effects:
##  Formula: ~1 | id
##          (Intercept) Residual
## StdDev: 0.0002080774 1.407446
## 
## Correlation Structure: Compound symmetry
##  Formula: ~1 | id 
##  Parameter estimate(s):
##          Rho 
## -0.001858116 
## Fixed effects: list(fixed) 
##                                                   Value  Std.Error   DF
## (Intercept)                                  -1.3195677 0.03418591 7388
## genotypePos/Pos                              -0.1229449 0.05186767   40
## chowYellow (C8)                              -0.1348985 0.03962217   40
## metabolite_typeAmino Acids                    1.4723764 0.04954797 7388
## metabolite_typeorganic acids                  1.0586409 0.05650869 7388
## genotypePos/Pos:chowYellow (C8)               0.0664999 0.05951222   40
## genotypePos/Pos:metabolite_typeAmino Acids    0.0409294 0.07211022 7388
## genotypePos/Pos:metabolite_typeorganic acids  0.0865997 0.08171472 7388
##                                                t-value p-value
## (Intercept)                                  -38.59975  0.0000
## genotypePos/Pos                               -2.37036  0.0227
## chowYellow (C8)                               -3.40462  0.0015
## metabolite_typeAmino Acids                    29.71618  0.0000
## metabolite_typeorganic acids                  18.73412  0.0000
## genotypePos/Pos:chowYellow (C8)                1.11742  0.2705
## genotypePos/Pos:metabolite_typeAmino Acids     0.56759  0.5703
## genotypePos/Pos:metabolite_typeorganic acids   1.05978  0.2893
##  Correlation: 
##                                              (Intr) gntP/P cY(C8) mtb_AA
## genotypePos/Pos                              -0.658                     
## chowYellow (C8)                              -0.546  0.362              
## metabolite_typeAmino Acids                   -0.508  0.330  0.002       
## metabolite_typeorganic acids                 -0.441  0.287 -0.002  0.231
## genotypePos/Pos:chowYellow (C8)               0.363 -0.570 -0.667  0.000
## genotypePos/Pos:metabolite_typeAmino Acids    0.347 -0.499 -0.003 -0.681
## genotypePos/Pos:metabolite_typeorganic acids  0.304 -0.434  0.000 -0.157
##                                              mtbl_a gP/P:( gP/P:A
## genotypePos/Pos                                                  
## chowYellow (C8)                                                  
## metabolite_typeAmino Acids                                       
## metabolite_typeorganic acids                                     
## genotypePos/Pos:chowYellow (C8)               0.003              
## genotypePos/Pos:metabolite_typeAmino Acids   -0.157  0.001       
## genotypePos/Pos:metabolite_typeorganic acids -0.690 -0.003  0.269
## 
## Standardized Within-Group Residuals:
##         Min          Q1         Med          Q3         Max 
## -5.04991618 -0.47538509  0.05242169  0.51131403  4.29105732 
## 
## Number of Observations: 7436
## Number of Groups: 44
```
