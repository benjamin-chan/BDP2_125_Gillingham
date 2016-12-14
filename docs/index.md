---
title: "Metabolomics of very long-chain aclCoA dehydrogenase knockout mice"
date: "2016-12-14 09:35:19"
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
library(broom)
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
## [1] ggplot2_2.1.0     broom_0.4.1       nlme_3.1-128      dplyr_0.5.0      
## [5] magrittr_1.5      readxl_0.1.1      rmarkdown_1.0     knitr_1.14       
## [9] checkpoint_0.3.16
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.7      munsell_0.4.3    mnormt_1.5-4     colorspace_1.2-6
##  [5] lattice_0.20-34  R6_2.1.3         stringr_1.1.0    plyr_1.8.4      
##  [9] tools_3.3.1      parallel_3.3.1   grid_3.3.1       gtable_0.2.0    
## [13] psych_1.6.9      DBI_0.5-1        htmltools_0.3.5  assertthat_0.1  
## [17] digest_0.6.10    tibble_1.2       reshape2_1.4.1   formatR_1.4     
## [21] tidyr_0.6.0      evaluate_0.9     stringi_1.1.1    scales_0.4.0    
## [25] methods_3.3.1    foreign_0.8-67
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
kable(anova(M))
```



|                             | numDF| denDF|      F-value|   p-value|
|:----------------------------|-----:|-----:|------------:|---------:|
|(Intercept)                  |     1|   524| 1.493015e+05| 0.0000000|
|genotype                     |     1|    38| 9.891158e+04| 0.0000000|
|activity                     |     1|    38| 6.313460e+04| 0.0000000|
|metabolite                   |    14|   524| 2.824393e+06| 0.0000000|
|genotype:activity            |     1|    38| 6.335000e-04| 0.9800518|
|genotype:metabolite          |    14|   524| 3.382267e+00| 0.0000297|
|activity:metabolite          |    14|   524| 7.714754e-01| 0.7006795|
|genotype:activity:metabolite |    14|   524| 1.786415e+00| 0.0375486|

```r
kable(tidy(M, effects = "fixed"))
```



|term                                                 |   estimate| std.error|  statistic|   p.value|
|:----------------------------------------------------|----------:|---------:|----------:|---------:|
|(Intercept)                                          |  0.7800250| 0.2429359|  3.2108260| 0.0014047|
|genotypeKO                                           |  0.1111659| 0.3348671|  0.3319702| 0.7417359|
|activityExercise                                     |  0.2776343| 0.3275776|  0.8475375| 0.4020024|
|metabolitearginine                                   | -0.3534360| 0.3336614| -1.0592657| 0.2899668|
|metaboliteCITRIC                                     |  0.2699043| 0.3422640|  0.7885853| 0.4307111|
|metaboliteFUMARIC                                    | -1.7086321| 0.3422640| -4.9921470| 0.0000008|
|metaboliteglutamine                                  |  1.2558975| 0.3336614|  3.7639882| 0.0001861|
|metaboliteisoleucine                                 | -0.5368730| 0.3336614| -1.6090355| 0.1082108|
|metaboliteLACTIC                                     |  1.6755101| 0.3422640|  4.8953735| 0.0000013|
|metaboliteLCAC total                                 |  0.1169499| 0.3336614|  0.3505047| 0.7261008|
|metaboliteleucine                                    | -0.3706435| 0.3336614| -1.1108374| 0.2671478|
|metaboliteMALIC                                      |  0.1387884| 0.3422640|  0.4055011| 0.6852747|
|metaboliteMCAC Total                                 | -1.2621611| 0.3336614| -3.7827607| 0.0001730|
|metaboliteMETHYLSUCCINIC                             | -2.0968883| 0.3422640| -6.1265233| 0.0000000|
|metabolitePYRUVIC_P2P                                | -0.9377425| 0.3422640| -2.7398224| 0.0063568|
|metaboliteSUCCINIC-2                                 |  0.3515290| 0.3422640|  1.0270699| 0.3048614|
|metabolitevaline                                     | -0.3999569| 0.3336614| -1.1986911| 0.2311899|
|genotypeKO:activityExercise                          | -0.2725876| 0.4271739| -0.6381185| 0.5272236|
|genotypeKO:metabolitearginine                        | -0.0845930| 0.4655750| -0.1816958| 0.8558917|
|genotypeKO:metaboliteCITRIC                          | -0.1808884| 0.4717783| -0.3834183| 0.7015652|
|genotypeKO:metaboliteFUMARIC                         | -0.1459773| 0.4717783| -0.3094192| 0.7571257|
|genotypeKO:metaboliteglutamine                       | -0.0934134| 0.4655750| -0.2006409| 0.8410572|
|genotypeKO:metaboliteisoleucine                      | -0.1551589| 0.4655750| -0.3332630| 0.7390691|
|genotypeKO:metaboliteLACTIC                          | -0.2130032| 0.4717783| -0.4514900| 0.6518232|
|genotypeKO:metaboliteLCAC total                      |  0.2585912| 0.4655750|  0.5554233| 0.5788422|
|genotypeKO:metaboliteleucine                         | -0.1519655| 0.4655750| -0.3264038| 0.7442492|
|genotypeKO:metaboliteMALIC                           | -0.1775990| 0.4717783| -0.3764459| 0.7067378|
|genotypeKO:metaboliteMCAC Total                      | -0.2590161| 0.4655750| -0.5563359| 0.5782188|
|genotypeKO:metaboliteMETHYLSUCCINIC                  | -0.3177490| 0.4717783| -0.6735132| 0.5009178|
|genotypeKO:metabolitePYRUVIC_P2P                     | -2.5513187| 0.4717783| -5.4078758| 0.0000001|
|genotypeKO:metaboliteSUCCINIC-2                      | -0.1502577| 0.4717783| -0.3184922| 0.7502385|
|genotypeKO:metabolitevaline                          | -0.1344774| 0.4655750| -0.2888414| 0.7728169|
|activityExercise:metabolitearginine                  | -0.5175138| 0.4551654| -1.1369796| 0.2560665|
|activityExercise:metaboliteCITRIC                    | -0.2035335| 0.4615087| -0.4410177| 0.6593821|
|activityExercise:metaboliteFUMARIC                   | -0.2704327| 0.4615087| -0.5859754| 0.5581444|
|activityExercise:metaboliteglutamine                 | -0.3376519| 0.4551654| -0.7418225| 0.4585272|
|activityExercise:metaboliteisoleucine                | -0.2834567| 0.4551654| -0.6227554| 0.5337162|
|activityExercise:metaboliteLACTIC                    | -0.2581078| 0.4615087| -0.5592697| 0.5762166|
|activityExercise:metaboliteLCAC total                | -0.4170637| 0.4551654| -0.9162905| 0.3599360|
|activityExercise:metaboliteleucine                   | -0.3160980| 0.4551654| -0.6944685| 0.4876962|
|activityExercise:metaboliteMALIC                     | -0.1439485| 0.4615087| -0.3119086| 0.7552340|
|activityExercise:metaboliteMCAC Total                | -0.3739947| 0.4551654| -0.8216678| 0.4116397|
|activityExercise:metaboliteMETHYLSUCCINIC            | -0.3997036| 0.4615087| -0.8660804| 0.3868424|
|activityExercise:metabolitePYRUVIC_P2P               | -1.3300154| 0.4615087| -2.8818860| 0.0041152|
|activityExercise:metaboliteSUCCINIC-2                | -0.3087569| 0.4615087| -0.6690164| 0.5037797|
|activityExercise:metabolitevaline                    | -0.2501754| 0.4551654| -0.5496363| 0.5828030|
|genotypeKO:activityExercise:metabolitearginine       |  0.2727438| 0.5742045|  0.4749941| 0.6349891|
|genotypeKO:activityExercise:metaboliteCITRIC         |  0.0996485| 0.6048647|  0.1647452| 0.8692081|
|genotypeKO:activityExercise:metaboliteFUMARIC        |  0.2290113| 0.6014938|  0.3807376| 0.7035523|
|genotypeKO:activityExercise:metaboliteglutamine      |  0.0838461| 0.5947977|  0.1409657| 0.8879512|
|genotypeKO:activityExercise:metaboliteisoleucine     |  0.2896583| 0.5969625|  0.4852202| 0.6277229|
|genotypeKO:activityExercise:metaboliteLACTIC         |  0.2524058| 0.5989676|  0.4214014| 0.6736348|
|genotypeKO:activityExercise:metaboliteLCAC total     |  0.3834262| 0.6047703|  0.6340031| 0.5263555|
|genotypeKO:activityExercise:metaboliteleucine        | -0.4830162| 0.6010027| -0.8036839| 0.4219440|
|genotypeKO:activityExercise:metaboliteMALIC          |  0.1327393| 0.6047662|  0.2194886| 0.8263549|
|genotypeKO:activityExercise:metaboliteMCAC Total     |  0.2062359| 0.6050828|  0.3408392| 0.7333613|
|genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC |  0.5375277| 0.6001742|  0.8956195| 0.3708672|
|genotypeKO:activityExercise:metabolitePYRUVIC_P2P    |  2.1895025| 0.6137569|  3.5673774| 0.0003936|
|genotypeKO:activityExercise:metaboliteSUCCINIC-2     |  0.0678668| 0.6165921|  0.1100676| 0.9123979|
|genotypeKO:activityExercise:metabolitevaline         |  0.5488350| 0.5997617|  0.9150884| 0.3605661|


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
kable(anova(M))
```



|                         | numDF| denDF|      F-value|   p-value|
|:------------------------|-----:|-----:|------------:|---------:|
|(Intercept)              |     1|   546| 1.110882e+07| 0.0000000|
|genotype                 |     1|    39| 8.586642e+11| 0.0000000|
|chow                     |     1|    39| 4.026049e+10| 0.0000000|
|metabolite               |    14|   546| 2.897170e+11| 0.0000000|
|genotype:chow            |     1|    39| 1.647336e+00| 0.2068961|
|genotype:metabolite      |    14|   546| 3.537650e+00| 0.0000135|
|chow:metabolite          |    14|   546| 3.383727e+00| 0.0000289|
|genotype:chow:metabolite |    14|   546| 1.312332e+00| 0.1948527|

```r
kable(tidy(M, effects = "fixed"))
```



|term                                                  |   estimate| std.error|  statistic|   p.value|
|:-----------------------------------------------------|----------:|---------:|----------:|---------:|
|(Intercept)                                           |  0.2557123| 0.1585866|  1.6124461| 0.1074427|
|genotypeKO                                            |  0.6626366| 0.2242753|  2.9545680| 0.0052861|
|chowYellow (C8)                                       |  0.7549455| 0.2009892|  3.7561504| 0.0005635|
|metabolitearginine                                    |  0.1236579| 0.2221485|  0.5566453| 0.5779979|
|metaboliteCITRIC                                      |  0.0395968| 0.2221485|  0.1782448| 0.8585968|
|metaboliteFUMARIC                                     | -1.0066586| 0.2221485| -4.5314680| 0.0000072|
|metaboliteglutamine                                   |  1.7844750| 0.2221485|  8.0328039| 0.0000000|
|metaboliteisoleucine                                  | -0.0276526| 0.2221485| -0.1244782| 0.9009824|
|metaboliteLACTIC                                      |  2.1516434| 0.2221485|  9.6856102| 0.0000000|
|metaboliteLC even AC total                            |  0.4067190| 0.2221485|  1.8308434| 0.0676688|
|metaboliteLC odd AC total                             | -0.5647188| 0.2221485| -2.5420786| 0.0112947|
|metaboliteleucine                                     |  0.0958975| 0.2221485|  0.4316820| 0.6661429|
|metaboliteMALIC                                       |  0.7454015| 0.2221485|  3.3554206| 0.0008475|
|metaboliteMCAC total                                  | -0.7176324| 0.2221485| -3.2304179| 0.0013105|
|metaboliteMETHYLSUCCINIC                              | -1.4975941| 0.2221485| -6.7414110| 0.0000000|
|metaboliteSUCCINIC-2                                  |  0.9674087| 0.2221485|  4.3547845| 0.0000159|
|metabolitevaline                                      |  0.1464322| 0.2221485|  0.6591639| 0.5100684|
|genotypeKO:chowYellow (C8)                            | -0.6978098| 0.3056066| -2.2833600| 0.0279394|
|genotypeKO:metabolitearginine                         | -0.6736076| 0.3141654| -2.1441179| 0.0324643|
|genotypeKO:metaboliteCITRIC                           |  0.0088516| 0.3141654|  0.0281751| 0.9775328|
|genotypeKO:metaboliteFUMARIC                          | -0.9692841| 0.3141654| -3.0852672| 0.0021367|
|genotypeKO:metaboliteglutamine                        | -0.7667286| 0.3141654| -2.4405257| 0.0149827|
|genotypeKO:metaboliteisoleucine                       | -0.7247245| 0.3141654| -2.3068249| 0.0214380|
|genotypeKO:metaboliteLACTIC                           | -0.8237407| 0.3141654| -2.6219972| 0.0089855|
|genotypeKO:metaboliteLC even AC total                 | -0.1083540| 0.3141654| -0.3448948| 0.7303063|
|genotypeKO:metaboliteLC odd AC total                  |  0.0839471| 0.3141654|  0.2672068| 0.7894109|
|genotypeKO:metaboliteleucine                          | -0.6996743| 0.3141654| -2.2270892| 0.0263483|
|genotypeKO:metaboliteMALIC                            | -0.7985721| 0.3141654| -2.5418847| 0.0113009|
|genotypeKO:metaboliteMCAC total                       | -0.5409875| 0.3141654| -1.7219833| 0.0856389|
|genotypeKO:metaboliteMETHYLSUCCINIC                   | -0.9184289| 0.3141654| -2.9233932| 0.0036060|
|genotypeKO:metaboliteSUCCINIC-2                       | -0.7749538| 0.3141654| -2.4667068| 0.0139421|
|genotypeKO:metabolitevaline                           | -0.7162648| 0.3141654| -2.2798973| 0.0229986|
|chowYellow (C8):metabolitearginine                    | -1.2626957| 0.2673925| -4.7222552| 0.0000030|
|chowYellow (C8):metaboliteCITRIC                      | -0.5439097| 0.2949119| -1.8443124| 0.0656792|
|chowYellow (C8):metaboliteFUMARIC                     | -0.9714318| 0.2891552| -3.3595509| 0.0008352|
|chowYellow (C8):metaboliteglutamine                   | -0.9631569| 0.2771722| -3.4749409| 0.0005516|
|chowYellow (C8):metaboliteisoleucine                  | -1.1705246| 0.2876369| -4.0694523| 0.0000541|
|chowYellow (C8):metaboliteLACTIC                      | -0.8315275| 0.2854145| -2.9134029| 0.0037215|
|chowYellow (C8):metaboliteLC even AC total            | -0.7875069| 0.2840107| -2.7728067| 0.0057474|
|chowYellow (C8):metaboliteLC odd AC total             | -0.8928525| 0.2937578| -3.0394168| 0.0024839|
|chowYellow (C8):metaboliteleucine                     | -1.0879909| 0.2979676| -3.6513736| 0.0002859|
|chowYellow (C8):metaboliteMALIC                       | -0.8293624| 0.2825788| -2.9349773| 0.0034761|
|chowYellow (C8):metaboliteMCAC total                  | -0.9477404| 0.2855944| -3.3184832| 0.0009654|
|chowYellow (C8):metaboliteMETHYLSUCCINIC              | -0.3178135| 0.3020975| -1.0520231| 0.2932543|
|chowYellow (C8):metaboliteSUCCINIC-2                  | -0.7365276| 0.2795860| -2.6343507| 0.0086691|
|chowYellow (C8):metabolitevaline                      | -0.5978008| 0.2843059| -2.1026676| 0.0359520|
|genotypeKO:chowYellow (C8):metabolitearginine         |  1.1183598| 0.4199741|  2.6629259| 0.0079750|
|genotypeKO:chowYellow (C8):metaboliteCITRIC           |  0.5398332| 0.4351838|  1.2404716| 0.2153338|
|genotypeKO:chowYellow (C8):metaboliteFUMARIC          |  0.9596491| 0.4323726|  2.2194956| 0.0268631|
|genotypeKO:chowYellow (C8):metaboliteglutamine        |  0.9859728| 0.4269021|  2.3095995| 0.0212825|
|genotypeKO:chowYellow (C8):metaboliteisoleucine       |  1.0857875| 0.4296033|  2.5274192| 0.0117715|
|genotypeKO:chowYellow (C8):metaboliteLACTIC           |  0.6497277| 0.4302308|  1.5101839| 0.1315750|
|genotypeKO:chowYellow (C8):metaboliteLC even AC total |  0.7968450| 0.4309968|  1.8488423| 0.0650210|
|genotypeKO:chowYellow (C8):metaboliteLC odd AC total  |  0.1361594| 0.4360711|  0.3122412| 0.7549764|
|genotypeKO:chowYellow (C8):metaboliteleucine          |  1.0353980| 0.4380953|  2.3634079| 0.0184571|
|genotypeKO:chowYellow (C8):metaboliteMALIC            |  0.6371135| 0.4288386|  1.4856718| 0.1379431|
|genotypeKO:chowYellow (C8):metaboliteMCAC total       |  0.8934399| 0.4326781|  2.0649068| 0.0394034|
|genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   |  0.9976601| 0.4415882|  2.2592546| 0.0242610|
|genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       |  0.6302176| 0.4266922|  1.4769839| 0.1402565|
|genotypeKO:chowYellow (C8):metabolitevaline           |  0.7417809| 0.4293746|  1.7275843| 0.0846282|
