---
title: "Metabolomics of very long-chain aclCoA dehydrogenase knockout mice"
date: "2016-12-16 11:37:15"
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
library(broom)
library(ggplot2)
library(knitr)
```

Reproducibility steps.


```r
sessionInfo()
```

```
## R version 3.3.2 (2016-10-31)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows >= 8 x64 (build 9200)
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
## [1] ggplot2_2.1.0     broom_0.4.1       nlme_3.1-128      dplyr_0.5.0      
## [5] magrittr_1.5      readxl_0.1.1      rmarkdown_1.0     knitr_1.14       
## [9] checkpoint_0.3.18
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.7        munsell_0.4.3      mnormt_1.5-4      
##  [4] colorspace_1.2-6   lattice_0.20-34    R6_2.1.3          
##  [7] highr_0.6          stringr_1.1.0      plyr_1.8.4        
## [10] tools_3.3.2        parallel_3.3.2     grid_3.3.2        
## [13] gtable_0.2.0       psych_1.6.9        DBI_0.5-1         
## [16] htmltools_0.3.5    lazyeval_0.2.0     assertthat_0.1    
## [19] digest_0.6.10      tibble_1.2         RColorBrewer_1.1-2
## [22] reshape2_1.4.1     formatR_1.4        tidyr_0.6.0       
## [25] evaluate_0.9       labeling_0.3       stringi_1.1.1     
## [28] scales_0.4.0       foreign_0.8-67
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
##   <chr>   <fctr>   <fctr>  <fctr>          <fctr>     <fctr>   <dbl>
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
##   <chr>   <fctr>   <fctr>      <fctr>          <fctr>           <fctr>
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
D2 <- L2[["data"]] %>% filter(metabolite != "LCAC total") %>% mutate(metabolite = droplevels(metabolite))
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
##      <fctr>           <fctr> <int>
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
##      <fctr>           <fctr> <int>
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
##      <fctr>           <fctr> <int>
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
##         <fctr>           <fctr> <int>
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
M1 <- D1 %>% lme(fixed, data = ., random = random, correlation = NULL, control = ctrl)
M1 %>% plot
```

![plot of chunk lmeDiagnosticAim1](../figures/lmeDiagnosticAim1-1.png)

```r
M1 %>% ranef %>% plot
```

![plot of chunk lmeDiagnosticAim1](../figures/lmeDiagnosticAim1-2.png)

```r
M1 %>% anova %>% kable
```



|                             | numDF| denDF|     F-value|   p-value|
|:----------------------------|-----:|-----:|-----------:|---------:|
|(Intercept)                  |     1|   524| 173.2318729| 0.0000000|
|genotype                     |     1|    38|   6.2464288| 0.0168826|
|activity                     |     1|    38|   0.9103141| 0.3460614|
|metabolite                   |    14|   524|  94.5005176| 0.0000000|
|genotype:activity            |     1|    38|   0.1426504| 0.7077617|
|genotype:metabolite          |    14|   524|   3.3163227| 0.0000409|
|activity:metabolite          |    14|   524|   0.2932042| 0.9945942|
|genotype:activity:metabolite |    14|   524|   0.9079309| 0.5500735|

```r
M1 %>% tidy(effects = "fixed") %>% kable
```



|term                                                 |   estimate| std.error|  statistic|   p.value|
|:----------------------------------------------------|----------:|---------:|----------:|---------:|
|(Intercept)                                          |  0.7797332| 0.2535225|  3.0755970| 0.0022105|
|genotypeKO                                           |  0.1114578| 0.3494567|  0.3189459| 0.7515142|
|activityExercise                                     |  0.2779262| 0.3418497|  0.8130068| 0.4212789|
|metabolitearginine                                   | -0.3531442| 0.3494387| -1.0106040| 0.3126724|
|metaboliteCITRIC                                     |  0.2699043| 0.3585155|  0.7528386| 0.4518850|
|metaboliteFUMARIC                                    | -1.7086321| 0.3585155| -4.7658524| 0.0000024|
|metaboliteglutamine                                  |  1.2561893| 0.3494387|  3.5948773| 0.0003552|
|metaboliteisoleucine                                 | -0.5365811| 0.3494387| -1.5355514| 0.1252520|
|metaboliteLACTIC                                     |  1.6755101| 0.3585155|  4.6734656| 0.0000038|
|metaboliteLCAC total                                 |  0.1172417| 0.3494387|  0.3355145| 0.7373713|
|metaboliteleucine                                    | -0.3703517| 0.3494387| -1.0598472| 0.2897024|
|metaboliteMALIC                                      |  0.1387884| 0.3585155|  0.3871197| 0.6988248|
|metaboliteMCAC Total                                 | -1.2618693| 0.3494387| -3.6111317| 0.0003342|
|metaboliteMETHYLSUCCINIC                             | -2.0968883| 0.3585155| -5.8488072| 0.0000000|
|metabolitePYRUVIC_P2P                                | -0.9377425| 0.3585155| -2.6156259| 0.0091632|
|metaboliteSUCCINIC-2                                 |  0.3515290| 0.3585155|  0.9805126| 0.3272856|
|metabolitevaline                                     | -0.3996651| 0.3494387| -1.1437343| 0.2532561|
|genotypeKO:activityExercise                          | -0.2622578| 0.4767548| -0.5500896| 0.5854771|
|genotypeKO:metabolitearginine                        | -0.0848849| 0.4876346| -0.1740748| 0.8618739|
|genotypeKO:metaboliteCITRIC                          | -0.1808884| 0.4941796| -0.3660379| 0.7144845|
|genotypeKO:metaboliteFUMARIC                         | -0.1459773| 0.4941796| -0.2953932| 0.7678105|
|genotypeKO:metaboliteglutamine                       | -0.0937053| 0.4876346| -0.1921629| 0.8476891|
|genotypeKO:metaboliteisoleucine                      | -0.1554508| 0.4876346| -0.3187854| 0.7500163|
|genotypeKO:metaboliteLACTIC                          | -0.2130032| 0.4941796| -0.4310239| 0.6666282|
|genotypeKO:metaboliteLCAC total                      |  0.2582994| 0.4876346|  0.5296986| 0.5965453|
|genotypeKO:metaboliteleucine                         | -0.1522573| 0.4876346| -0.3122365| 0.7549850|
|genotypeKO:metaboliteMALIC                           | -0.1775990| 0.4941796| -0.3593816| 0.7194544|
|genotypeKO:metaboliteMCAC Total                      | -0.2593079| 0.4876346| -0.5317669| 0.5951128|
|genotypeKO:metaboliteMETHYLSUCCINIC                  | -0.3177490| 0.4941796| -0.6429828| 0.5205165|
|genotypeKO:metabolitePYRUVIC_P2P                     | -2.5513187| 0.4941796| -5.1627361| 0.0000003|
|genotypeKO:metaboliteSUCCINIC-2                      | -0.1502577| 0.4941796| -0.3040549| 0.7612067|
|genotypeKO:metabolitevaline                          | -0.1347692| 0.4876346| -0.2763734| 0.7823703|
|activityExercise:metabolitearginine                  | -0.5178056| 0.4767296| -1.0861621| 0.2779067|
|activityExercise:metaboliteCITRIC                    | -0.2035335| 0.4834223| -0.4210263| 0.6739085|
|activityExercise:metaboliteFUMARIC                   | -0.2704327| 0.4834223| -0.5594130| 0.5761189|
|activityExercise:metaboliteglutamine                 | -0.3379438| 0.4767296| -0.7088794| 0.4787146|
|activityExercise:metaboliteisoleucine                | -0.2837486| 0.4767296| -0.5951981| 0.5519680|
|activityExercise:metaboliteLACTIC                    | -0.2581078| 0.4834223| -0.5339179| 0.5936247|
|activityExercise:metaboliteLCAC total                | -0.4173556| 0.4767296| -0.8754556| 0.3817272|
|activityExercise:metaboliteleucine                   | -0.3163899| 0.4767296| -0.6636674| 0.5071951|
|activityExercise:metaboliteMALIC                     | -0.1439485| 0.4834223| -0.2977698| 0.7659969|
|activityExercise:metaboliteMCAC Total                | -0.3742866| 0.4767296| -0.7851130| 0.4327422|
|activityExercise:metaboliteMETHYLSUCCINIC            | -0.3997036| 0.4834223| -0.8268209| 0.4087149|
|activityExercise:metabolitePYRUVIC_P2P               | -1.3300154| 0.4834223| -2.7512497| 0.0061423|
|activityExercise:metaboliteSUCCINIC-2                | -0.3087569| 0.4834223| -0.6386898| 0.5233038|
|activityExercise:metabolitevaline                    | -0.2504673| 0.4767296| -0.5253864| 0.5995368|
|genotypeKO:activityExercise:metabolitearginine       |  0.3551782| 0.6694138|  0.5305810| 0.5959339|
|genotypeKO:activityExercise:metaboliteCITRIC         |  0.1822831| 0.6741964|  0.2703709| 0.7869814|
|genotypeKO:activityExercise:metaboliteFUMARIC        |  0.2315081| 0.6741964|  0.3433837| 0.7314476|
|genotypeKO:activityExercise:metaboliteglutamine      |  0.2463394| 0.6694138|  0.3679927| 0.7130273|
|genotypeKO:activityExercise:metaboliteisoleucine     |  0.2951461| 0.6694138|  0.4409024| 0.6594655|
|genotypeKO:activityExercise:metaboliteLACTIC         |  0.2101961| 0.6741964|  0.3117728| 0.7553372|
|genotypeKO:activityExercise:metaboliteLCAC total     |  0.3722348| 0.6694138|  0.5560609| 0.5784066|
|genotypeKO:activityExercise:metaboliteleucine        | -0.2996051| 0.6694138| -0.4475633| 0.6546533|
|genotypeKO:activityExercise:metaboliteMALIC          |  0.1829736| 0.6741964|  0.2713950| 0.7861941|
|genotypeKO:activityExercise:metaboliteMCAC Total     |  0.2012857| 0.6694138|  0.3006895| 0.7637705|
|genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC |  0.3903953| 0.6741964|  0.5790527| 0.5628024|
|genotypeKO:activityExercise:metabolitePYRUVIC_P2P    |  1.8254609| 0.6741964|  2.7076099| 0.0069983|
|genotypeKO:activityExercise:metaboliteSUCCINIC-2     |  0.1282931| 0.6741964|  0.1902903| 0.8491553|
|genotypeKO:activityExercise:metabolitevaline         |  0.3640955| 0.6694138|  0.5439021| 0.5867402|

```r
M1 %>% tidy(effects = "fixed") %>% write.csv(file = "../data/processed/lmeFixedCoefAim1.csv", row.names = FALSE)
```

Calculate contrasts of *Exercise vs Rest* given metabolite and genotype.


```r
Ftests <- rbind(contrastGenotype("WT", D1, "activity", "Exercise"),
                contrastGenotype("KO", D1, "activity", "Exercise"))
Ftests %>% kable
```



|contrast |metabolite       |genotype |       beta| numDF| denDF|   F.value|   p.value|
|:--------|:----------------|:--------|----------:|-----:|-----:|---------:|---------:|
|Exercise |3-HYDROXYBUTYRIC |WT       |  0.2779262|     1|    38| 0.6609801| 0.4212789|
|Exercise |arginine         |WT       | -0.2398795|     1|    38| 0.5210553| 0.4748105|
|Exercise |CITRIC           |WT       |  0.0743927|     1|    38| 0.0473576| 0.8288911|
|Exercise |FUMARIC          |WT       |  0.0074935|     1|    38| 0.0004805| 0.9826262|
|Exercise |glutamine        |WT       | -0.0600176|     1|    38| 0.0326178| 0.8576381|
|Exercise |isoleucine       |WT       | -0.0058224|     1|    38| 0.0003070| 0.9861129|
|Exercise |LACTIC           |WT       |  0.0198184|     1|    38| 0.0033610| 0.9540732|
|Exercise |LCAC total       |WT       | -0.1394294|     1|    38| 0.1760380| 0.6771631|
|Exercise |leucine          |WT       | -0.0384637|     1|    38| 0.0133968| 0.9084646|
|Exercise |MALIC            |WT       |  0.1339776|     1|    38| 0.1536010| 0.6973070|
|Exercise |MCAC Total       |WT       | -0.0963604|     1|    38| 0.0840805| 0.7734184|
|Exercise |METHYLSUCCINIC   |WT       | -0.1217775|     1|    38| 0.1269005| 0.7236379|
|Exercise |PYRUVIC_P2P      |WT       | -1.0520892|     1|    38| 9.4718467| 0.0038602|
|Exercise |SUCCINIC-2       |WT       | -0.0308307|     1|    38| 0.0081339| 0.9286116|
|Exercise |valine           |WT       |  0.0274589|     1|    38| 0.0068275| 0.9345804|
|Exercise |3-HYDROXYBUTYRIC |KO       |  0.0156683|     1|    38| 0.0022230| 0.9626415|
|Exercise |arginine         |KO       | -0.1469590|     1|    38| 0.1955646| 0.6608319|
|Exercise |CITRIC           |KO       | -0.0055821|     1|    38| 0.0002822| 0.9866861|
|Exercise |FUMARIC          |KO       | -0.0232563|     1|    38| 0.0048976| 0.9445745|
|Exercise |glutamine        |KO       | -0.0759361|     1|    38| 0.0522148| 0.8204784|
|Exercise |isoleucine       |KO       |  0.0270659|     1|    38| 0.0066335| 0.9355145|
|Exercise |LACTIC           |KO       | -0.0322433|     1|    38| 0.0094141| 0.9232156|
|Exercise |LCAC total       |KO       | -0.0294524|     1|    38| 0.0078549| 0.9298432|
|Exercise |leucine          |KO       | -0.6003266|     1|    38| 3.2634195| 0.0787656|
|Exercise |MALIC            |KO       |  0.0546934|     1|    38| 0.0270874| 0.8701450|
|Exercise |MCAC Total       |KO       | -0.1573325|     1|    38| 0.2241480| 0.6386060|
|Exercise |METHYLSUCCINIC   |KO       |  0.0063600|     1|    38| 0.0003663| 0.9848308|
|Exercise |PYRUVIC_P2P      |KO       |  0.5111139|     1|    38| 2.3655552| 0.1323259|
|Exercise |SUCCINIC-2       |KO       | -0.1647955|     1|    38| 0.2459168| 0.6228211|
|Exercise |valine           |KO       |  0.1292966|     1|    38| 0.1513812| 0.6993907|

```r
Ftests %>% write.csv(file = "../data/processed/contrastsAim1.csv", row.names = FALSE)
```

Click [link to figure]("../figures/plotDataAim1.png").


## Aim 2: Chow

**MCT chow versus experimental chow with triheptanoin**

Estimate model.
Specify the correlation structure using `cs`.
Use `corSymm`, *general correlation matrix, with no additional structure*.


```r
rm(M)
```

```
## Warning in rm(M): object 'M' not found
```

```r
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
M2 <- D2 %>% lme(fixed, data = ., random = random, correlation = NULL, control = ctrl)
M2 %>% plot   
```

![plot of chunk lmeDiagnosticAim2](../figures/lmeDiagnosticAim2-1.png)

```r
M2 %>% ranef %>% plot
```

![plot of chunk lmeDiagnosticAim2](../figures/lmeDiagnosticAim2-2.png)

```r
M2 %>% anova %>% kable
```



|                         | numDF| denDF|     F-value|   p-value|
|:------------------------|-----:|-----:|-----------:|---------:|
|(Intercept)              |     1|   546| 491.0054200| 0.0000000|
|genotype                 |     1|    39|   3.1749786| 0.0825621|
|chow                     |     1|    39|   0.4012764| 0.5301285|
|metabolite               |    14|   546| 117.9870484| 0.0000000|
|genotype:chow            |     1|    39|   0.2250611| 0.6378560|
|genotype:metabolite      |    14|   546|   2.9759057| 0.0002065|
|chow:metabolite          |    14|   546|   1.9324058| 0.0210998|
|genotype:chow:metabolite |    14|   546|   0.5017797| 0.9325101|

```r
M2 %>% tidy(effects = "fixed") %>% kable
```



|term                                                  |   estimate| std.error|  statistic|   p.value|
|:-----------------------------------------------------|----------:|---------:|----------:|---------:|
|(Intercept)                                           |  0.2557123| 0.1704581|  1.5001475| 0.1341542|
|genotypeKO                                            |  0.6626366| 0.2410642|  2.7487976| 0.0090147|
|chowYellow (C8)                                       |  0.7527158| 0.2410642|  3.1224708| 0.0033733|
|metabolitearginine                                    |  0.1236579| 0.2409803|  0.5131452| 0.6080573|
|metaboliteCITRIC                                      |  0.0395968| 0.2409803|  0.1643155| 0.8695436|
|metaboliteFUMARIC                                     | -1.0066586| 0.2409803| -4.1773484| 0.0000343|
|metaboliteglutamine                                   |  1.7844750| 0.2409803|  7.4050661| 0.0000000|
|metaboliteisoleucine                                  | -0.0276526| 0.2409803| -0.1147507| 0.9086849|
|metaboliteLACTIC                                      |  2.1516434| 0.2409803|  8.9287110| 0.0000000|
|metaboliteLC even AC total                            |  0.4067190| 0.2409803|  1.6877689| 0.0920267|
|metaboliteLC odd AC total                             | -0.5647188| 0.2409803| -2.3434234| 0.0194655|
|metaboliteleucine                                     |  0.0958975| 0.2409803|  0.3979475| 0.6908245|
|metaboliteMALIC                                       |  0.7454015| 0.2409803|  3.0932053| 0.0020813|
|metaboliteMCAC total                                  | -0.7176324| 0.2409803| -2.9779712| 0.0030305|
|metaboliteMETHYLSUCCINIC                              | -1.4975941| 0.2409803| -6.2145914| 0.0000000|
|metaboliteSUCCINIC-2                                  |  0.9674087| 0.2409803|  4.0144721| 0.0000679|
|metabolitevaline                                      |  0.1464322| 0.2409803|  0.6076524| 0.5436707|
|genotypeKO:chowYellow (C8)                            | -0.7024715| 0.3451513| -2.0352565| 0.0486612|
|genotypeKO:metabolitearginine                         | -0.6736076| 0.3407976| -1.9765620| 0.0485939|
|genotypeKO:metaboliteCITRIC                           |  0.0088516| 0.3407976|  0.0259733| 0.9792881|
|genotypeKO:metaboliteFUMARIC                          | -0.9692841| 0.3407976| -2.8441635| 0.0046198|
|genotypeKO:metaboliteglutamine                        | -0.7667286| 0.3407976| -2.2498065| 0.0248587|
|genotypeKO:metaboliteisoleucine                       | -0.7247245| 0.3407976| -2.1265539| 0.0339050|
|genotypeKO:metaboliteLACTIC                           | -0.8237407| 0.3407976| -2.4170965| 0.0159715|
|genotypeKO:metaboliteLC even AC total                 | -0.1083540| 0.3407976| -0.3179423| 0.7506502|
|genotypeKO:metaboliteLC odd AC total                  |  0.0839471| 0.3407976|  0.2463254| 0.8055229|
|genotypeKO:metaboliteleucine                          | -0.6996743| 0.3407976| -2.0530493| 0.0405437|
|genotypeKO:metaboliteMALIC                            | -0.7985721| 0.3407976| -2.3432446| 0.0194748|
|genotypeKO:metaboliteMCAC total                       | -0.5409875| 0.3407976| -1.5874158| 0.1129974|
|genotypeKO:metaboliteMETHYLSUCCINIC                   | -0.9184289| 0.3407976| -2.6949395| 0.0072570|
|genotypeKO:metaboliteSUCCINIC-2                       | -0.7749538| 0.3407976| -2.2739416| 0.0233568|
|genotypeKO:metabolitevaline                           | -0.7162648| 0.3407976| -2.1017307| 0.0360344|
|chowYellow (C8):metabolitearginine                    | -0.8936187| 0.3407976| -2.6221391| 0.0089818|
|chowYellow (C8):metaboliteCITRIC                      | -0.6917746| 0.3407976| -2.0298694| 0.0428542|
|chowYellow (C8):metaboliteFUMARIC                     | -0.9370847| 0.3407976| -2.7496810| 0.0061630|
|chowYellow (C8):metaboliteglutamine                   | -0.7969906| 0.3407976| -2.3386039| 0.0197158|
|chowYellow (C8):metaboliteisoleucine                  | -0.7794798| 0.3407976| -2.2872221| 0.0225645|
|chowYellow (C8):metaboliteLACTIC                      | -0.7969796| 0.3407976| -2.3385717| 0.0197175|
|chowYellow (C8):metaboliteLC even AC total            | -0.7582319| 0.3407976| -2.2248746| 0.0264975|
|chowYellow (C8):metaboliteLC odd AC total             | -1.0009945| 0.3407976| -2.9372112| 0.0034515|
|chowYellow (C8):metaboliteleucine                     | -0.7143303| 0.3407976| -2.0960544| 0.0365371|
|chowYellow (C8):metaboliteMALIC                       | -0.8078578| 0.3407976| -2.3704915| 0.0181108|
|chowYellow (C8):metaboliteMCAC total                  | -0.7605770| 0.3407976| -2.2317557| 0.0260361|
|chowYellow (C8):metaboliteMETHYLSUCCINIC              | -0.2792540| 0.3407976| -0.8194129| 0.4129084|
|chowYellow (C8):metaboliteSUCCINIC-2                  | -0.7361425| 0.3407976| -2.1600579| 0.0312027|
|chowYellow (C8):metabolitevaline                      | -0.6012405| 0.3407976| -1.7642157| 0.0782548|
|genotypeKO:chowYellow (C8):metabolitearginine         |  0.7595912| 0.4879479|  1.5567056| 0.1201198|
|genotypeKO:chowYellow (C8):metaboliteCITRIC           |  0.6915345| 0.4879479|  1.4172302| 0.1569860|
|genotypeKO:chowYellow (C8):metaboliteFUMARIC          |  0.9505151| 0.4879479|  1.9479849| 0.0519288|
|genotypeKO:chowYellow (C8):metaboliteglutamine        |  0.8106335| 0.4879479|  1.6613116| 0.0972249|
|genotypeKO:chowYellow (C8):metaboliteisoleucine       |  0.7003375| 0.4879479|  1.4352711| 0.1517825|
|genotypeKO:chowYellow (C8):metaboliteLACTIC           |  0.6222640| 0.4879479|  1.2752674| 0.2027570|
|genotypeKO:chowYellow (C8):metaboliteLC even AC total |  0.7517479| 0.4879479|  1.5406315| 0.1239858|
|genotypeKO:chowYellow (C8):metaboliteLC odd AC total  |  0.2546537| 0.4879479|  0.5218871| 0.6019606|
|genotypeKO:chowYellow (C8):metaboliteleucine          |  0.6677446| 0.4879479|  1.3684752| 0.1717262|
|genotypeKO:chowYellow (C8):metaboliteMALIC            |  0.6230042| 0.4879479|  1.2767842| 0.2022212|
|genotypeKO:chowYellow (C8):metaboliteMCAC total       |  0.7140072| 0.4879479|  1.4632858| 0.1439646|
|genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   |  0.9664959| 0.4879479|  1.9807358| 0.0481223|
|genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       |  0.6372278| 0.4879479|  1.3059342| 0.1921249|
|genotypeKO:chowYellow (C8):metabolitevaline           |  0.7500654| 0.4879479|  1.5371834| 0.1248277|

```r
M2 %>% tidy(effects = "fixed") %>% write.csv(file = "../data/processed/lmeFixedCoefAim2.csv", row.names = FALSE)
```

Calculate contrasts of *Yellow (C8) vs White (C7)* given metabolite and genotype.


```r
Ftests <- rbind(contrastGenotype("WT", D2, "chow", "Yellow (C8)"),
                contrastGenotype("KO", D2, "chow", "Yellow (C8)"))
Ftests %>% kable
```



|contrast    |metabolite       |genotype |       beta| numDF| denDF|   F.value|   p.value|
|:-----------|:----------------|:--------|----------:|-----:|-----:|---------:|---------:|
|Yellow (C8) |3-HYDROXYBUTYRIC |WT       |  0.7527158|     1|    39| 9.7498236| 0.0033733|
|Yellow (C8) |arginine         |WT       | -0.1409029|     1|    39| 0.3416447| 0.5622477|
|Yellow (C8) |CITRIC           |WT       |  0.0609412|     1|    39| 0.0639082| 0.8017500|
|Yellow (C8) |FUMARIC          |WT       | -0.1843689|     1|    39| 0.5849383| 0.4489865|
|Yellow (C8) |glutamine        |WT       | -0.0442748|     1|    39| 0.0337325| 0.8552286|
|Yellow (C8) |isoleucine       |WT       | -0.0267640|     1|    39| 0.0123264| 0.9121664|
|Yellow (C8) |LACTIC           |WT       | -0.0442638|     1|    39| 0.0337157| 0.8552641|
|Yellow (C8) |LC even AC total |WT       | -0.0055162|     1|    39| 0.0005236| 0.9818607|
|Yellow (C8) |LC odd AC total  |WT       | -0.2482787|     1|    39| 1.0607518| 0.3093891|
|Yellow (C8) |leucine          |WT       |  0.0383855|     1|    39| 0.0253553| 0.8743067|
|Yellow (C8) |MALIC            |WT       | -0.0551420|     1|    39| 0.0523239| 0.8202630|
|Yellow (C8) |MCAC total       |WT       | -0.0078612|     1|    39| 0.0010634| 0.9741516|
|Yellow (C8) |METHYLSUCCINIC   |WT       |  0.4734618|     1|    39| 3.8574887| 0.0566804|
|Yellow (C8) |SUCCINIC-2       |WT       |  0.0165732|     1|    39| 0.0047266| 0.9455397|
|Yellow (C8) |valine           |WT       |  0.1514753|     1|    39| 0.3948375| 0.5334304|
|Yellow (C8) |3-HYDROXYBUTYRIC |KO       |  0.0502442|     1|    39| 0.0413731| 0.8398770|
|Yellow (C8) |arginine         |KO       | -0.0837832|     1|    39| 0.1150428| 0.7362941|
|Yellow (C8) |CITRIC           |KO       |  0.0500042|     1|    39| 0.0409787| 0.8406315|
|Yellow (C8) |FUMARIC          |KO       |  0.0636747|     1|    39| 0.0664476| 0.7979365|
|Yellow (C8) |glutamine        |KO       |  0.0638872|     1|    39| 0.0668917| 0.7972775|
|Yellow (C8) |isoleucine       |KO       | -0.0288980|     1|    39| 0.0136862| 0.9074700|
|Yellow (C8) |LACTIC           |KO       | -0.1244713|     1|    39| 0.2539124| 0.6171684|
|Yellow (C8) |LC even AC total |KO       |  0.0437602|     1|    39| 0.0313837| 0.8603039|
|Yellow (C8) |LC odd AC total  |KO       | -0.6960966|     1|    39| 7.9411684| 0.0075494|
|Yellow (C8) |leucine          |KO       |  0.0036586|     1|    39| 0.0002194| 0.9882585|
|Yellow (C8) |MALIC            |KO       | -0.1346094|     1|    39| 0.2969586| 0.5888983|
|Yellow (C8) |MCAC total       |KO       |  0.0036745|     1|    39| 0.0002213| 0.9882073|
|Yellow (C8) |METHYLSUCCINIC   |KO       |  0.7374861|     1|    39| 8.9136005| 0.0048697|
|Yellow (C8) |SUCCINIC-2       |KO       | -0.0486705|     1|    39| 0.0388219| 0.8448252|
|Yellow (C8) |valine           |KO       |  0.1990692|     1|    39| 0.6494616| 0.4251920|

```r
Ftests %>% write.csv(file = "../data/processed/contrastsAim2.csv", row.names = FALSE)
```

Click [link to figure]("../figures/plotDataAim2.png").


Save `lme` objects for interactive work.


```r
save(M1, M2, file = "../data/processed/lmeObjects.RData")
```
