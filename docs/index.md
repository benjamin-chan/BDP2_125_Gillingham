---
title: "Metabolomics of very long-chain aclCoA dehydrogenase knockout mice"
date: "2016-12-29 09:12:53"
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
library(doParallel)
```

```
## Loading required package: foreach
```

```
## Loading required package: iterators
```

```
## Loading required package: parallel
```

Reproducibility steps.


```r
sessionInfo()
```

```
## R version 3.3.2 (2016-10-31)
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
## [1] parallel  stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
##  [1] doParallel_1.0.10 iterators_1.0.8   foreach_1.4.3    
##  [4] ggplot2_2.1.0     broom_0.4.1       nlme_3.1-128     
##  [7] dplyr_0.5.0       magrittr_1.5      readxl_0.1.1     
## [10] rmarkdown_1.0     knitr_1.14        checkpoint_0.3.16
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.7      munsell_0.4.3    mnormt_1.5-4     colorspace_1.2-6
##  [5] lattice_0.20-34  R6_2.1.3         stringr_1.1.0    plyr_1.8.4      
##  [9] tools_3.3.2      grid_3.3.2       gtable_0.2.0     psych_1.6.9     
## [13] DBI_0.5-1        htmltools_0.3.5  assertthat_0.1   digest_0.6.10   
## [17] tibble_1.2       reshape2_1.4.1   formatR_1.4      tidyr_0.6.0     
## [21] codetools_0.2-14 evaluate_0.9     stringi_1.1.1    scales_0.4.0    
## [25] methods_3.3.2    foreign_0.8-67
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

Import Aim 1: Exercise data.


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

* Metabolite `LCAC total` is the sum of `LC even AC total` and `LC odd AC total`; remove `LCAC total` data


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


## Exclusions

> From: Benjamin Chan <chanb@ohsu.edu>  
> Date: Wednesday, December 28, 2016 at 11:30 AM  
> To: Melanie Gillingham <gillingm@ohsu.edu>  
> Cc: Garen Gaston <gastong@ohsu.edu>  
> Subject: RE: Gillingham spreadsheets for stats  
> 
> Hi Garen,
>  
> Can you check these data points? The raw values look suspicious to me, like
> they're the "floor" of the detectable range. I'm wondering if they're real or
> if they should be excluded.
 
```
> D1 %>% filter(value < 1e-6)
# AIM 1
      id genotype activity    chow metabolite_type  metabolite    value  logValue     zValue  zLogValue important
   <chr>   <fctr>   <fctr>  <fctr>          <fctr>      <fctr>    <dbl>     <dbl>      <dbl>      <dbl>     <lgl>
1   1060       KO Exercise Regular     Amino acids     leucine 1.11e-07 -6.954677 -10.830077 -190.02270      TRUE
2   1030       WT Exercise Regular   Organic acids PYRUVIC_P2P 1.11e-07 -6.954677  -1.413972  -13.00909      TRUE
3   1094       WT Exercise Regular   Organic acids PYRUVIC_P2P 1.11e-07 -6.954677  -1.413972  -13.00909      TRUE
4   1134       KO     Rest Regular   Organic acids PYRUVIC_P2P 1.11e-07 -6.954677  -1.413972  -13.00909      TRUE
5 1204B        KO     Rest Regular   Organic acids PYRUVIC_P2P 1.11e-07 -6.954677  -1.413972  -13.00909      TRUE
6 1205B        KO     Rest Regular   Organic acids PYRUVIC_P2P 1.11e-07 -6.954677  -1.413972  -13.00909      TRUE
7   1014       KO Exercise Regular   Organic acids PYRUVIC_P2P 1.11e-07 -6.954677  -1.413972  -13.00909      TRUE
8   1017       KO Exercise Regular   Organic acids PYRUVIC_P2P 1.11e-07 -6.954677  -1.413972  -13.00909      TRUE
9   1019       KO Exercise Regular   Organic acids PYRUVIC_P2P 1.11e-07 -6.954677  -1.413972  -13.00909      TRUE

> D2 %>% filter(value < 1e-6)
# AIM 2
     id genotype activity        chow metabolite_type       metabolite    value  logValue    zValue  zLogValue important
  <chr>   <fctr>   <fctr>      <fctr>          <fctr>           <fctr>    <dbl>     <dbl>     <dbl>      <dbl>     <lgl>
1  1199       WT Exercise  White (C7)   Organic acids 3-HYDROXYBUTYRIC 1.11e-07 -6.954677 -2.516025 -42.518943      TRUE
2  1101       WT Exercise Yellow (C8)   Organic acids           CITRIC 1.11e-07 -6.954677 -2.601795  -3.013728      TRUE
3  1194       WT Exercise  White (C7)   Organic acids           CITRIC 1.11e-07 -6.954677 -2.601795  -3.013728      TRUE
```

### Aim 1

> From: Melanie Gillingham   
> Sent: Wednesday, December 28, 2016 1:15 PM  
> To: Benjamin Chan <chanb@ohsu.edu>  
> Cc: Garen Gaston <gastong@ohsu.edu>  
> Subject: Re: Gillingham spreadsheets for stats  
> 
> Garen and I are going over the min and max for aim 1. We both agree the low
> pyruvic values listed below and the low leucine value are most likely
> undetectable peaks. This could be due to peak misalignment, sample processing,
> or a whole list of other issues and they should probably be deleted.
> 
> The range of log z-scores for the key metabolites in aim 1 is -14 to 9.8 which
> seems to suggest some variability in the data.


```r
D1 <- D1 %>% mutate(exclude = value < 1e-6)
message(sprintf("Excluding %d data points", D1 %>% filter(exclude) %>% tally %>% as.integer))
```

```
## Excluding 9 data points
```

```r
D1 %>% filter(exclude) %>% kable
```



|id    |genotype |activity |chow    |metabolite_type |metabolite  | value|  logValue|     zValue|  zLogValue|important |exclude |
|:-----|:--------|:--------|:-------|:---------------|:-----------|-----:|---------:|----------:|----------:|:---------|:-------|
|1060  |KO       |Exercise |Regular |Amino acids     |leucine     | 1e-07| -6.954677| -10.830077| -190.02270|TRUE      |TRUE    |
|1030  |WT       |Exercise |Regular |Organic acids   |PYRUVIC_P2P | 1e-07| -6.954677|  -1.413972|  -13.00909|TRUE      |TRUE    |
|1094  |WT       |Exercise |Regular |Organic acids   |PYRUVIC_P2P | 1e-07| -6.954677|  -1.413972|  -13.00909|TRUE      |TRUE    |
|1134  |KO       |Rest     |Regular |Organic acids   |PYRUVIC_P2P | 1e-07| -6.954677|  -1.413972|  -13.00909|TRUE      |TRUE    |
|1204B |KO       |Rest     |Regular |Organic acids   |PYRUVIC_P2P | 1e-07| -6.954677|  -1.413972|  -13.00909|TRUE      |TRUE    |
|1205B |KO       |Rest     |Regular |Organic acids   |PYRUVIC_P2P | 1e-07| -6.954677|  -1.413972|  -13.00909|TRUE      |TRUE    |
|1014  |KO       |Exercise |Regular |Organic acids   |PYRUVIC_P2P | 1e-07| -6.954677|  -1.413972|  -13.00909|TRUE      |TRUE    |
|1017  |KO       |Exercise |Regular |Organic acids   |PYRUVIC_P2P | 1e-07| -6.954677|  -1.413972|  -13.00909|TRUE      |TRUE    |
|1019  |KO       |Exercise |Regular |Organic acids   |PYRUVIC_P2P | 1e-07| -6.954677|  -1.413972|  -13.00909|TRUE      |TRUE    |

```r
D1 %>% filter(!exclude) %>% summarize(min(zLogValue), max(zLogValue))
```

```
## # A tibble: 1 × 2
##   `min(zLogValue)` `max(zLogValue)`
##              <dbl>            <dbl>
## 1        -14.42533         9.823728
```

```r
D1 <- D1 %>% filter(!exclude)
```

### Aim 2

> From: Melanie Gillingham   
> Sent: Wednesday, December 28, 2016 11:50 AM  
> To: Benjamin Chan <chanb@ohsu.edu>  
> Cc: Garen Gaston <gastong@ohsu.edu>  
> Subject: Re: Gillingham spreadsheets for stats  
> 
> Garen and I are looking at the min and max right now. Yes, the 1.11e-7 is a
> place holder for below the limit of detection --- or zero. Those values could be
> excluded. We think the 1199 beta hydroxybutyrate is a machine error. The
> fumeric value could also be an error.
> 
> For specific aim 2 the range of the key metabolites for the log z scores looks
> like -9 to 11. That seems like a spread of values around the normal to me.


```r
D2 <- D2 %>% mutate(exclude = value < 1e-6 | zLogValue > 20)
message(sprintf("Excluding %d data points", D2 %>% filter(exclude) %>% tally %>% as.integer))
```

```
## Excluding 4 data points
```

```r
D2 %>% filter(exclude) %>% kable
```



|id   |genotype |activity |chow        |metabolite_type |metabolite       |      value|  logValue|     zValue|  zLogValue|important |exclude |
|:----|:--------|:--------|:-----------|:---------------|:----------------|----------:|---------:|----------:|----------:|:---------|:-------|
|1199 |WT       |Exercise |White (C7)  |Organic acids   |3-HYDROXYBUTYRIC |  0.0000001| -6.954677|  -2.516025| -42.518943|TRUE      |TRUE    |
|1101 |WT       |Exercise |Yellow (C8) |Organic acids   |CITRIC           |  0.0000001| -6.954677|  -2.601795|  -3.013728|TRUE      |TRUE    |
|1194 |WT       |Exercise |White (C7)  |Organic acids   |CITRIC           |  0.0000001| -6.954677|  -2.601795|  -3.013728|TRUE      |TRUE    |
|1199 |WT       |Exercise |White (C7)  |Organic acids   |FUMARIC          | 30.5975143|  1.485686| 774.289007|  20.038743|TRUE      |TRUE    |

```r
D2 %>% filter(!exclude) %>% summarize(min(zLogValue), max(zLogValue))
```

```
## # A tibble: 1 × 2
##   `min(zLogValue)` `max(zLogValue)`
##              <dbl>            <dbl>
## 1        -9.673172         11.97689
```

```r
D2 <- D2 %>% filter(!exclude)
```


## Check data

Check the `value` and `zLogValue`.


```r
kable(summarizeOutcome(D1))
```



|                         |       Min.| X1st.Qu.|  Median|     Mean| X3rd.Qu.|    Max.|
|:------------------------|----------:|--------:|-------:|--------:|--------:|-------:|
|raw                      |   0.006574|   1.2630|  3.1020| 28.84000|  12.3400| 438.700|
|log-transform            |  -2.182000|   0.1014|  0.4917|  0.50970|   1.0910|   2.642|
|normalized raw           |  -6.384000|  -0.9750| -0.3241|  0.01552|   0.6346|  15.100|
|normalized log-transform | -14.430000|  -1.0640| -0.2196| -0.16020|   0.6920|   9.824|

```r
D1 %>%
  ggplot +
  aes(x = zLogValue, y = metabolite, color = activity, fill = activity) +
  geom_jitter(alpha = 1/2) +
  facet_wrap(~ genotype) +
  scale_x_continuous("Normalized log-transformed values") +
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



|                         |     Min.| X1st.Qu.|  Median|     Mean| X3rd.Qu.|    Max.|
|:------------------------|--------:|--------:|-------:|--------:|--------:|-------:|
|raw                      |  0.02303|  0.91940|  3.3490| 24.88000|  12.0700| 408.800|
|log-transform            | -1.63800| -0.03649|  0.5249|  0.52900|   1.0820|   2.612|
|normalized raw           | -6.59600| -0.91660| -0.3149|  0.80960|   0.6589|  56.330|
|normalized log-transform | -9.67300| -0.98980| -0.1182|  0.02061|   0.7084|  11.980|

```r
D2 %>%
  ggplot +
  aes(x = zLogValue, y = metabolite, color = chow, fill = chow) +
  geom_jitter(alpha = 1/2) +
  facet_wrap(~ genotype) +
  scale_x_continuous("Normalized log-transformed values") +
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
## 2       WT Exercise   163
## 3       KO     Rest   147
## 4       KO Exercise   161
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
## 13       WT      PYRUVIC_P2P    18
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
## 24       KO          leucine    20
## 25       KO            MALIC    21
## 26       KO       MCAC Total    21
## 27       KO   METHYLSUCCINIC    21
## 28       KO      PYRUVIC_P2P    15
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
## 13     Rest      PYRUVIC_P2P    16
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
## 24 Exercise          leucine    21
## 25 Exercise            MALIC    22
## 26 Exercise       MCAC Total    22
## 27 Exercise   METHYLSUCCINIC    22
## 28 Exercise      PYRUVIC_P2P    17
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
## 1       WT  White (C7)   162
## 2       WT Yellow (C8)   164
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
## 1        WT 3-HYDROXYBUTYRIC    21
## 2        WT         arginine    22
## 3        WT           CITRIC    20
## 4        WT          FUMARIC    21
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
## 1   White (C7) 3-HYDROXYBUTYRIC    21
## 2   White (C7)         arginine    22
## 3   White (C7)           CITRIC    21
## 4   White (C7)          FUMARIC    21
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
## 18 Yellow (C8)           CITRIC    20
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


## Methods

A mixed linear effects model was estimated for each aim.
Metabolomic values were log-transformed, then normalized.
Fixed effects for Aim 1 were activity, genotype, and metabolite.
Fixed effects for Aim 2 were chow, genotype, and metabolite.
All 2-way and 3-way interactions between fixed effects were included in the models.
Animal ID was the random effect.
A general correlation structure was assumed.
Estimates for the contrasts comparing exercise versus rest for Aim 1, and yellow C8 chow versus white C7 chow for Aim 2, for each combination of genotype and metabolite are presented.
P-values were adjusted to control the false discovery rate, the expected proportion of false discoveries amongst the rejected hypotheses.
The data was analyzed using R version 3.3.2 (2016-10-31) and the `nlme` package version 3.1.128.


## References

Benjamini, Y., and Hochberg, Y.
(1995).
Controlling the false discovery rate: a practical and powerful approach to multiple testing.
*Journal of the Royal Statistical Society Series B* 57, 289â€“300.


```r
citation()
```

```
## 
## To cite R in publications use:
## 
##   R Core Team (2016). R: A language and environment for
##   statistical computing. R Foundation for Statistical Computing,
##   Vienna, Austria. URL https://www.R-project.org/.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {R: A Language and Environment for Statistical Computing},
##     author = {{R Core Team}},
##     organization = {R Foundation for Statistical Computing},
##     address = {Vienna, Austria},
##     year = {2016},
##     url = {https://www.R-project.org/},
##   }
## 
## We have invested a lot of time and effort in creating R, please
## cite it when using it for data analysis. See also
## 'citation("pkgname")' for citing R packages.
```

```r
citation("nlme")
```

```
## 
## Pinheiro J, Bates D, DebRoy S, Sarkar D and R Core Team (2016).
## _nlme: Linear and Nonlinear Mixed Effects Models_. R package
## version 3.1-128, <URL: http://CRAN.R-project.org/package=nlme>.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {{nlme}: Linear and Nonlinear Mixed Effects Models},
##     author = {Jose Pinheiro and Douglas Bates and Saikat DebRoy and Deepayan Sarkar and {R Core Team}},
##     year = {2016},
##     note = {R package version 3.1-128},
##     url = {http://CRAN.R-project.org/package=nlme},
##   }
```


## Aim 1: Exercise

**Rest versus exhaustive exercise**

Estimate model.
Specify the correlation structure using `cs`.
Use `corSymm`, *general correlation matrix, with no additional structure*.


```r
fixed <- formula(zLogValue ~
                   genotype +
                   activity +
                   metabolite +
                   genotype * activity +
                   genotype * metabolite +
                   activity * metabolite +
                   activity * metabolite * genotype)
random <- formula(~ 1 | id)
ctrl <- lmeControl(opt = "optim",
                   maxIter = 500, msMaxIter = 500,
                   tolerance = 1e-6, niterEM = 25, msMaxEval = 200, msTol = 1e-7)
cs <-
  corSymm(form = random, fixed = FALSE) %>%
  Initialize(data = D1)
Dim(cs)
```

```
## $N
## [1] 613
## 
## $M
## [1] 42
## 
## $maxLen
## [1] 15
## 
## $sumLenSq
## [1] 9013
## 
## $len
## groups
##   1170   1171   1172   1201   1202   1209   1210   1211   1212 1208B  
##     15     15     15     15      7     15     15     15     15     15 
##   1028   1029   1030   1033   1034   1046   1090   1091   1092   1094 
##     15     15     14     15     15     15     15     15     15     14 
##   1095   1134   1135   1150   1151   1163   1164 1204B  1205B  1206B  
##     15     14     15     15     15     15     15     14     14     15 
## 1207B    1010   1012   1014   1017   1018   1019   1060   1066   1073 
##     15     15     15     14     14     15     14     14     15     15 
##   1076   1077 
##     15     15 
## 
## $start
##  [1]  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
## [24] 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41
```

```r
M1 <- D1 %>% lme(fixed, data = ., random = random, correlation = cs, control = ctrl)
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
|(Intercept)                  |     1|   515|   26.651214| 0.0000003|
|genotype                     |     1|    38|   17.795900| 0.0001469|
|activity                     |     1|    38|   35.350458| 0.0000007|
|metabolite                   |    14|   515| 2546.828335| 0.0000000|
|genotype:activity            |     1|    38|    1.169781| 0.2862624|
|genotype:metabolite          |    14|   515|    9.538856| 0.0000000|
|activity:metabolite          |    14|   515|    5.681778| 0.0000000|
|genotype:activity:metabolite |    14|   515|    3.519445| 0.0000153|

```r
M1 %>% tidy(effects = "fixed") %>% kable
```



|term                                                 |   estimate| std.error|  statistic|   p.value|
|:----------------------------------------------------|----------:|---------:|----------:|---------:|
|(Intercept)                                          |  0.0227883| 0.5011002|  0.0454765| 0.9637451|
|genotypeKO                                           |  1.4316904| 0.6908930|  2.0722317| 0.0450791|
|activityExercise                                     |  3.6039489| 0.6758451|  5.3325069| 0.0000047|
|metabolitearginine                                   | -0.0227883| 0.6763983| -0.0336907| 0.9731369|
|metaboliteCITRIC                                     |  0.0000000| 0.6933375|  0.0000000| 1.0000000|
|metaboliteFUMARIC                                    |  0.0000000| 0.6933375|  0.0000000| 1.0000000|
|metaboliteglutamine                                  | -0.0227883| 0.6763983| -0.0336907| 0.9731369|
|metaboliteisoleucine                                 | -0.0227883| 0.6763983| -0.0336907| 0.9731369|
|metaboliteLACTIC                                     |  0.0000000| 0.6933375|  0.0000000| 1.0000000|
|metaboliteLCAC total                                 | -0.0227883| 0.6763983| -0.0336907| 0.9731369|
|metaboliteleucine                                    | -0.0227883| 0.6763983| -0.0336907| 0.9731369|
|metaboliteMALIC                                      |  0.0000000| 0.6933375|  0.0000000| 1.0000000|
|metaboliteMCAC Total                                 | -0.0227883| 0.6763983| -0.0336907| 0.9731369|
|metaboliteMETHYLSUCCINIC                             |  0.0000000| 0.6933375|  0.0000000| 1.0000000|
|metabolitePYRUVIC_P2P                                |  0.0000000| 0.6933375|  0.0000000| 1.0000000|
|metaboliteSUCCINIC-2                                 |  0.0000000| 0.6933375|  0.0000000| 1.0000000|
|metabolitevaline                                     | -0.0227883| 0.6763983| -0.0336907| 0.9731369|
|genotypeKO:activityExercise                          | -3.1567736| 0.8955769| -3.5248495| 0.0011231|
|genotypeKO:metabolitearginine                        | -1.2620040| 0.9434828| -1.3376015| 0.1816168|
|genotypeKO:metaboliteCITRIC                          | -2.2469853| 0.9556997| -2.3511415| 0.0190919|
|genotypeKO:metaboliteFUMARIC                         | -1.9350772| 0.9556997| -2.0247753| 0.0434059|
|genotypeKO:metaboliteglutamine                       | -1.2013529| 0.9434828| -1.2733172| 0.2034799|
|genotypeKO:metaboliteisoleucine                      | -2.0554242| 0.9434828| -2.1785496| 0.0298179|
|genotypeKO:metaboliteLACTIC                          | -2.3190401| 0.9556997| -2.4265363| 0.0155859|
|genotypeKO:metaboliteLCAC total                      |  0.7282474| 0.9434828|  0.7718714| 0.4405446|
|genotypeKO:metaboliteleucine                         | -2.4844849| 0.9434828| -2.6333123| 0.0087098|
|genotypeKO:metaboliteMALIC                           | -2.1763951| 0.9556997| -2.2772791| 0.0231787|
|genotypeKO:metaboliteMCAC Total                      | -2.7767584| 0.9434828| -2.9430937| 0.0033961|
|genotypeKO:metaboliteMETHYLSUCCINIC                  | -2.3095282| 0.9556997| -2.4165834| 0.0160136|
|genotypeKO:metabolitePYRUVIC_P2P                     | -2.5314310| 1.0043202| -2.5205416| 0.0120185|
|genotypeKO:metaboliteSUCCINIC-2                      | -2.1014789| 0.9556997| -2.1988903| 0.0283296|
|genotypeKO:metabolitevaline                          | -1.8096226| 0.9434828| -1.9180239| 0.0556606|
|activityExercise:metabolitearginine                  | -5.1357470| 0.9224036| -5.5677870| 0.0000000|
|activityExercise:metaboliteCITRIC                    | -2.7774879| 0.9348960| -2.9709057| 0.0031079|
|activityExercise:metaboliteFUMARIC                   | -3.5223344| 0.9348960| -3.7676215| 0.0001839|
|activityExercise:metaboliteglutamine                 | -4.3826734| 0.9224036| -4.7513618| 0.0000026|
|activityExercise:metaboliteisoleucine                | -3.6864990| 0.9224036| -3.9966223| 0.0000736|
|activityExercise:metaboliteLACTIC                    | -3.4579580| 0.9348960| -3.6987620| 0.0002399|
|activityExercise:metaboliteLCAC total                | -4.4184264| 0.9224036| -4.7901225| 0.0000022|
|activityExercise:metaboliteleucine                   | -4.5964697| 0.9224036| -4.9831435| 0.0000009|
|activityExercise:metaboliteMALIC                     | -2.1642567| 0.9348960| -2.3149704| 0.0210068|
|activityExercise:metaboliteMCAC Total                | -4.4805887| 0.9224036| -4.8575141| 0.0000016|
|activityExercise:metaboliteMETHYLSUCCINIC            | -4.1314810| 0.9348960| -4.4191876| 0.0000121|
|activityExercise:metabolitePYRUVIC_P2P               | -3.2270171| 0.9586197| -3.3663161| 0.0008187|
|activityExercise:metaboliteSUCCINIC-2                | -4.1408326| 0.9348960| -4.4291904| 0.0000116|
|activityExercise:metabolitevaline                    | -3.1587770| 0.9224036| -3.4245062| 0.0006650|
|genotypeKO:activityExercise:metabolitearginine       |  3.6599936| 1.2133089|  3.0165391| 0.0026831|
|genotypeKO:activityExercise:metaboliteCITRIC         |  2.6228001| 1.2246618|  2.1416526| 0.0326900|
|genotypeKO:activityExercise:metaboliteFUMARIC        |  2.6495961| 1.2113833|  2.1872483| 0.0291734|
|genotypeKO:activityExercise:metaboliteglutamine      |  2.8835679| 1.2321285|  2.3403142| 0.0196484|
|genotypeKO:activityExercise:metaboliteisoleucine     |  3.6371315| 1.2285679|  2.9604644| 0.0032134|
|genotypeKO:activityExercise:metaboliteLACTIC         |  2.7345248| 1.2355967|  2.2131208| 0.0273270|
|genotypeKO:activityExercise:metaboliteLCAC total     |  3.8855450| 1.2102018|  3.2106588| 0.0014070|
|genotypeKO:activityExercise:metaboliteleucine        |  6.7437147| 1.2006905|  5.6165305| 0.0000000|
|genotypeKO:activityExercise:metaboliteMALIC          |  2.4028648| 1.2414248|  1.9355701| 0.0534674|
|genotypeKO:activityExercise:metaboliteMCAC Total     |  2.4351299| 1.2471626|  1.9525361| 0.0514160|
|genotypeKO:activityExercise:metaboliteMETHYLSUCCINIC |  3.6381650| 1.2607577|  2.8856971| 0.0040693|
|genotypeKO:activityExercise:metabolitePYRUVIC_P2P    |  3.6853010| 1.3326654|  2.7653611| 0.0058900|
|genotypeKO:activityExercise:metaboliteSUCCINIC-2     |  0.7841181| 1.2318453|  0.6365394| 0.5247077|
|genotypeKO:activityExercise:metabolitevaline         |  5.5745391| 1.2269501|  4.5434113| 0.0000069|

```r
M1 %>% tidy(effects = "fixed") %>% write.csv(file = "../data/processed/lmeFixedCoefAim1.csv", row.names = FALSE)
```

Calculate contrasts of *Exercise vs Rest* given metabolite and genotype.


```r
metabolites <- c("3-HYDROXYBUTYRIC",
                 "arginine",
                 "CITRIC",  # Returns "Singularity in backsolve at level 0, block 1" error
                 "FUMARIC",
                 "glutamine",
                 "isoleucine",
                 "LACTIC",
                 "LCAC total",  # Returns "Singularity in backsolve at level 0, block 1" error
                 "leucine",
                 "MALIC",
                 "MCAC Total",
                 "METHYLSUCCINIC",
                 "PYRUVIC_P2P",
                 "SUCCINIC-2",
                 "valine")
Ftests <- runClusters(D1, metabolites, fixed, "activity", "Exercise", ctrl)
```

```
## Loading required package: data.table
```

```
## -------------------------------------------------------------------------
```

```
## data.table + dplyr code now lives in dtplyr.
## Please library(dtplyr)!
```

```
## -------------------------------------------------------------------------
```

```
## 
## Attaching package: 'data.table'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     between, last
```

```r
Ftests %>% kable
```



|contrast |metabolite       |genotype |       beta|        se| numDF| denDF|    F.value|   p.value|    lowerCL|    upperCL| p.adjustBH|
|:--------|:----------------|:--------|----------:|---------:|-----:|-----:|----------:|---------:|----------:|----------:|----------:|
|Exercise |3-HYDROXYBUTYRIC |WT       |  3.6039489| 0.6758451|     1|    38| 28.4356299| 0.0000047|  2.2793168|  4.9285811|  0.0001404|
|Exercise |arginine         |WT       | -1.5317981| 0.6572672|     1|    38|  5.4314922| 0.0251848| -2.8200181| -0.2435781|  0.1259238|
|Exercise |CITRIC           |WT       |  0.8263977| 0.6759296|     1|    38|  1.4947731| 0.2290064| -0.4984000|  2.1511954|  0.4565051|
|Exercise |FUMARIC          |WT       |  0.0815774| 0.6758389|     1|    38|  0.0145698| 0.9045604| -1.2430425|  1.4061973|  0.9357521|
|Exercise |glutamine        |WT       | -0.7787245| 0.6572845|     1|    38|  1.4036566| 0.2434694| -2.0669784|  0.5095295|  0.4565051|
|Exercise |isoleucine       |WT       | -0.0825501| 0.6572500|     1|    38|  0.0157752| 0.9007114| -1.3707365|  1.2056364|  0.9357521|
|Exercise |LACTIC           |WT       |  0.1459493| 0.6758836|     1|    38|  0.0466294| 0.8301906| -1.1787583|  1.4706569|  0.9351152|
|Exercise |LCAC total       |WT       | -0.8144775| 0.6572852|     1|    38|  1.5355023| 0.2228878| -2.1027329|  0.4737778|  0.4565051|
|Exercise |leucine          |WT       | -0.9925208| 0.6573031|     1|    38|  2.2800682| 0.1393173| -2.2808112|  0.2957697|  0.4179520|
|Exercise |MALIC            |WT       |  1.4396267| 0.6759943|     1|    38|  4.5353810| 0.0397381|  0.1147023|  2.7645511|  0.1703061|
|Exercise |MCAC Total       |WT       | -0.8766398| 0.6572449|     1|    38|  1.7790488| 0.1902059| -2.1648162|  0.4115366|  0.4565051|
|Exercise |METHYLSUCCINIC   |WT       | -0.5275962| 0.6759014|     1|    38|  0.6093077| 0.4398843| -1.8523385|  0.7971461|  0.6850117|
|Exercise |PYRUVIC_P2P      |WT       |  0.3768007| 0.7084068|     1|    38|  0.2829164| 0.5978922| -1.0116510|  1.7652525|  0.7798594|
|Exercise |SUCCINIC-2       |WT       | -0.5369476| 0.6759606|     1|    38|  0.6309881| 0.4319262| -1.8618060|  0.7879107|  0.6850117|
|Exercise |valine           |WT       |  0.4451719| 0.6573105|     1|    38|  0.4586849| 0.5023419| -0.8431329|  1.7334767|  0.6850117|
|Exercise |3-HYDROXYBUTYRIC |KO       |  0.4515957| 0.5521287|     1|    38|  0.6689888| 0.4185037| -0.6305567|  1.5337480|  0.6850117|
|Exercise |arginine         |KO       | -1.0602767| 0.5394379|     1|    38|  3.8632739| 0.0566934| -2.1175556| -0.0029978|  0.2126003|
|Exercise |CITRIC           |KO       |  0.2281115| 0.6316464|     1|    38|  0.1304206| 0.7199964| -1.0098927|  1.4661157|  0.8639957|
|Exercise |FUMARIC          |KO       | -0.4320436| 0.6143010|     1|    38|  0.4946438| 0.4861502| -1.6360514|  0.7719643|  0.6850117|
|Exercise |glutamine        |KO       | -1.0679252| 0.5962023|     1|    38|  3.2084435| 0.0812268| -2.2364602|  0.1006098|  0.2707561|
|Exercise |isoleucine       |KO       |  0.3907228| 0.5768242|     1|    38|  0.4588287| 0.5022753| -0.7398318|  1.5212774|  0.6850117|
|Exercise |LACTIC           |KO       | -0.2369854| 0.6081970|     1|    38|  0.1518289| 0.6989690| -1.4290295|  0.9550588|  0.8639957|
|Exercise |LCAC total       |KO       | -0.1202451| 0.5975952|     1|    38|  0.0404875| 0.8416037| -1.2915102|  1.0510201|  0.9351152|
|Exercise |leucine          |KO       |  2.5656006| 0.5927838|     1|    38| 18.7320553| 0.0001055|  1.4037657|  3.7274355|  0.0007910|
|Exercise |MALIC            |KO       |  0.7272148| 0.5845137|     1|    38|  1.5478755| 0.2210691| -0.4184110|  1.8728407|  0.4565051|
|Exercise |MCAC Total       |KO       | -1.5975113| 0.6044981|     1|    38|  6.9838986| 0.0118800| -2.7823058| -0.4127167|  0.0712800|
|Exercise |METHYLSUCCINIC   |KO       | -0.0099667| 0.6217004|     1|    38|  0.0002570| 0.9872933| -1.2284770|  1.2085436|  0.9872933|
|Exercise |PYRUVIC_P2P      |KO       |  0.9473895| 0.7195225|     1|    38|  1.7336777| 0.1958281| -0.4628487|  2.3576276|  0.4565051|
|Exercise |SUCCINIC-2       |KO       | -2.8280667| 0.5799819|     1|    38| 23.7766346| 0.0000195| -3.9648103| -1.6913230|  0.0002930|
|Exercise |valine           |KO       |  2.8924597| 0.6257376|     1|    38| 21.3673273| 0.0000429|  1.6660367|  4.1188828|  0.0004285|

```r
Ftests %>% write.csv(file = "../data/processed/contrastsAim1.csv", row.names = FALSE)
```

Click [link to figure](../figures/plotDataAim1.png).


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
fixed <- formula(zLogValue ~
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
## [1] 641
## 
## $M
## [1] 43
## 
## $maxLen
## [1] 15
## 
## $sumLenSq
## [1] 9561
## 
## $len
## groups
##       1101       1102       1103       1184       1185       1186 
##         14         15         15         15         15         15 
##       1195       1196       1197       1203 1192/1198B       1176 
##         15         15         15         15         15         15 
##       1177       1179       1180       1190       1191       1193 
##         15         15         15         15         15         15 
##       1194       1199       1200 1192A/1198       1107       1113 
##         14         13         15         15         15         15 
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
M2 <- D2 %>% lme(fixed, data = ., random = random, correlation = cs, control = ctrl)
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
|(Intercept)              |     1|   542|   0.0475051| 0.8275449|
|genotype                 |     1|    39| 156.1233007| 0.0000000|
|chow                     |     1|    39|   2.8792447| 0.0976967|
|metabolite               |    14|   542| 131.1492272| 0.0000000|
|genotype:chow            |     1|    39|  13.0119134| 0.0008684|
|genotype:metabolite      |    14|   542| 151.7395112| 0.0000000|
|chow:metabolite          |    14|   542|  39.1729742| 0.0000000|
|genotype:chow:metabolite |    14|   542|   8.7331225| 0.0000000|

```r
M2 %>% tidy(effects = "fixed") %>% kable
```



|term                                                  |   estimate| std.error|  statistic|   p.value|
|:-----------------------------------------------------|----------:|---------:|----------:|---------:|
|(Intercept)                                           | -0.2243894| 0.3226303| -0.6955001| 0.4870400|
|genotypeKO                                            | -0.2565878| 0.4476548| -0.5731823| 0.5698117|
|chowYellow (C8)                                       |  0.5318680| 0.4082139|  1.3029149| 0.2002446|
|metabolitearginine                                    |  0.8608379| 0.3929874|  2.1904976| 0.0289134|
|metaboliteCITRIC                                      |  0.4778077| 0.4028581|  1.1860444| 0.2361244|
|metaboliteFUMARIC                                     | -0.1561029| 0.4016448| -0.3886590| 0.6976811|
|metaboliteglutamine                                   |  0.7780084| 0.3929874|  1.9797287| 0.0482394|
|metaboliteisoleucine                                  |  0.4985815| 0.3929874|  1.2686959| 0.2050942|
|metaboliteLACTIC                                      |  0.7313983| 0.3929874|  1.8611242| 0.0632678|
|metaboliteLC even AC total                            |  0.2750335| 0.3929874|  0.6998532| 0.4843193|
|metaboliteLC odd AC total                             |  2.7660864| 0.3929874|  7.0386139| 0.0000000|
|metaboliteleucine                                     | -0.1180236| 0.3929874| -0.3003241| 0.7640451|
|metaboliteMALIC                                       |  1.3214848| 0.3929874|  3.3626646| 0.0008264|
|metaboliteMCAC total                                  |  0.3089526| 0.3929874|  0.7861642| 0.4321148|
|metaboliteMETHYLSUCCINIC                              | -1.8490995| 0.3929874| -4.7052389| 0.0000032|
|metaboliteSUCCINIC-2                                  | -0.2194593| 0.3929874| -0.5584386| 0.5767755|
|metabolitevaline                                      | -1.2555356| 0.3929874| -3.1948497| 0.0014804|
|genotypeKO:chowYellow (C8)                            | -0.3101297| 0.5968976| -0.5195694| 0.6063003|
|genotypeKO:metabolitearginine                         |  0.2070325| 0.5487188|  0.3773017| 0.7060971|
|genotypeKO:metaboliteCITRIC                           |  0.2548507| 0.5558308|  0.4585041| 0.6467741|
|genotypeKO:metaboliteFUMARIC                          | -0.3750257| 0.5549521| -0.6757803| 0.4994684|
|genotypeKO:metaboliteglutamine                        | -1.0449949| 0.5487188| -1.9044269| 0.0573838|
|genotypeKO:metaboliteisoleucine                       | -0.3794905| 0.5487188| -0.6915938| 0.4894886|
|genotypeKO:metaboliteLACTIC                           | -1.5887387| 0.5487188| -2.8953604| 0.0039400|
|genotypeKO:metaboliteLC even AC total                 |  5.3454817| 0.5487188|  9.7417506| 0.0000000|
|genotypeKO:metaboliteLC odd AC total                  |  7.8995681| 0.5487188| 14.3963869| 0.0000000|
|genotypeKO:metaboliteleucine                          | -0.0738024| 0.5487188| -0.1344995| 0.8930575|
|genotypeKO:metaboliteMALIC                            | -2.4479605| 0.5487188| -4.4612295| 0.0000099|
|genotypeKO:metaboliteMCAC total                       |  1.5651727| 0.5487188|  2.8524131| 0.0045045|
|genotypeKO:metaboliteMETHYLSUCCINIC                   | -0.8636346| 0.5487188| -1.5739110| 0.1160915|
|genotypeKO:metaboliteSUCCINIC-2                       | -2.7513876| 0.5487188| -5.0142033| 0.0000007|
|genotypeKO:metabolitevaline                           | -0.2673637| 0.5487188| -0.4872508| 0.6262776|
|chowYellow (C8):metabolitearginine                    | -1.1087712| 0.4531028| -2.4470634| 0.0147190|
|chowYellow (C8):metaboliteCITRIC                      | -0.4813775| 0.5211884| -0.9236151| 0.3560977|
|chowYellow (C8):metaboliteFUMARIC                     |  0.0418750| 0.4756038|  0.0880460| 0.9298726|
|chowYellow (C8):metaboliteglutamine                   | -0.8844583| 0.4599115| -1.9231052| 0.0549904|
|chowYellow (C8):metaboliteisoleucine                  | -0.6455048| 0.4762594| -1.3553640| 0.1758662|
|chowYellow (C8):metaboliteLACTIC                      | -1.0179170| 0.4526461| -2.2488142| 0.0249251|
|chowYellow (C8):metaboliteLC even AC total            | -0.6519009| 0.4636841| -1.4059160| 0.1603221|
|chowYellow (C8):metaboliteLC odd AC total             | -3.0582665| 0.4869092| -6.2809787| 0.0000000|
|chowYellow (C8):metaboliteleucine                     | -0.0309266| 0.4953740| -0.0624307| 0.9502428|
|chowYellow (C8):metaboliteMALIC                       | -1.3289661| 0.4649127| -2.8585280| 0.0044199|
|chowYellow (C8):metaboliteMCAC total                  | -0.5570152| 0.5106332| -1.0908325| 0.2758313|
|chowYellow (C8):metaboliteMETHYLSUCCINIC              |  1.4054188| 0.4541898|  3.0943427| 0.0020743|
|chowYellow (C8):metaboliteSUCCINIC-2                  | -0.1554527| 0.4659413| -0.3336315| 0.7387867|
|chowYellow (C8):metabolitevaline                      |  1.1449130| 0.4675004|  2.4490098| 0.0146406|
|genotypeKO:chowYellow (C8):metabolitearginine         |  0.4276371| 0.7092275|  0.6029618| 0.5467863|
|genotypeKO:chowYellow (C8):metaboliteCITRIC           |  0.3040930| 0.7591676|  0.4005611| 0.6889012|
|genotypeKO:chowYellow (C8):metaboliteFUMARIC          |  0.1930033| 0.7183640|  0.2686706| 0.7882854|
|genotypeKO:chowYellow (C8):metaboliteglutamine        |  1.4090557| 0.7142934|  1.9726568| 0.0490425|
|genotypeKO:chowYellow (C8):metaboliteisoleucine       |  0.0612917| 0.7226224|  0.0848184| 0.9324371|
|genotypeKO:chowYellow (C8):metaboliteLACTIC           | -0.6660769| 0.7156668| -0.9307081| 0.3524189|
|genotypeKO:chowYellow (C8):metaboliteLC even AC total |  0.8510769| 0.7055637|  1.2062369| 0.2282524|
|genotypeKO:chowYellow (C8):metaboliteLC odd AC total  | -4.2357780| 0.7372659| -5.7452517| 0.0000000|
|genotypeKO:chowYellow (C8):metaboliteleucine          | -0.2091283| 0.7307366| -0.2861883| 0.7748432|
|genotypeKO:chowYellow (C8):metaboliteMALIC            | -1.6026369| 0.7216110| -2.2209153| 0.0267693|
|genotypeKO:chowYellow (C8):metaboliteMCAC total       |  0.5885308| 0.7508774|  0.7837908| 0.4335052|
|genotypeKO:chowYellow (C8):metaboliteMETHYLSUCCINIC   |  1.5785435| 0.7116464|  2.2181570| 0.0269578|
|genotypeKO:chowYellow (C8):metaboliteSUCCINIC-2       | -1.3937944| 0.7227139| -1.9285561| 0.0543079|
|genotypeKO:chowYellow (C8):metabolitevaline           |  0.5722234| 0.7175979|  0.7974150| 0.4255592|

```r
M2 %>% tidy(effects = "fixed") %>% write.csv(file = "../data/processed/lmeFixedCoefAim2.csv", row.names = FALSE)
```

Calculate contrasts of *Yellow (C8) vs White (C7)* given metabolite and genotype.


```r
metabolites <- c("3-HYDROXYBUTYRIC",
                 "arginine",
                 "CITRIC",
                 "FUMARIC",
                 "glutamine",
                 "isoleucine",
                 "LACTIC",
                 "LC even AC total",
                 "LC odd AC total",
                 "leucine",
                 "MALIC",
                 "MCAC total",
                 "METHYLSUCCINIC",
                 "SUCCINIC-2",
                 "valine")
Ftests <- runClusters(D2, metabolites, fixed, "chow", "Yellow (C8)", ctrl)
Ftests %>% kable
```



|contrast    |metabolite       |genotype |       beta|        se| numDF| denDF|     F.value|   p.value|    lowerCL|    upperCL| p.adjustBH|
|:-----------|:----------------|:--------|----------:|---------:|-----:|-----:|-----------:|---------:|----------:|----------:|----------:|
|Yellow (C8) |3-HYDROXYBUTYRIC |WT       |  0.5318680| 0.4082139|     1|    39|   1.6975872| 0.2002446| -0.2682166|  1.3319526|  0.4151485|
|Yellow (C8) |arginine         |WT       | -0.6096388| 0.3842302|     1|    39|   2.5174583| 0.1206677| -1.3627161|  0.1434385|  0.3016694|
|Yellow (C8) |CITRIC           |WT       |  0.0888875| 0.4206346|     1|    39|   0.0446551| 0.8337399| -0.7355413|  0.9133162|  0.9263776|
|Yellow (C8) |FUMARIC          |WT       |  0.4970684| 0.3878706|     1|    39|   1.6423227| 0.2075743| -0.2631441|  1.2572808|  0.4151485|
|Yellow (C8) |glutamine        |WT       | -0.3535584| 0.3729375|     1|    39|   0.8987733| 0.3489486| -1.0845024|  0.3773857|  0.4984980|
|Yellow (C8) |isoleucine       |WT       | -0.0842644| 0.3935005|     1|    39|   0.0458562| 0.8315531| -0.8555112|  0.6869823|  0.9263776|
|Yellow (C8) |LACTIC           |WT       | -0.5752375| 0.3815658|     1|    39|   2.2727693| 0.1397230| -1.3230928|  0.1726178|  0.3224376|
|Yellow (C8) |LC even AC total |WT       | -0.1003251| 0.4008493|     1|    39|   0.0626407| 0.8036836| -0.8859752|  0.6853251|  0.9263776|
|Yellow (C8) |LC odd AC total  |WT       | -2.5381786| 0.3904299|     1|    39|  42.2627855| 0.0000001| -3.3034072| -1.7729499|  0.0000010|
|Yellow (C8) |leucine          |WT       |  0.4741568| 0.3832777|     1|    39|   1.5304418| 0.2234436| -0.2770537|  1.2253673|  0.4189568|
|Yellow (C8) |MALIC            |WT       | -0.8999270| 0.3834465|     1|    39|   5.5081403| 0.0241002| -1.6514684| -0.1483857|  0.0723006|
|Yellow (C8) |MCAC total       |WT       | -0.0088246| 0.3865260|     1|    39|   0.0005212| 0.9819019| -0.7664017|  0.7487525|  0.9838357|
|Yellow (C8) |METHYLSUCCINIC   |WT       |  2.0535162| 0.3754463|     1|    39|  29.9157923| 0.0000028|  1.3176549|  2.7893775|  0.0000169|
|Yellow (C8) |SUCCINIC-2       |WT       |  0.4175039| 0.3805110|     1|    39|   1.2038895| 0.2792743| -0.3282840|  1.1632917|  0.4654572|
|Yellow (C8) |valine           |WT       |  1.6317464| 0.3825581|     1|    39|  18.1932639| 0.0001229|  0.8819464|  2.3815464|  0.0005266|
|Yellow (C8) |3-HYDROXYBUTYRIC |KO       |  0.2225476| 0.4417253|     1|    39|   0.2538290| 0.6172260| -0.6432181|  1.0883133|  0.7715325|
|Yellow (C8) |arginine         |KO       | -0.4838852| 0.4298054|     1|    39|   1.2674797| 0.2671205| -1.3262883|  0.3585179|  0.4654572|
|Yellow (C8) |CITRIC           |KO       |  0.0594225| 0.4444267|     1|    39|   0.0178773| 0.8943228| -0.8116379|  0.9304828|  0.9582030|
|Yellow (C8) |FUMARIC          |KO       |  0.4678990| 0.4435791|     1|    39|   1.1126588| 0.2979958| -0.4015001|  1.3372981|  0.4705197|
|Yellow (C8) |glutamine        |KO       |  0.7479619| 0.4282112|     1|    39|   3.0510052| 0.0885595| -0.0913166|  1.5872403|  0.2415260|
|Yellow (C8) |isoleucine       |KO       | -0.3639556| 0.4421544|     1|    39|   0.6775619| 0.4154327| -1.2305624|  0.5026511|  0.5664991|
|Yellow (C8) |LACTIC           |KO       | -1.4623440| 0.4419777|     1|    39|  10.9470808| 0.0020232| -2.3286045| -0.5960836|  0.0075871|
|Yellow (C8) |LC even AC total |KO       |  0.4270631| 0.4208155|     1|    39|   1.0299132| 0.3164339| -0.3977201|  1.2518462|  0.4746509|
|Yellow (C8) |LC odd AC total  |KO       | -7.0721373| 0.4459037|     1|    39| 251.5470626| 0.0000000| -7.9460924| -6.1981821|  0.0000000|
|Yellow (C8) |leucine          |KO       | -0.0088549| 0.4342615|     1|    39|   0.0004158| 0.9838357| -0.8599918|  0.8422820|  0.9838357|
|Yellow (C8) |MALIC            |KO       | -2.7071369| 0.4457004|     1|    39|  36.8921826| 0.0000004| -3.5806937| -1.8335802|  0.0000031|
|Yellow (C8) |MCAC total       |KO       |  0.2688824| 0.4379490|     1|    39|   0.3769449| 0.5428078| -0.5894818|  1.1272466|  0.7080102|
|Yellow (C8) |METHYLSUCCINIC   |KO       |  3.2061670| 0.4463944|     1|    39|  51.5863513| 0.0000000|  2.3312501|  4.0810839|  0.0000002|
|Yellow (C8) |SUCCINIC-2       |KO       | -1.3216006| 0.4461850|     1|    39|   8.7734506| 0.0051833| -2.1961072| -0.4470940|  0.0172777|
|Yellow (C8) |valine           |KO       |  1.9253432| 0.4398920|     1|    39|  19.1568567| 0.0000874|  1.0631708|  2.7875156|  0.0004368|

```r
Ftests %>% write.csv(file = "../data/processed/contrastsAim2.csv", row.names = FALSE)
```

Click [link to figure](../figures/plotDataAim2.png).


## Save

Save `lme` objects for interactive work.


```r
save(M1, M2, file = "../data/processed/lmeObjects.RData")
```
