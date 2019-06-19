---
title: "Metabolomics of very long-chain aclCoA dehydrogenase knockout mice"
date: "2019-06-18 23:31:40"
author: Benjamin Chan (chanb@ohsu.edu)
output:
  html_document:
    toc: true
    theme: simplex
---

---

# Preamble

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
```

```
## Registered S3 methods overwritten by 'ggplot2':
##   method         from 
##   [.quosures     rlang
##   c.quosures     rlang
##   print.quosures rlang
```

```r
library(svglite)
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
## R version 3.6.0 (2019-04-26)
## Platform: x86_64-w64-mingw32/x64 (64-bit)
## Running under: Windows 10 x64 (build 17134)
## 
## Matrix products: default
## 
## locale:
## [1] LC_COLLATE=English_United States.1252 
## [2] LC_CTYPE=English_United States.1252   
## [3] LC_MONETARY=English_United States.1252
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.1252    
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] doParallel_1.0.14 iterators_1.0.10  foreach_1.4.4    
##  [4] svglite_1.2.1     ggplot2_3.1.0     broom_0.5.1      
##  [7] nlme_3.1-137      dplyr_0.8.0.1     magrittr_1.5     
## [10] readxl_1.3.1      rmarkdown_1.12    knitr_1.22       
## [13] checkpoint_0.4.5 
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.1       cellranger_1.1.0 pillar_1.3.1     compiler_3.6.0  
##  [5] plyr_1.8.4       tools_3.6.0      digest_0.6.18    evaluate_0.13   
##  [9] tibble_2.1.1     gtable_0.3.0     lattice_0.20-38  pkgconfig_2.0.2 
## [13] rlang_0.3.3      xfun_0.5         withr_2.1.2      stringr_1.4.0   
## [17] gdtools_0.1.7    generics_0.0.2   grid_3.6.0       tidyselect_0.2.5
## [21] glue_1.3.1       R6_2.4.0         purrr_0.3.2      tidyr_0.8.3     
## [25] codetools_0.2-16 backports_1.1.3  scales_1.0.0     htmltools_0.3.6 
## [29] assertthat_0.2.1 colorspace_1.4-1 stringi_1.4.3    lazyeval_0.2.2  
## [33] munsell_0.5.0    crayon_1.3.4
```

```r
set.seed(as.integer(as.Date("2016-11-18")))
```

Source user-defined functions.


```r
sapply(list.files("lib", full.names = TRUE), source)
```

```
##         lib/library.R
## value   ?            
## visible FALSE
```

---

# Read data

Import the data.
Data files are locally stored.

## Aim 1-a

Import Aim 1 Amino Acids and TCA Cycle intermediates (Neg).


```r
f <- "data/raw/Z scores Aim 1 Neg-aa.xlsx"
L1a <- importDataToList(f)
L1a[["file"]]
```

```
## [1] "data/raw/Z scores Aim 1 Neg-aa.xlsx"
```

```r
L1a[["dim"]]
```

```
## [1] 984   5
```

```r
L1a[["names"]]
```

```
## [1] "condition"  "id"         "metabolite" "z_value"    "genotype"
```

```r
L1a[["head"]]
```

```
## # A tibble: 6 x 5
##   condition    id metabolite z_value genotype
##   <fct>     <dbl> <fct>        <dbl> <fct>   
## 1 Exercise   1010 FUMARIC     -0.934 KO      
## 2 Exercise   1012 FUMARIC     -1.21  KO      
## 3 Exercise   1014 FUMARIC      0.196 KO      
## 4 Exercise   1017 FUMARIC      0.455 KO      
## 5 Exercise   1018 FUMARIC     -0.580 KO      
## 6 Exercise   1019 FUMARIC     -0.593 KO
```

```r
L1a[["data"]] %>% pull(condition) %>% levels()
```

```
## [1] "Rest"     "Exercise"
```

```r
L1a[["data"]] %>% pull(genotype) %>% levels()
```

```
## [1] "WT" "KO"
```

```r
L1a[["data"]] %>% pull(metabolite) %>% levels()
```

```
##  [1] "alaninesarcosine" "arginine"         "asparagine"      
##  [4] "AsparticAcid"     "CITRIC"           "FUMARIC"         
##  [7] "glutamic"         "glutamine"        "glycine"         
## [10] "histidine"        "isoleucine"       "LACTIC"          
## [13] "leucine"          "lysine"           "MALIC"           
## [16] "methionine"       "METHYLSUCCINIC"   "Phenylalanine"   
## [19] "serine"           "SUCCINIC2"        "threonine"       
## [22] "Tryptophan"       "Tyrosine"         "valine"
```

```r
D1a <- L1a[["data"]]
```

## Aim 1-b

Import Aim 1 Acylcarnitines.


```r
f <- "data/raw/Z scores Aim 1 AC.xlsx"
L1b <- importDataToList(f)
L1b[["file"]]
```

```
## [1] "data/raw/Z scores Aim 1 AC.xlsx"
```

```r
L1b[["dim"]]
```

```
## [1] 899   5
```

```r
L1b[["names"]]
```

```
## [1] "condition"  "id"         "metabolite" "z_value"    "genotype"
```

```r
L1b[["head"]]
```

```
## # A tibble: 6 x 5
##   condition    id metabolite      z_value genotype
##   <fct>     <dbl> <fct>             <dbl> <fct>   
## 1 Exercise   1010 acetylcarnitine   -1.84 KO      
## 2 Exercise   1012 acetylcarnitine   -2.15 KO      
## 3 Exercise   1014 acetylcarnitine   -2.04 KO      
## 4 Exercise   1017 acetylcarnitine   -2.43 KO      
## 5 Exercise   1018 acetylcarnitine   -2.49 KO      
## 6 Exercise   1019 acetylcarnitine   -1.83 KO
```

```r
L1b[["data"]] %>% pull(condition) %>% levels()
```

```
## [1] "Rest"     "Exercise"
```

```r
L1b[["data"]] %>% pull(genotype) %>% levels()
```

```
## [1] "WT" "KO"
```

```r
L1b[["data"]] %>% pull(metabolite) %>% levels()
```

```
##  [1] "2methylbutyrylcarnitine" "3HMG"                   
##  [3] "acetylcarnitine"         "butyrylcarnitine"       
##  [5] "C101total"               "C12"                    
##  [7] "C141total"               "C142total"              
##  [9] "C151total"               "C161total"              
## [11] "C171total"               "C181total"              
## [13] "C182total"               "C191total"              
## [15] "carnitine"               "ethylmalonylcarnitine"  
## [17] "isobutyrylcarnitine"     "methylsuccinylcarnitine"
## [19] "ndecanoylcarnitine"      "nhexanoylcarnitine"     
## [21] "noctanoylcarnitine"      "propionylcarnitine"
```

```r
D1b <- L1b[["data"]]
```

## Aim 2-a

Import Aim 2 Amino Acids and TCA Cycle intermediates (Neg).


```r
f <- "data/raw/Z scores Aim 2 Neg-aa.xlsx"
L2a <- importDataToList(f)
L2a[["file"]]
```

```
## [1] "data/raw/Z scores Aim 2 Neg-aa.xlsx"
```

```r
L2a[["dim"]]
```

```
## [1] 1035    5
```

```r
L2a[["names"]]
```

```
## [1] "condition"  "id"         "metabolite" "z_value"    "genotype"
```

```r
L2a[["head"]]
```

```
## # A tibble: 6 x 5
##   condition    id metabolite z_value genotype
##   <fct>     <dbl> <fct>        <dbl> <fct>   
## 1 C7         1120 FUMARIC      -3.23 KO      
## 2 C7         1126 FUMARIC      -3.22 KO      
## 3 C7         1127 FUMARIC      -3.49 KO      
## 4 C7         1128 FUMARIC      -2.31 KO      
## 5 C7         1142 FUMARIC      -1.73 KO      
## 6 C7         1144 FUMARIC      -1.71 KO
```

```r
L2a[["data"]] %>% pull(condition) %>% levels()
```

```
## [1] "C7"       "C8"       "Exercise"
```

```r
L2a[["data"]] %>% pull(genotype) %>% levels()
```

```
## [1] "WT" "KO"
```

```r
L2a[["data"]] %>% pull(metabolite) %>% levels()
```

```
##  [1] "alaninesarcosine" "arginine"         "asparagine"      
##  [4] "AsparticAcid"     "CITRIC"           "FUMARIC"         
##  [7] "glutamic"         "glutamine"        "glycine"         
## [10] "histidine"        "isoleucine"       "LACTIC"          
## [13] "leucine"          "lysine"           "MALIC"           
## [16] "methionine"       "METHYLSUCCINIC"   "Phenylalanine"   
## [19] "serine"           "SUCCINIC2"        "threonine"       
## [22] "Tryptophan"       "Tyrosine"         "valine"
```

```r
D2a <- L2a[["data"]]
```

## Aim 2-b

Import Aim 2 Acylcarnitines.


```r
f <- "data/raw/Z scores Aim 2 AC.xlsx"
L2b <- importDataToList(f)
L2b[["file"]]
```

```
## [1] "data/raw/Z scores Aim 2 AC.xlsx"
```

```r
L2b[["dim"]]
```

```
## [1] 948   5
```

```r
L2b[["names"]]
```

```
## [1] "condition"  "id"         "metabolite" "z_value"    "genotype"
```

```r
L2b[["head"]]
```

```
## # A tibble: 6 x 5
##   condition    id metabolite      z_value genotype
##   <fct>     <dbl> <fct>             <dbl> <fct>   
## 1 C7         1120 acetylcarnitine   -2.57 KO      
## 2 C7         1126 acetylcarnitine   -2.43 KO      
## 3 C7         1127 acetylcarnitine   -2.58 KO      
## 4 C7         1128 acetylcarnitine   -2.11 KO      
## 5 C7         1142 acetylcarnitine   -2.67 KO      
## 6 C7         1143 acetylcarnitine   -2.55 KO
```

```r
L2b[["data"]] %>% pull(condition) %>% levels()
```

```
## [1] "C7"       "C8"       "Exercise"
```

```r
L2b[["data"]] %>% pull(genotype) %>% levels()
```

```
## [1] "WT" "KO"
```

```r
L2b[["data"]] %>% pull(metabolite) %>% levels()
```

```
##  [1] "2methylbutyrylcarnitine" "3HMG"                   
##  [3] "acetylcarnitine"         "butyrylcarnitine"       
##  [5] "C101total"               "C12"                    
##  [7] "C141total"               "C142total"              
##  [9] "C151total"               "C161total"              
## [11] "C171total"               "C181total"              
## [13] "C182total"               "C191total"              
## [15] "carnitine"               "ethylmalonylcarnitine"  
## [17] "isobutyrylcarnitine"     "methylsuccinylcarnitine"
## [19] "ndecanoylcarnitine"      "nhexanoylcarnitine"     
## [21] "noctanoylcarnitine"      "propionylcarnitine"
```

```r
D2b <- L2b[["data"]]
```

---

# Model

Basic data preprocessing steps:

1. Find and remove outliers using Grubbs.
2. Standardize by Z score.
3. Check for normality.
4. Transform and recheck normality. We decided not to use transformed data as it did not correct the few conditions with multiple skewed groups.


## Methods

A mixed linear effects model was estimated for each aim.
Fixed effects for Aim 1 were condition (wildtype-rest (ref), wildtype-exercise, knockout-rest, knockout-exercise), and metabolite.
Fixed effects for Aim 2 were condition (wildtype-exercise (ref), knockout-exercise, knockout-C7, knockout-C8), and metabolite.
All 2-way interactions between fixed effects were included in the models.
Animal ID was the random effect.
A general correlation structure was assumed.
Estimates for the contrasts comparing each combination of condition and metabolite are presented.
P-values were adjusted to control the false discovery rate, the expected proportion of false discoveries amongst the rejected hypotheses.
The data was analyzed using R version 3.6.0 (2019-04-26) and the `nlme` package version 3.1.137.

Estimate model.
Specify the correlation structure using `cs`.
Use `corCompSymm`, *compound symmetry structure corresponding to a constant correlation*.

**References**

Benjamini, Y., and Hochberg, Y.
(1995).
Controlling the false discovery rate: a practical and powerful approach to multiple testing.
*Journal of the Royal Statistical Society Series B* 57, 289–300.


```r
citation()
```

```
## 
## To cite R in publications use:
## 
##   R Core Team (2019). R: A language and environment for
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
##     year = {2019},
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
## Pinheiro J, Bates D, DebRoy S, Sarkar D, R Core Team (2018).
## _nlme: Linear and Nonlinear Mixed Effects Models_. R package
## version 3.1-137, <URL: https://CRAN.R-project.org/package=nlme>.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {{nlme}: Linear and Nonlinear Mixed Effects Models},
##     author = {Jose Pinheiro and Douglas Bates and Saikat DebRoy and Deepayan Sarkar and {R Core Team}},
##     year = {2018},
##     note = {R package version 3.1-137},
##     url = {https://CRAN.R-project.org/package=nlme},
##   }
```

Set some constants


```r
fixed <- formula(z_value ~
                   genotype +
                   condition +
                   metabolite +
                   genotype * condition +
                   genotype * metabolite +
                   condition * metabolite +
                   genotype * condition * metabolite)
random <- formula(~ 1 | id)
ctrl <- lmeControl(opt = "optim",
                   maxIter = 500, msMaxIter = 500,
                   tolerance = 1e-6, niterEM = 25, msMaxEval = 200, msTol = 1e-7)
```


## Aim 1

* WT-rest (ref)
* WT-exercise
* KO-rest
* KO-exercise

### Aim 1-a: Amino Acids and TCA Cycle intermediates (Neg)


```r
t0 <- Sys.time()
M <- estimateModel(data = D1a, fixed, random)
Sys.time() - t0
```

```
## Time difference of 1.744354 secs
```

```r
M %>% plot()
```

![plot of chunk lmeDiagnosticAim1a](figures/lmeDiagnosticAim1a-1.png)

```r
M %>% ranef() %>% plot()
```

![plot of chunk lmeDiagnosticAim1a](figures/lmeDiagnosticAim1a-2.png)

```r
M %>% anova() %>% kable()
```



|                              | numDF| denDF|    F-value|   p-value|
|:-----------------------------|-----:|-----:|----------:|---------:|
|(Intercept)                   |     1|   850|  0.6111944| 0.4345562|
|genotype                      |     1|    38|  0.0258576| 0.8731004|
|condition                     |     1|    38|  0.0136006| 0.9077740|
|metabolite                    |    23|   850| 10.5205817| 0.0000000|
|genotype:condition            |     1|    38|  1.1042985| 0.2999619|
|genotype:metabolite           |    23|   850|  3.2640881| 0.0000004|
|condition:metabolite          |    23|   850|  7.1521169| 0.0000000|
|genotype:condition:metabolite |    23|   850|  1.4379894| 0.0838042|

```r
fixef <- M %>% tidy(effects = "fixed") %>% filter(grepl("condition", term)) %>% adjustPvalue()
fixef %>% filter(sig) %>% kable()
```



|term                                   | estimate| std.error| statistic|   p.value| p.adjustBH|sig  |
|:--------------------------------------|--------:|---------:|---------:|---------:|----------:|:----|
|conditionExercise:metaboliteTryptophan | 3.233007| 0.8846537|  3.654545| 0.0002734|  0.0131252|TRUE |

```r
fixef %>% write.csv(file = "data/processed/lmeFixedCoefAim1a.csv", row.names = FALSE)
M1a <- M
```

Calculate contrasts of *Exercise vs Rest* given metabolite and genotype.


```r
metabolites <- L1a[["data"]] %>% pull(metabolite) %>% levels()
Ftests <- runClusters(D1a, metabolites, fixed, random, "condition", "Exercise", ctrl)
```

```
## Loading required package: data.table
```

```
## 
## Attaching package: 'data.table'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     between, first, last
```

```r
Ftests %>% kable(digits = 5)
```



|contrast |metabolite       |genotype |     beta|      se| numDF| denDF|  F.value| p.value|  lowerCL|  upperCL| p.adjustBH|sig   |
|:--------|:----------------|:--------|--------:|-------:|-----:|-----:|--------:|-------:|--------:|--------:|----------:|:-----|
|Exercise |alaninesarcosine |WT       | -0.37426| 0.65644|     1|    38|  0.32507| 0.57194| -1.66086|  0.91233|    0.80434|FALSE |
|Exercise |arginine         |WT       | -1.20541| 0.65644|     1|    38|  3.37197| 0.07415| -2.49200|  0.08118|    0.23727|FALSE |
|Exercise |asparagine       |WT       | -0.34403| 0.65644|     1|    38|  0.27468| 0.60326| -1.63063|  0.94256|    0.80434|FALSE |
|Exercise |AsparticAcid     |WT       | -0.63171| 0.67086|     1|    38|  0.88670| 0.35232| -1.94656|  0.68315|    0.67646|FALSE |
|Exercise |CITRIC           |WT       |  1.08549| 0.67417|     1|    38|  2.59245| 0.11565| -0.23586|  2.40685|    0.32910|FALSE |
|Exercise |FUMARIC          |WT       |  0.03760| 0.67417|     1|    38|  0.00311| 0.95582| -1.28376|  1.35896|    0.97615|FALSE |
|Exercise |glutamic         |WT       | -0.00921| 0.65644|     1|    38|  0.00020| 0.98888| -1.29580|  1.27738|    0.98888|FALSE |
|Exercise |glutamine        |WT       | -1.90178| 0.67411|     1|    38|  7.95893| 0.00757| -3.22302| -0.58054|    0.04540|TRUE  |
|Exercise |glycine          |WT       | -0.63369| 0.65644|     1|    38|  0.93191| 0.34047| -1.92029|  0.65290|    0.67646|FALSE |
|Exercise |histidine        |WT       | -1.40926| 0.65644|     1|    38|  4.60890| 0.03825| -2.69585| -0.12267|    0.14124|FALSE |
|Exercise |isoleucine       |WT       | -0.07105| 0.65595|     1|    38|  0.01173| 0.91432| -1.35669|  1.21460|    0.97528|FALSE |
|Exercise |LACTIC           |WT       |  0.05173| 0.67417|     1|    38|  0.00589| 0.93924| -1.26963|  1.37309|    0.97615|FALSE |
|Exercise |leucine          |WT       | -0.73412| 0.65595|     1|    38|  1.25254| 0.27009| -2.01977|  0.55152|    0.56367|FALSE |
|Exercise |lysine           |WT       | -2.17233| 0.67411|     1|    38| 10.38443| 0.00261| -3.49357| -0.85109|    0.01789|TRUE  |
|Exercise |MALIC            |WT       |  1.74236| 0.67417|     1|    38|  6.67928| 0.01372|  0.42100|  3.06371|    0.07260|FALSE |
|Exercise |methionine       |WT       | -1.06207| 0.67408|     1|    38|  2.48247| 0.12341| -2.38325|  0.25910|    0.32910|FALSE |
|Exercise |METHYLSUCCINIC   |WT       | -0.58806| 0.70970|     1|    38|  0.68658| 0.41250| -1.97904|  0.80293|    0.74413|FALSE |
|Exercise |Phenylalanine    |WT       |  0.40919| 0.68760|     1|    38|  0.35414| 0.55531| -0.93848|  1.75685|    0.80434|FALSE |
|Exercise |serine           |WT       | -0.86079| 0.65644|     1|    38|  1.71951| 0.19762| -2.14738|  0.42581|    0.49926|FALSE |
|Exercise |SUCCINIC2        |WT       | -0.35796| 0.67417|     1|    38|  0.28192| 0.59853| -1.67932|  0.96340|    0.80434|FALSE |
|Exercise |threonine        |WT       | -1.52864| 0.67086|     1|    38|  5.19220| 0.02840| -2.84350| -0.21379|    0.12391|FALSE |
|Exercise |Tryptophan       |WT       |  2.85874| 0.65644|     1|    38| 18.96551| 0.00010|  1.57215|  4.14534|    0.00233|TRUE  |
|Exercise |Tyrosine         |WT       |  2.27737| 0.65595|     1|    38| 12.05363| 0.00130|  0.99172|  3.56301|    0.01044|TRUE  |
|Exercise |valine           |WT       |  0.51861| 0.65595|     1|    38|  0.62508| 0.43407| -0.76704|  1.80426|    0.74413|FALSE |
|Exercise |alaninesarcosine |KO       | -0.46080| 0.65395|     1|    38|  0.49653| 0.48533| -1.74251|  0.82091|    0.78137|FALSE |
|Exercise |arginine         |KO       | -0.80571| 0.65367|     1|    38|  1.51930| 0.22530| -2.08688|  0.47546|    0.51497|FALSE |
|Exercise |asparagine       |KO       |  0.12812| 0.65312|     1|    38|  0.03848| 0.84552| -1.15197|  1.40821|    0.92239|FALSE |
|Exercise |AsparticAcid     |KO       |  0.21643| 0.65362|     1|    38|  0.10964| 0.74237| -1.06465|  1.49751|    0.89543|FALSE |
|Exercise |CITRIC           |KO       | -0.41523| 0.66783|     1|    38|  0.38659| 0.53781| -1.72414|  0.89368|    0.80434|FALSE |
|Exercise |FUMARIC          |KO       |  0.18169| 0.67117|     1|    38|  0.07328| 0.78809| -1.13379|  1.49716|    0.90067|FALSE |
|Exercise |glutamic         |KO       | -0.45677| 0.65278|     1|    38|  0.48962| 0.48836| -1.73621|  0.82266|    0.78137|FALSE |
|Exercise |glutamine        |KO       | -1.65687| 0.65114|     1|    38|  6.47490| 0.01513| -2.93307| -0.38067|    0.07260|FALSE |
|Exercise |glycine          |KO       | -0.20159| 0.66912|     1|    38|  0.09076| 0.76485| -1.51304|  1.10986|    0.89543|FALSE |
|Exercise |histidine        |KO       | -1.02948| 0.64606|     1|    38|  2.53920| 0.11934| -2.29573|  0.23676|    0.32910|FALSE |
|Exercise |isoleucine       |KO       |  0.51308| 0.64786|     1|    38|  0.62720| 0.43330| -0.75670|  1.78285|    0.74413|FALSE |
|Exercise |LACTIC           |KO       | -0.13608| 0.65540|     1|    38|  0.04311| 0.83663| -1.42065|  1.14849|    0.92239|FALSE |
|Exercise |leucine          |KO       |  2.65159| 0.66366|     1|    38| 15.96305| 0.00029|  1.35083|  3.95235|    0.00321|TRUE  |
|Exercise |lysine           |KO       | -1.41055| 0.64700|     1|    38|  4.75299| 0.03551| -2.67864| -0.14245|    0.14124|FALSE |
|Exercise |MALIC            |KO       |  0.30216| 0.66593|     1|    38|  0.20588| 0.65259| -1.00303|  1.60736|    0.82433|FALSE |
|Exercise |methionine       |KO       | -0.82512| 0.64705|     1|    38|  1.62613| 0.20998| -2.09332|  0.44308|    0.50395|FALSE |
|Exercise |METHYLSUCCINIC   |KO       |  0.21199| 0.65473|     1|    38|  0.10483| 0.74788| -1.07125|  1.49522|    0.89543|FALSE |
|Exercise |Phenylalanine    |KO       |  0.75605| 0.66271|     1|    38|  1.30156| 0.26107| -0.54283|  2.05494|    0.56367|FALSE |
|Exercise |serine           |KO       | -0.31942| 0.66518|     1|    38|  0.23059| 0.63384| -1.62315|  0.98431|    0.82227|FALSE |
|Exercise |SUCCINIC2        |KO       | -1.25729| 0.65510|     1|    38|  3.68346| 0.06249| -2.54126|  0.02668|    0.21424|FALSE |
|Exercise |threonine        |KO       |  0.37703| 0.64778|     1|    38|  0.33876| 0.56398| -0.89260|  1.64666|    0.80434|FALSE |
|Exercise |Tryptophan       |KO       |  2.55778| 0.64869|     1|    38| 15.54712| 0.00033|  1.28637|  3.82919|    0.00321|TRUE  |
|Exercise |Tyrosine         |KO       |  3.80816| 0.68056|     1|    38| 31.31058| 0.00000|  2.47427|  5.14204|    0.00010|TRUE  |
|Exercise |valine           |KO       |  2.69917| 0.64774|     1|    38| 17.36429| 0.00017|  1.42962|  3.96872|    0.00274|TRUE  |

```r
Ftests %>% write.csv(file = "data/processed/contrastsAim1a.csv", row.names = FALSE)
```

### Aim 1-b: Acylcarnitines


```r
t0 <- Sys.time()
M <- estimateModel(data = D1b, fixed, random)
Sys.time() - t0
```

```
## Time difference of 0.5964379 secs
```

```r
M %>% plot()
```

![plot of chunk lmeDiagnosticAim1b](figures/lmeDiagnosticAim1b-1.png)

```r
M %>% ranef() %>% plot()
```

![plot of chunk lmeDiagnosticAim1b](figures/lmeDiagnosticAim1b-2.png)

```r
M %>% anova() %>% kable()
```



|                              | numDF| denDF|    F-value|   p-value|
|:-----------------------------|-----:|-----:|----------:|---------:|
|(Intercept)                   |     1|   773| 28.4004530| 0.0000001|
|genotype                      |     1|    38| 13.4275328| 0.0007537|
|condition                     |     1|    38|  0.0005077| 0.9821416|
|metabolite                    |    21|   773| 87.1190391| 0.0000000|
|genotype:condition            |     1|    38|  1.9545142| 0.1702088|
|genotype:metabolite           |    21|   773| 96.9402594| 0.0000000|
|condition:metabolite          |    21|   773|  2.2555214| 0.0010851|
|genotype:condition:metabolite |    21|   773|  1.3496214| 0.1352562|

```r
fixef <- M %>% tidy(effects = "fixed") %>% filter(grepl("condition", term)) %>% adjustPvalue()
fixef %>% filter(sig) %>% kable()
```



|term | estimate| std.error| statistic| p.value| p.adjustBH|sig |
|:----|--------:|---------:|---------:|-------:|----------:|:---|

```r
fixef %>% write.csv(file = "data/processed/lmeFixedCoefAim1b.csv", row.names = FALSE)
M1b <- M
```

Calculate contrasts of *Exercise vs Rest* given metabolite and genotype.


```r
metabolites <- L1b[["data"]] %>% pull(metabolite) %>% levels()
Ftests <- runClusters(D1b, metabolites, fixed, random, "condition", "Exercise", ctrl)
Ftests %>% kable(digits = 5)
```



|contrast |metabolite              |genotype |     beta|      se| numDF| denDF|  F.value| p.value|  lowerCL|  upperCL| p.adjustBH|sig   |
|:--------|:-----------------------|:--------|--------:|-------:|-----:|-----:|--------:|-------:|--------:|--------:|----------:|:-----|
|Exercise |2methylbutyrylcarnitine |WT       |  2.34938| 0.85655|     1|    38|  7.52319| 0.00924|  0.67057|  4.02819|    0.13552|FALSE |
|Exercise |3HMG                    |WT       | -0.38240| 0.80788|     1|    38|  0.22405| 0.63868| -1.96581|  1.20101|    0.95754|FALSE |
|Exercise |acetylcarnitine         |WT       | -0.22584| 0.80788|     1|    38|  0.07815| 0.78134| -1.80925|  1.35757|    0.98001|FALSE |
|Exercise |butyrylcarnitine        |WT       | -0.64051| 0.80788|     1|    38|  0.62858| 0.43280| -2.22392|  0.94290|    0.94452|FALSE |
|Exercise |C101total               |WT       |  0.36503| 0.84741|     1|    38|  0.18555| 0.66908| -1.29586|  2.02592|    0.95754|FALSE |
|Exercise |C12                     |WT       | -0.15053| 0.80788|     1|    38|  0.03472| 0.85318| -1.73394|  1.43288|    0.98001|FALSE |
|Exercise |C141total               |WT       | -0.27016| 0.80788|     1|    38|  0.11182| 0.73992| -1.85357|  1.31326|    0.95754|FALSE |
|Exercise |C142total               |WT       | -0.49882| 0.80788|     1|    38|  0.38124| 0.54062| -2.08223|  1.08459|    0.95150|FALSE |
|Exercise |C151total               |WT       | -0.38852| 0.80788|     1|    38|  0.23128| 0.63333| -1.97193|  1.19489|    0.95754|FALSE |
|Exercise |C161total               |WT       | -0.90280| 0.82991|     1|    38|  1.18337| 0.28352| -2.52939|  0.72380|    0.90928|FALSE |
|Exercise |C171total               |WT       | -0.56444| 0.80788|     1|    38|  0.48814| 0.48901| -2.14786|  1.01897|    0.94452|FALSE |
|Exercise |C181total               |WT       | -0.54518| 0.82991|     1|    38|  0.43154| 0.51519| -2.17178|  1.08141|    0.94452|FALSE |
|Exercise |C182total               |WT       | -1.30465| 0.82991|     1|    38|  2.47128| 0.12423| -2.93124|  0.32195|    0.60736|FALSE |
|Exercise |C191total               |WT       | -0.06808| 0.82991|     1|    38|  0.00673| 0.93505| -1.69468|  1.55851|    0.98001|FALSE |
|Exercise |carnitine               |WT       |  1.39532| 0.80788|     1|    38|  2.98300| 0.09226| -0.18810|  2.97873|    0.50746|FALSE |
|Exercise |ethylmalonylcarnitine   |WT       |  2.72676| 0.80788|     1|    38| 11.39206| 0.00171|  1.14335|  4.31018|    0.03764|TRUE  |
|Exercise |isobutyrylcarnitine     |WT       |  1.02763| 0.84741|     1|    38|  1.47056| 0.23274| -0.63327|  2.68852|    0.90274|FALSE |
|Exercise |methylsuccinylcarnitine |WT       |  3.62902| 0.80788|     1|    38| 20.17833| 0.00006|  2.04560|  5.21243|    0.00281|TRUE  |
|Exercise |ndecanoylcarnitine      |WT       | -0.76046| 0.80788|     1|    38|  0.88606| 0.35249| -2.34388|  0.82295|    0.94452|FALSE |
|Exercise |nhexanoylcarnitine      |WT       |  0.07801| 0.80788|     1|    38|  0.00932| 0.92358| -1.50540|  1.66142|    0.98001|FALSE |
|Exercise |noctanoylcarnitine      |WT       |  0.63598| 0.80788|     1|    38|  0.61972| 0.43603| -0.94743|  2.21939|    0.94452|FALSE |
|Exercise |propionylcarnitine      |WT       |  1.78886| 0.80788|     1|    38|  4.90297| 0.03288|  0.20545|  3.37227|    0.28938|FALSE |
|Exercise |2methylbutyrylcarnitine |KO       | -0.09179| 0.89069|     1|    38|  0.01062| 0.91846| -1.83751|  1.65393|    0.98001|FALSE |
|Exercise |3HMG                    |KO       | -0.18368| 0.83102|     1|    38|  0.04886| 0.82625| -1.81246|  1.44510|    0.98001|FALSE |
|Exercise |acetylcarnitine         |KO       | -0.83419| 0.80908|     1|    38|  1.06304| 0.30904| -2.41996|  0.75158|    0.90928|FALSE |
|Exercise |butyrylcarnitine        |KO       | -1.48886| 0.80905|     1|    38|  3.38660| 0.07355| -3.07456|  0.09684|    0.50746|FALSE |
|Exercise |C101total               |KO       |  0.30529| 0.80790|     1|    38|  0.14280| 0.70762| -1.27817|  1.88876|    0.95754|FALSE |
|Exercise |C12                     |KO       |  0.08062| 0.83091|     1|    38|  0.00941| 0.92321| -1.54793|  1.70917|    0.98001|FALSE |
|Exercise |C141total               |KO       |  0.04435| 0.83145|     1|    38|  0.00285| 0.95774| -1.58527|  1.67396|    0.98001|FALSE |
|Exercise |C142total               |KO       | -0.53978| 0.80912|     1|    38|  0.44504| 0.50873| -2.12562|  1.04607|    0.94452|FALSE |
|Exercise |C151total               |KO       | -0.55002| 0.80806|     1|    38|  0.46330| 0.50021| -2.13380|  1.03376|    0.94452|FALSE |
|Exercise |C161total               |KO       |  0.46848| 0.80808|     1|    38|  0.33611| 0.56550| -1.11532|  2.05228|    0.95701|FALSE |
|Exercise |C171total               |KO       | -0.61092| 0.80969|     1|    38|  0.56928| 0.45519| -2.19788|  0.97605|    0.94452|FALSE |
|Exercise |C181total               |KO       |  0.83388| 0.81037|     1|    38|  1.05886| 0.30998| -0.75442|  2.42217|    0.90928|FALSE |
|Exercise |C182total               |KO       | -1.89566| 0.80969|     1|    38|  5.48127| 0.02457| -3.48263| -0.30869|    0.27024|FALSE |
|Exercise |C191total               |KO       | -0.95253| 0.80874|     1|    38|  1.38720| 0.24620| -2.53762|  0.63257|    0.90274|FALSE |
|Exercise |carnitine               |KO       | -0.57328| 0.82654|     1|    38|  0.48108| 0.49215| -2.19326|  1.04670|    0.94452|FALSE |
|Exercise |ethylmalonylcarnitine   |KO       | -0.05218| 0.80915|     1|    38|  0.00416| 0.94892| -1.63809|  1.53374|    0.98001|FALSE |
|Exercise |isobutyrylcarnitine     |KO       |  1.56599| 0.89089|     1|    38|  3.08982| 0.08684| -0.18012|  3.31209|    0.50746|FALSE |
|Exercise |methylsuccinylcarnitine |KO       |  0.01727| 0.82712|     1|    38|  0.00044| 0.98345| -1.60386|  1.63839|    0.98345|FALSE |
|Exercise |ndecanoylcarnitine      |KO       | -0.32063| 0.82639|     1|    38|  0.15053| 0.70019| -1.94031|  1.29906|    0.95754|FALSE |
|Exercise |nhexanoylcarnitine      |KO       | -1.00752| 0.80899|     1|    38|  1.55103| 0.22061| -2.59312|  0.57808|    0.90274|FALSE |
|Exercise |noctanoylcarnitine      |KO       | -0.30141| 0.84883|     1|    38|  0.12609| 0.72448| -1.96509|  1.36226|    0.95754|FALSE |
|Exercise |propionylcarnitine      |KO       | -0.34830| 0.82620|     1|    38|  0.17772| 0.67571| -1.96762|  1.27102|    0.95754|FALSE |

```r
Ftests %>% write.csv(file = "data/processed/contrastsAim1b.csv", row.names = FALSE)
```


## Aim 2

* WT-exer (ref)
* KO-C7
* KO-C8
* KO-exer

### Aim 2-a: Amino Acids and TCA Cycle intermediates (Neg)


```r
t0 <- Sys.time()
M <- estimateModel(data = D2a, fixed, random)
```

```
## Error in MEEM(object, conLin, control$niterEM): Singularity in backsolve at level 0, block 1
```

```r
Sys.time() - t0
```

```
## Time difference of 0.117888 secs
```

```r
M %>% plot()
```

![plot of chunk lmeDiagnosticAim2a](figures/lmeDiagnosticAim2a-1.png)

```r
M %>% ranef() %>% plot()
```

![plot of chunk lmeDiagnosticAim2a](figures/lmeDiagnosticAim2a-2.png)

```r
M %>% anova() %>% kable()
```



|                              | numDF| denDF|    F-value|   p-value|
|:-----------------------------|-----:|-----:|----------:|---------:|
|(Intercept)                   |     1|   773| 28.4004530| 0.0000001|
|genotype                      |     1|    38| 13.4275328| 0.0007537|
|condition                     |     1|    38|  0.0005077| 0.9821416|
|metabolite                    |    21|   773| 87.1190391| 0.0000000|
|genotype:condition            |     1|    38|  1.9545142| 0.1702088|
|genotype:metabolite           |    21|   773| 96.9402594| 0.0000000|
|condition:metabolite          |    21|   773|  2.2555214| 0.0010851|
|genotype:condition:metabolite |    21|   773|  1.3496214| 0.1352562|

```r
fixef <- M %>% tidy(effects = "fixed") %>% filter(grepl("condition", term)) %>% adjustPvalue()
fixef %>% filter(sig) %>% kable()
```



|term | estimate| std.error| statistic| p.value| p.adjustBH|sig |
|:----|--------:|---------:|---------:|-------:|----------:|:---|

```r
fixef %>% write.csv(file = "data/processed/lmeFixedCoefAim2a.csv", row.names = FALSE)
M2a <- M
```

Calculate contrasts of *C7 vs Exercise* given metabolite and genotype.


```r
metabolites <- L2a[["data"]] %>% pull(metabolite) %>% levels()
Ftests <- runClusters(D2a, metabolites, fixed, random, "condition", "C7", ctrl)
```

```
## Error in {: task 1 failed - "Singularity in backsolve at level 0, block 1"
```

```r
Ftests %>% kable(digits = 5)
```



|contrast |metabolite              |genotype |     beta|      se| numDF| denDF|  F.value| p.value|  lowerCL|  upperCL| p.adjustBH|sig   |
|:--------|:-----------------------|:--------|--------:|-------:|-----:|-----:|--------:|-------:|--------:|--------:|----------:|:-----|
|Exercise |2methylbutyrylcarnitine |WT       |  2.34938| 0.85655|     1|    38|  7.52319| 0.00924|  0.67057|  4.02819|    0.13552|FALSE |
|Exercise |3HMG                    |WT       | -0.38240| 0.80788|     1|    38|  0.22405| 0.63868| -1.96581|  1.20101|    0.95754|FALSE |
|Exercise |acetylcarnitine         |WT       | -0.22584| 0.80788|     1|    38|  0.07815| 0.78134| -1.80925|  1.35757|    0.98001|FALSE |
|Exercise |butyrylcarnitine        |WT       | -0.64051| 0.80788|     1|    38|  0.62858| 0.43280| -2.22392|  0.94290|    0.94452|FALSE |
|Exercise |C101total               |WT       |  0.36503| 0.84741|     1|    38|  0.18555| 0.66908| -1.29586|  2.02592|    0.95754|FALSE |
|Exercise |C12                     |WT       | -0.15053| 0.80788|     1|    38|  0.03472| 0.85318| -1.73394|  1.43288|    0.98001|FALSE |
|Exercise |C141total               |WT       | -0.27016| 0.80788|     1|    38|  0.11182| 0.73992| -1.85357|  1.31326|    0.95754|FALSE |
|Exercise |C142total               |WT       | -0.49882| 0.80788|     1|    38|  0.38124| 0.54062| -2.08223|  1.08459|    0.95150|FALSE |
|Exercise |C151total               |WT       | -0.38852| 0.80788|     1|    38|  0.23128| 0.63333| -1.97193|  1.19489|    0.95754|FALSE |
|Exercise |C161total               |WT       | -0.90280| 0.82991|     1|    38|  1.18337| 0.28352| -2.52939|  0.72380|    0.90928|FALSE |
|Exercise |C171total               |WT       | -0.56444| 0.80788|     1|    38|  0.48814| 0.48901| -2.14786|  1.01897|    0.94452|FALSE |
|Exercise |C181total               |WT       | -0.54518| 0.82991|     1|    38|  0.43154| 0.51519| -2.17178|  1.08141|    0.94452|FALSE |
|Exercise |C182total               |WT       | -1.30465| 0.82991|     1|    38|  2.47128| 0.12423| -2.93124|  0.32195|    0.60736|FALSE |
|Exercise |C191total               |WT       | -0.06808| 0.82991|     1|    38|  0.00673| 0.93505| -1.69468|  1.55851|    0.98001|FALSE |
|Exercise |carnitine               |WT       |  1.39532| 0.80788|     1|    38|  2.98300| 0.09226| -0.18810|  2.97873|    0.50746|FALSE |
|Exercise |ethylmalonylcarnitine   |WT       |  2.72676| 0.80788|     1|    38| 11.39206| 0.00171|  1.14335|  4.31018|    0.03764|TRUE  |
|Exercise |isobutyrylcarnitine     |WT       |  1.02763| 0.84741|     1|    38|  1.47056| 0.23274| -0.63327|  2.68852|    0.90274|FALSE |
|Exercise |methylsuccinylcarnitine |WT       |  3.62902| 0.80788|     1|    38| 20.17833| 0.00006|  2.04560|  5.21243|    0.00281|TRUE  |
|Exercise |ndecanoylcarnitine      |WT       | -0.76046| 0.80788|     1|    38|  0.88606| 0.35249| -2.34388|  0.82295|    0.94452|FALSE |
|Exercise |nhexanoylcarnitine      |WT       |  0.07801| 0.80788|     1|    38|  0.00932| 0.92358| -1.50540|  1.66142|    0.98001|FALSE |
|Exercise |noctanoylcarnitine      |WT       |  0.63598| 0.80788|     1|    38|  0.61972| 0.43603| -0.94743|  2.21939|    0.94452|FALSE |
|Exercise |propionylcarnitine      |WT       |  1.78886| 0.80788|     1|    38|  4.90297| 0.03288|  0.20545|  3.37227|    0.28938|FALSE |
|Exercise |2methylbutyrylcarnitine |KO       | -0.09179| 0.89069|     1|    38|  0.01062| 0.91846| -1.83751|  1.65393|    0.98001|FALSE |
|Exercise |3HMG                    |KO       | -0.18368| 0.83102|     1|    38|  0.04886| 0.82625| -1.81246|  1.44510|    0.98001|FALSE |
|Exercise |acetylcarnitine         |KO       | -0.83419| 0.80908|     1|    38|  1.06304| 0.30904| -2.41996|  0.75158|    0.90928|FALSE |
|Exercise |butyrylcarnitine        |KO       | -1.48886| 0.80905|     1|    38|  3.38660| 0.07355| -3.07456|  0.09684|    0.50746|FALSE |
|Exercise |C101total               |KO       |  0.30529| 0.80790|     1|    38|  0.14280| 0.70762| -1.27817|  1.88876|    0.95754|FALSE |
|Exercise |C12                     |KO       |  0.08062| 0.83091|     1|    38|  0.00941| 0.92321| -1.54793|  1.70917|    0.98001|FALSE |
|Exercise |C141total               |KO       |  0.04435| 0.83145|     1|    38|  0.00285| 0.95774| -1.58527|  1.67396|    0.98001|FALSE |
|Exercise |C142total               |KO       | -0.53978| 0.80912|     1|    38|  0.44504| 0.50873| -2.12562|  1.04607|    0.94452|FALSE |
|Exercise |C151total               |KO       | -0.55002| 0.80806|     1|    38|  0.46330| 0.50021| -2.13380|  1.03376|    0.94452|FALSE |
|Exercise |C161total               |KO       |  0.46848| 0.80808|     1|    38|  0.33611| 0.56550| -1.11532|  2.05228|    0.95701|FALSE |
|Exercise |C171total               |KO       | -0.61092| 0.80969|     1|    38|  0.56928| 0.45519| -2.19788|  0.97605|    0.94452|FALSE |
|Exercise |C181total               |KO       |  0.83388| 0.81037|     1|    38|  1.05886| 0.30998| -0.75442|  2.42217|    0.90928|FALSE |
|Exercise |C182total               |KO       | -1.89566| 0.80969|     1|    38|  5.48127| 0.02457| -3.48263| -0.30869|    0.27024|FALSE |
|Exercise |C191total               |KO       | -0.95253| 0.80874|     1|    38|  1.38720| 0.24620| -2.53762|  0.63257|    0.90274|FALSE |
|Exercise |carnitine               |KO       | -0.57328| 0.82654|     1|    38|  0.48108| 0.49215| -2.19326|  1.04670|    0.94452|FALSE |
|Exercise |ethylmalonylcarnitine   |KO       | -0.05218| 0.80915|     1|    38|  0.00416| 0.94892| -1.63809|  1.53374|    0.98001|FALSE |
|Exercise |isobutyrylcarnitine     |KO       |  1.56599| 0.89089|     1|    38|  3.08982| 0.08684| -0.18012|  3.31209|    0.50746|FALSE |
|Exercise |methylsuccinylcarnitine |KO       |  0.01727| 0.82712|     1|    38|  0.00044| 0.98345| -1.60386|  1.63839|    0.98345|FALSE |
|Exercise |ndecanoylcarnitine      |KO       | -0.32063| 0.82639|     1|    38|  0.15053| 0.70019| -1.94031|  1.29906|    0.95754|FALSE |
|Exercise |nhexanoylcarnitine      |KO       | -1.00752| 0.80899|     1|    38|  1.55103| 0.22061| -2.59312|  0.57808|    0.90274|FALSE |
|Exercise |noctanoylcarnitine      |KO       | -0.30141| 0.84883|     1|    38|  0.12609| 0.72448| -1.96509|  1.36226|    0.95754|FALSE |
|Exercise |propionylcarnitine      |KO       | -0.34830| 0.82620|     1|    38|  0.17772| 0.67571| -1.96762|  1.27102|    0.95754|FALSE |

```r
Ftests %>% write.csv(file = "data/processed/contrastsAim2ai.csv", row.names = FALSE)
```

Calculate contrasts of *C8 vs Exercise* given metabolite and genotype.


```r
metabolites <- L2a[["data"]] %>% pull(metabolite) %>% levels()
Ftests <- runClusters(D2a, metabolites, fixed, random, "condition", "C8", ctrl)
```

```
## Error in {: task 1 failed - "Singularity in backsolve at level 0, block 1"
```

```r
Ftests %>% kable(digits = 5)
```



|contrast |metabolite              |genotype |     beta|      se| numDF| denDF|  F.value| p.value|  lowerCL|  upperCL| p.adjustBH|sig   |
|:--------|:-----------------------|:--------|--------:|-------:|-----:|-----:|--------:|-------:|--------:|--------:|----------:|:-----|
|Exercise |2methylbutyrylcarnitine |WT       |  2.34938| 0.85655|     1|    38|  7.52319| 0.00924|  0.67057|  4.02819|    0.13552|FALSE |
|Exercise |3HMG                    |WT       | -0.38240| 0.80788|     1|    38|  0.22405| 0.63868| -1.96581|  1.20101|    0.95754|FALSE |
|Exercise |acetylcarnitine         |WT       | -0.22584| 0.80788|     1|    38|  0.07815| 0.78134| -1.80925|  1.35757|    0.98001|FALSE |
|Exercise |butyrylcarnitine        |WT       | -0.64051| 0.80788|     1|    38|  0.62858| 0.43280| -2.22392|  0.94290|    0.94452|FALSE |
|Exercise |C101total               |WT       |  0.36503| 0.84741|     1|    38|  0.18555| 0.66908| -1.29586|  2.02592|    0.95754|FALSE |
|Exercise |C12                     |WT       | -0.15053| 0.80788|     1|    38|  0.03472| 0.85318| -1.73394|  1.43288|    0.98001|FALSE |
|Exercise |C141total               |WT       | -0.27016| 0.80788|     1|    38|  0.11182| 0.73992| -1.85357|  1.31326|    0.95754|FALSE |
|Exercise |C142total               |WT       | -0.49882| 0.80788|     1|    38|  0.38124| 0.54062| -2.08223|  1.08459|    0.95150|FALSE |
|Exercise |C151total               |WT       | -0.38852| 0.80788|     1|    38|  0.23128| 0.63333| -1.97193|  1.19489|    0.95754|FALSE |
|Exercise |C161total               |WT       | -0.90280| 0.82991|     1|    38|  1.18337| 0.28352| -2.52939|  0.72380|    0.90928|FALSE |
|Exercise |C171total               |WT       | -0.56444| 0.80788|     1|    38|  0.48814| 0.48901| -2.14786|  1.01897|    0.94452|FALSE |
|Exercise |C181total               |WT       | -0.54518| 0.82991|     1|    38|  0.43154| 0.51519| -2.17178|  1.08141|    0.94452|FALSE |
|Exercise |C182total               |WT       | -1.30465| 0.82991|     1|    38|  2.47128| 0.12423| -2.93124|  0.32195|    0.60736|FALSE |
|Exercise |C191total               |WT       | -0.06808| 0.82991|     1|    38|  0.00673| 0.93505| -1.69468|  1.55851|    0.98001|FALSE |
|Exercise |carnitine               |WT       |  1.39532| 0.80788|     1|    38|  2.98300| 0.09226| -0.18810|  2.97873|    0.50746|FALSE |
|Exercise |ethylmalonylcarnitine   |WT       |  2.72676| 0.80788|     1|    38| 11.39206| 0.00171|  1.14335|  4.31018|    0.03764|TRUE  |
|Exercise |isobutyrylcarnitine     |WT       |  1.02763| 0.84741|     1|    38|  1.47056| 0.23274| -0.63327|  2.68852|    0.90274|FALSE |
|Exercise |methylsuccinylcarnitine |WT       |  3.62902| 0.80788|     1|    38| 20.17833| 0.00006|  2.04560|  5.21243|    0.00281|TRUE  |
|Exercise |ndecanoylcarnitine      |WT       | -0.76046| 0.80788|     1|    38|  0.88606| 0.35249| -2.34388|  0.82295|    0.94452|FALSE |
|Exercise |nhexanoylcarnitine      |WT       |  0.07801| 0.80788|     1|    38|  0.00932| 0.92358| -1.50540|  1.66142|    0.98001|FALSE |
|Exercise |noctanoylcarnitine      |WT       |  0.63598| 0.80788|     1|    38|  0.61972| 0.43603| -0.94743|  2.21939|    0.94452|FALSE |
|Exercise |propionylcarnitine      |WT       |  1.78886| 0.80788|     1|    38|  4.90297| 0.03288|  0.20545|  3.37227|    0.28938|FALSE |
|Exercise |2methylbutyrylcarnitine |KO       | -0.09179| 0.89069|     1|    38|  0.01062| 0.91846| -1.83751|  1.65393|    0.98001|FALSE |
|Exercise |3HMG                    |KO       | -0.18368| 0.83102|     1|    38|  0.04886| 0.82625| -1.81246|  1.44510|    0.98001|FALSE |
|Exercise |acetylcarnitine         |KO       | -0.83419| 0.80908|     1|    38|  1.06304| 0.30904| -2.41996|  0.75158|    0.90928|FALSE |
|Exercise |butyrylcarnitine        |KO       | -1.48886| 0.80905|     1|    38|  3.38660| 0.07355| -3.07456|  0.09684|    0.50746|FALSE |
|Exercise |C101total               |KO       |  0.30529| 0.80790|     1|    38|  0.14280| 0.70762| -1.27817|  1.88876|    0.95754|FALSE |
|Exercise |C12                     |KO       |  0.08062| 0.83091|     1|    38|  0.00941| 0.92321| -1.54793|  1.70917|    0.98001|FALSE |
|Exercise |C141total               |KO       |  0.04435| 0.83145|     1|    38|  0.00285| 0.95774| -1.58527|  1.67396|    0.98001|FALSE |
|Exercise |C142total               |KO       | -0.53978| 0.80912|     1|    38|  0.44504| 0.50873| -2.12562|  1.04607|    0.94452|FALSE |
|Exercise |C151total               |KO       | -0.55002| 0.80806|     1|    38|  0.46330| 0.50021| -2.13380|  1.03376|    0.94452|FALSE |
|Exercise |C161total               |KO       |  0.46848| 0.80808|     1|    38|  0.33611| 0.56550| -1.11532|  2.05228|    0.95701|FALSE |
|Exercise |C171total               |KO       | -0.61092| 0.80969|     1|    38|  0.56928| 0.45519| -2.19788|  0.97605|    0.94452|FALSE |
|Exercise |C181total               |KO       |  0.83388| 0.81037|     1|    38|  1.05886| 0.30998| -0.75442|  2.42217|    0.90928|FALSE |
|Exercise |C182total               |KO       | -1.89566| 0.80969|     1|    38|  5.48127| 0.02457| -3.48263| -0.30869|    0.27024|FALSE |
|Exercise |C191total               |KO       | -0.95253| 0.80874|     1|    38|  1.38720| 0.24620| -2.53762|  0.63257|    0.90274|FALSE |
|Exercise |carnitine               |KO       | -0.57328| 0.82654|     1|    38|  0.48108| 0.49215| -2.19326|  1.04670|    0.94452|FALSE |
|Exercise |ethylmalonylcarnitine   |KO       | -0.05218| 0.80915|     1|    38|  0.00416| 0.94892| -1.63809|  1.53374|    0.98001|FALSE |
|Exercise |isobutyrylcarnitine     |KO       |  1.56599| 0.89089|     1|    38|  3.08982| 0.08684| -0.18012|  3.31209|    0.50746|FALSE |
|Exercise |methylsuccinylcarnitine |KO       |  0.01727| 0.82712|     1|    38|  0.00044| 0.98345| -1.60386|  1.63839|    0.98345|FALSE |
|Exercise |ndecanoylcarnitine      |KO       | -0.32063| 0.82639|     1|    38|  0.15053| 0.70019| -1.94031|  1.29906|    0.95754|FALSE |
|Exercise |nhexanoylcarnitine      |KO       | -1.00752| 0.80899|     1|    38|  1.55103| 0.22061| -2.59312|  0.57808|    0.90274|FALSE |
|Exercise |noctanoylcarnitine      |KO       | -0.30141| 0.84883|     1|    38|  0.12609| 0.72448| -1.96509|  1.36226|    0.95754|FALSE |
|Exercise |propionylcarnitine      |KO       | -0.34830| 0.82620|     1|    38|  0.17772| 0.67571| -1.96762|  1.27102|    0.95754|FALSE |

```r
Ftests %>% write.csv(file = "data/processed/contrastsAim2aii.csv", row.names = FALSE)
```

### Aim 2-b: Acylcarnitines


```r
t0 <- Sys.time()
M <- estimateModel(data = D2b, fixed, random)
```

```
## Error in MEEM(object, conLin, control$niterEM): Singularity in backsolve at level 0, block 1
```

```r
Sys.time() - t0
```

```
## Time difference of 0.1178899 secs
```

```r
M %>% plot()
```

![plot of chunk lmeDiagnosticAim2b](figures/lmeDiagnosticAim2b-1.png)

```r
M %>% ranef() %>% plot()
```

![plot of chunk lmeDiagnosticAim2b](figures/lmeDiagnosticAim2b-2.png)

```r
M %>% anova() %>% kable()
```



|                              | numDF| denDF|    F-value|   p-value|
|:-----------------------------|-----:|-----:|----------:|---------:|
|(Intercept)                   |     1|   773| 28.4004530| 0.0000001|
|genotype                      |     1|    38| 13.4275328| 0.0007537|
|condition                     |     1|    38|  0.0005077| 0.9821416|
|metabolite                    |    21|   773| 87.1190391| 0.0000000|
|genotype:condition            |     1|    38|  1.9545142| 0.1702088|
|genotype:metabolite           |    21|   773| 96.9402594| 0.0000000|
|condition:metabolite          |    21|   773|  2.2555214| 0.0010851|
|genotype:condition:metabolite |    21|   773|  1.3496214| 0.1352562|

```r
fixef <- M %>% tidy(effects = "fixed") %>% filter(grepl("condition", term)) %>% adjustPvalue()
fixef %>% filter(sig) %>% kable()
```



|term | estimate| std.error| statistic| p.value| p.adjustBH|sig |
|:----|--------:|---------:|---------:|-------:|----------:|:---|

```r
fixef %>% write.csv(file = "data/processed/lmeFixedCoefAim2b.csv", row.names = FALSE)
M2b <- M
```

Calculate contrasts of *C7 vs Exercise* given metabolite and genotype.


```r
metabolites <- L2b[["data"]] %>% pull(metabolite) %>% levels()
Ftests <- runClusters(D2b, metabolites, fixed, random, "condition", "C7", ctrl)
```

```
## Error in {: task 1 failed - "Singularity in backsolve at level 0, block 1"
```

```r
Ftests %>% kable(digits = 5)
```



|contrast |metabolite              |genotype |     beta|      se| numDF| denDF|  F.value| p.value|  lowerCL|  upperCL| p.adjustBH|sig   |
|:--------|:-----------------------|:--------|--------:|-------:|-----:|-----:|--------:|-------:|--------:|--------:|----------:|:-----|
|Exercise |2methylbutyrylcarnitine |WT       |  2.34938| 0.85655|     1|    38|  7.52319| 0.00924|  0.67057|  4.02819|    0.13552|FALSE |
|Exercise |3HMG                    |WT       | -0.38240| 0.80788|     1|    38|  0.22405| 0.63868| -1.96581|  1.20101|    0.95754|FALSE |
|Exercise |acetylcarnitine         |WT       | -0.22584| 0.80788|     1|    38|  0.07815| 0.78134| -1.80925|  1.35757|    0.98001|FALSE |
|Exercise |butyrylcarnitine        |WT       | -0.64051| 0.80788|     1|    38|  0.62858| 0.43280| -2.22392|  0.94290|    0.94452|FALSE |
|Exercise |C101total               |WT       |  0.36503| 0.84741|     1|    38|  0.18555| 0.66908| -1.29586|  2.02592|    0.95754|FALSE |
|Exercise |C12                     |WT       | -0.15053| 0.80788|     1|    38|  0.03472| 0.85318| -1.73394|  1.43288|    0.98001|FALSE |
|Exercise |C141total               |WT       | -0.27016| 0.80788|     1|    38|  0.11182| 0.73992| -1.85357|  1.31326|    0.95754|FALSE |
|Exercise |C142total               |WT       | -0.49882| 0.80788|     1|    38|  0.38124| 0.54062| -2.08223|  1.08459|    0.95150|FALSE |
|Exercise |C151total               |WT       | -0.38852| 0.80788|     1|    38|  0.23128| 0.63333| -1.97193|  1.19489|    0.95754|FALSE |
|Exercise |C161total               |WT       | -0.90280| 0.82991|     1|    38|  1.18337| 0.28352| -2.52939|  0.72380|    0.90928|FALSE |
|Exercise |C171total               |WT       | -0.56444| 0.80788|     1|    38|  0.48814| 0.48901| -2.14786|  1.01897|    0.94452|FALSE |
|Exercise |C181total               |WT       | -0.54518| 0.82991|     1|    38|  0.43154| 0.51519| -2.17178|  1.08141|    0.94452|FALSE |
|Exercise |C182total               |WT       | -1.30465| 0.82991|     1|    38|  2.47128| 0.12423| -2.93124|  0.32195|    0.60736|FALSE |
|Exercise |C191total               |WT       | -0.06808| 0.82991|     1|    38|  0.00673| 0.93505| -1.69468|  1.55851|    0.98001|FALSE |
|Exercise |carnitine               |WT       |  1.39532| 0.80788|     1|    38|  2.98300| 0.09226| -0.18810|  2.97873|    0.50746|FALSE |
|Exercise |ethylmalonylcarnitine   |WT       |  2.72676| 0.80788|     1|    38| 11.39206| 0.00171|  1.14335|  4.31018|    0.03764|TRUE  |
|Exercise |isobutyrylcarnitine     |WT       |  1.02763| 0.84741|     1|    38|  1.47056| 0.23274| -0.63327|  2.68852|    0.90274|FALSE |
|Exercise |methylsuccinylcarnitine |WT       |  3.62902| 0.80788|     1|    38| 20.17833| 0.00006|  2.04560|  5.21243|    0.00281|TRUE  |
|Exercise |ndecanoylcarnitine      |WT       | -0.76046| 0.80788|     1|    38|  0.88606| 0.35249| -2.34388|  0.82295|    0.94452|FALSE |
|Exercise |nhexanoylcarnitine      |WT       |  0.07801| 0.80788|     1|    38|  0.00932| 0.92358| -1.50540|  1.66142|    0.98001|FALSE |
|Exercise |noctanoylcarnitine      |WT       |  0.63598| 0.80788|     1|    38|  0.61972| 0.43603| -0.94743|  2.21939|    0.94452|FALSE |
|Exercise |propionylcarnitine      |WT       |  1.78886| 0.80788|     1|    38|  4.90297| 0.03288|  0.20545|  3.37227|    0.28938|FALSE |
|Exercise |2methylbutyrylcarnitine |KO       | -0.09179| 0.89069|     1|    38|  0.01062| 0.91846| -1.83751|  1.65393|    0.98001|FALSE |
|Exercise |3HMG                    |KO       | -0.18368| 0.83102|     1|    38|  0.04886| 0.82625| -1.81246|  1.44510|    0.98001|FALSE |
|Exercise |acetylcarnitine         |KO       | -0.83419| 0.80908|     1|    38|  1.06304| 0.30904| -2.41996|  0.75158|    0.90928|FALSE |
|Exercise |butyrylcarnitine        |KO       | -1.48886| 0.80905|     1|    38|  3.38660| 0.07355| -3.07456|  0.09684|    0.50746|FALSE |
|Exercise |C101total               |KO       |  0.30529| 0.80790|     1|    38|  0.14280| 0.70762| -1.27817|  1.88876|    0.95754|FALSE |
|Exercise |C12                     |KO       |  0.08062| 0.83091|     1|    38|  0.00941| 0.92321| -1.54793|  1.70917|    0.98001|FALSE |
|Exercise |C141total               |KO       |  0.04435| 0.83145|     1|    38|  0.00285| 0.95774| -1.58527|  1.67396|    0.98001|FALSE |
|Exercise |C142total               |KO       | -0.53978| 0.80912|     1|    38|  0.44504| 0.50873| -2.12562|  1.04607|    0.94452|FALSE |
|Exercise |C151total               |KO       | -0.55002| 0.80806|     1|    38|  0.46330| 0.50021| -2.13380|  1.03376|    0.94452|FALSE |
|Exercise |C161total               |KO       |  0.46848| 0.80808|     1|    38|  0.33611| 0.56550| -1.11532|  2.05228|    0.95701|FALSE |
|Exercise |C171total               |KO       | -0.61092| 0.80969|     1|    38|  0.56928| 0.45519| -2.19788|  0.97605|    0.94452|FALSE |
|Exercise |C181total               |KO       |  0.83388| 0.81037|     1|    38|  1.05886| 0.30998| -0.75442|  2.42217|    0.90928|FALSE |
|Exercise |C182total               |KO       | -1.89566| 0.80969|     1|    38|  5.48127| 0.02457| -3.48263| -0.30869|    0.27024|FALSE |
|Exercise |C191total               |KO       | -0.95253| 0.80874|     1|    38|  1.38720| 0.24620| -2.53762|  0.63257|    0.90274|FALSE |
|Exercise |carnitine               |KO       | -0.57328| 0.82654|     1|    38|  0.48108| 0.49215| -2.19326|  1.04670|    0.94452|FALSE |
|Exercise |ethylmalonylcarnitine   |KO       | -0.05218| 0.80915|     1|    38|  0.00416| 0.94892| -1.63809|  1.53374|    0.98001|FALSE |
|Exercise |isobutyrylcarnitine     |KO       |  1.56599| 0.89089|     1|    38|  3.08982| 0.08684| -0.18012|  3.31209|    0.50746|FALSE |
|Exercise |methylsuccinylcarnitine |KO       |  0.01727| 0.82712|     1|    38|  0.00044| 0.98345| -1.60386|  1.63839|    0.98345|FALSE |
|Exercise |ndecanoylcarnitine      |KO       | -0.32063| 0.82639|     1|    38|  0.15053| 0.70019| -1.94031|  1.29906|    0.95754|FALSE |
|Exercise |nhexanoylcarnitine      |KO       | -1.00752| 0.80899|     1|    38|  1.55103| 0.22061| -2.59312|  0.57808|    0.90274|FALSE |
|Exercise |noctanoylcarnitine      |KO       | -0.30141| 0.84883|     1|    38|  0.12609| 0.72448| -1.96509|  1.36226|    0.95754|FALSE |
|Exercise |propionylcarnitine      |KO       | -0.34830| 0.82620|     1|    38|  0.17772| 0.67571| -1.96762|  1.27102|    0.95754|FALSE |

```r
Ftests %>% write.csv(file = "data/processed/contrastsAim2bi.csv", row.names = FALSE)
```

Calculate contrasts of *C8 vs Exercise* given metabolite and genotype.


```r
metabolites <- L2b[["data"]] %>% pull(metabolite) %>% levels()
Ftests <- runClusters(D2b, metabolites, fixed, random, "condition", "C8", ctrl)
```

```
## Error in {: task 1 failed - "Singularity in backsolve at level 0, block 1"
```

```r
Ftests %>% kable(digits = 5)
```



|contrast |metabolite              |genotype |     beta|      se| numDF| denDF|  F.value| p.value|  lowerCL|  upperCL| p.adjustBH|sig   |
|:--------|:-----------------------|:--------|--------:|-------:|-----:|-----:|--------:|-------:|--------:|--------:|----------:|:-----|
|Exercise |2methylbutyrylcarnitine |WT       |  2.34938| 0.85655|     1|    38|  7.52319| 0.00924|  0.67057|  4.02819|    0.13552|FALSE |
|Exercise |3HMG                    |WT       | -0.38240| 0.80788|     1|    38|  0.22405| 0.63868| -1.96581|  1.20101|    0.95754|FALSE |
|Exercise |acetylcarnitine         |WT       | -0.22584| 0.80788|     1|    38|  0.07815| 0.78134| -1.80925|  1.35757|    0.98001|FALSE |
|Exercise |butyrylcarnitine        |WT       | -0.64051| 0.80788|     1|    38|  0.62858| 0.43280| -2.22392|  0.94290|    0.94452|FALSE |
|Exercise |C101total               |WT       |  0.36503| 0.84741|     1|    38|  0.18555| 0.66908| -1.29586|  2.02592|    0.95754|FALSE |
|Exercise |C12                     |WT       | -0.15053| 0.80788|     1|    38|  0.03472| 0.85318| -1.73394|  1.43288|    0.98001|FALSE |
|Exercise |C141total               |WT       | -0.27016| 0.80788|     1|    38|  0.11182| 0.73992| -1.85357|  1.31326|    0.95754|FALSE |
|Exercise |C142total               |WT       | -0.49882| 0.80788|     1|    38|  0.38124| 0.54062| -2.08223|  1.08459|    0.95150|FALSE |
|Exercise |C151total               |WT       | -0.38852| 0.80788|     1|    38|  0.23128| 0.63333| -1.97193|  1.19489|    0.95754|FALSE |
|Exercise |C161total               |WT       | -0.90280| 0.82991|     1|    38|  1.18337| 0.28352| -2.52939|  0.72380|    0.90928|FALSE |
|Exercise |C171total               |WT       | -0.56444| 0.80788|     1|    38|  0.48814| 0.48901| -2.14786|  1.01897|    0.94452|FALSE |
|Exercise |C181total               |WT       | -0.54518| 0.82991|     1|    38|  0.43154| 0.51519| -2.17178|  1.08141|    0.94452|FALSE |
|Exercise |C182total               |WT       | -1.30465| 0.82991|     1|    38|  2.47128| 0.12423| -2.93124|  0.32195|    0.60736|FALSE |
|Exercise |C191total               |WT       | -0.06808| 0.82991|     1|    38|  0.00673| 0.93505| -1.69468|  1.55851|    0.98001|FALSE |
|Exercise |carnitine               |WT       |  1.39532| 0.80788|     1|    38|  2.98300| 0.09226| -0.18810|  2.97873|    0.50746|FALSE |
|Exercise |ethylmalonylcarnitine   |WT       |  2.72676| 0.80788|     1|    38| 11.39206| 0.00171|  1.14335|  4.31018|    0.03764|TRUE  |
|Exercise |isobutyrylcarnitine     |WT       |  1.02763| 0.84741|     1|    38|  1.47056| 0.23274| -0.63327|  2.68852|    0.90274|FALSE |
|Exercise |methylsuccinylcarnitine |WT       |  3.62902| 0.80788|     1|    38| 20.17833| 0.00006|  2.04560|  5.21243|    0.00281|TRUE  |
|Exercise |ndecanoylcarnitine      |WT       | -0.76046| 0.80788|     1|    38|  0.88606| 0.35249| -2.34388|  0.82295|    0.94452|FALSE |
|Exercise |nhexanoylcarnitine      |WT       |  0.07801| 0.80788|     1|    38|  0.00932| 0.92358| -1.50540|  1.66142|    0.98001|FALSE |
|Exercise |noctanoylcarnitine      |WT       |  0.63598| 0.80788|     1|    38|  0.61972| 0.43603| -0.94743|  2.21939|    0.94452|FALSE |
|Exercise |propionylcarnitine      |WT       |  1.78886| 0.80788|     1|    38|  4.90297| 0.03288|  0.20545|  3.37227|    0.28938|FALSE |
|Exercise |2methylbutyrylcarnitine |KO       | -0.09179| 0.89069|     1|    38|  0.01062| 0.91846| -1.83751|  1.65393|    0.98001|FALSE |
|Exercise |3HMG                    |KO       | -0.18368| 0.83102|     1|    38|  0.04886| 0.82625| -1.81246|  1.44510|    0.98001|FALSE |
|Exercise |acetylcarnitine         |KO       | -0.83419| 0.80908|     1|    38|  1.06304| 0.30904| -2.41996|  0.75158|    0.90928|FALSE |
|Exercise |butyrylcarnitine        |KO       | -1.48886| 0.80905|     1|    38|  3.38660| 0.07355| -3.07456|  0.09684|    0.50746|FALSE |
|Exercise |C101total               |KO       |  0.30529| 0.80790|     1|    38|  0.14280| 0.70762| -1.27817|  1.88876|    0.95754|FALSE |
|Exercise |C12                     |KO       |  0.08062| 0.83091|     1|    38|  0.00941| 0.92321| -1.54793|  1.70917|    0.98001|FALSE |
|Exercise |C141total               |KO       |  0.04435| 0.83145|     1|    38|  0.00285| 0.95774| -1.58527|  1.67396|    0.98001|FALSE |
|Exercise |C142total               |KO       | -0.53978| 0.80912|     1|    38|  0.44504| 0.50873| -2.12562|  1.04607|    0.94452|FALSE |
|Exercise |C151total               |KO       | -0.55002| 0.80806|     1|    38|  0.46330| 0.50021| -2.13380|  1.03376|    0.94452|FALSE |
|Exercise |C161total               |KO       |  0.46848| 0.80808|     1|    38|  0.33611| 0.56550| -1.11532|  2.05228|    0.95701|FALSE |
|Exercise |C171total               |KO       | -0.61092| 0.80969|     1|    38|  0.56928| 0.45519| -2.19788|  0.97605|    0.94452|FALSE |
|Exercise |C181total               |KO       |  0.83388| 0.81037|     1|    38|  1.05886| 0.30998| -0.75442|  2.42217|    0.90928|FALSE |
|Exercise |C182total               |KO       | -1.89566| 0.80969|     1|    38|  5.48127| 0.02457| -3.48263| -0.30869|    0.27024|FALSE |
|Exercise |C191total               |KO       | -0.95253| 0.80874|     1|    38|  1.38720| 0.24620| -2.53762|  0.63257|    0.90274|FALSE |
|Exercise |carnitine               |KO       | -0.57328| 0.82654|     1|    38|  0.48108| 0.49215| -2.19326|  1.04670|    0.94452|FALSE |
|Exercise |ethylmalonylcarnitine   |KO       | -0.05218| 0.80915|     1|    38|  0.00416| 0.94892| -1.63809|  1.53374|    0.98001|FALSE |
|Exercise |isobutyrylcarnitine     |KO       |  1.56599| 0.89089|     1|    38|  3.08982| 0.08684| -0.18012|  3.31209|    0.50746|FALSE |
|Exercise |methylsuccinylcarnitine |KO       |  0.01727| 0.82712|     1|    38|  0.00044| 0.98345| -1.60386|  1.63839|    0.98345|FALSE |
|Exercise |ndecanoylcarnitine      |KO       | -0.32063| 0.82639|     1|    38|  0.15053| 0.70019| -1.94031|  1.29906|    0.95754|FALSE |
|Exercise |nhexanoylcarnitine      |KO       | -1.00752| 0.80899|     1|    38|  1.55103| 0.22061| -2.59312|  0.57808|    0.90274|FALSE |
|Exercise |noctanoylcarnitine      |KO       | -0.30141| 0.84883|     1|    38|  0.12609| 0.72448| -1.96509|  1.36226|    0.95754|FALSE |
|Exercise |propionylcarnitine      |KO       | -0.34830| 0.82620|     1|    38|  0.17772| 0.67571| -1.96762|  1.27102|    0.95754|FALSE |

```r
Ftests %>% write.csv(file = "data/processed/contrastsAim2bii.csv", row.names = FALSE)
```


## Save

Save `lme` objects for interactive work.


```r
save(M1a, M1b, M2a, M2b, file = "data/processed/lmeObjects.RData")
```
