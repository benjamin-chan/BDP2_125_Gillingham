---
title: "Metabolomics of very long-chain aclCoA dehydrogenase knockout mice"
date: "2019-05-16 20:46:19"
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
```

```
## Warning: package 'readxl' was built under R version 3.5.3
```

```r
library(magrittr)
```

```
## Warning: package 'magrittr' was built under R version 3.5.3
```

```r
library(dplyr)
```

```
## Warning: package 'dplyr' was built under R version 3.5.3
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
## Warning: package 'nlme' was built under R version 3.5.3
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
```

```
## Warning: package 'broom' was built under R version 3.5.3
```

```r
library(ggplot2)
```

```
## Warning: package 'ggplot2' was built under R version 3.5.3
```

```r
library(svglite)
```

```
## Warning: package 'svglite' was built under R version 3.5.3
```

```r
library(knitr)
library(doParallel)
```

```
## Warning: package 'doParallel' was built under R version 3.5.3
```

```
## Loading required package: foreach
```

```
## Warning: package 'foreach' was built under R version 3.5.3
```

```
## Loading required package: iterators
```

```
## Warning: package 'iterators' was built under R version 3.5.3
```

```
## Loading required package: parallel
```

Reproducibility steps.


```r
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
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
##  [1] Rcpp_1.0.1       cellranger_1.1.0 pillar_1.3.1     compiler_3.5.1  
##  [5] plyr_1.8.4       tools_3.5.1      digest_0.6.18    evaluate_0.13   
##  [9] tibble_2.1.1     gtable_0.3.0     lattice_0.20-35  pkgconfig_2.0.2 
## [13] rlang_0.3.3      xfun_0.5         withr_2.1.2      stringr_1.4.0   
## [17] gdtools_0.1.7    generics_0.0.2   grid_3.5.1       tidyselect_0.2.5
## [21] glue_1.3.1       R6_2.4.0         purrr_0.3.2      tidyr_0.8.3     
## [25] codetools_0.2-15 backports_1.1.3  scales_1.0.0     htmltools_0.3.6 
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

Import Aim 1 Amino Acids and TCA Cycle intermediates (Neg).


```r
f <- "data/raw/Z scores Aim 1 Neg-aa.xlsx"
levels <- c("wtrest", "wtex", "korest", "koex")
L1a <- importDataToList(f, levels)
L1a[["file"]]
```

```
## [1] "data/raw/Z scores Aim 1 Neg-aa.xlsx"
```

```r
L1a[["dim"]]
```

```
## [1] 984   4
```

```r
L1a[["names"]]
```

```
## [1] "condition"  "id"         "metabolite" "z_value"
```

```r
L1a[["head"]]
```

```
## # A tibble: 6 x 4
##   condition    id metabolite z_value
##   <fct>     <dbl> <fct>        <dbl>
## 1 KOEX       1010 FUMARIC     -0.934
## 2 KOEX       1012 FUMARIC     -1.21 
## 3 KOEX       1014 FUMARIC      0.196
## 4 KOEX       1017 FUMARIC      0.455
## 5 KOEX       1018 FUMARIC     -0.580
## 6 KOEX       1019 FUMARIC     -0.593
```

```r
D1a <- L1a[["data"]]
```

Import Aim 1 Acylcarnitines.


```r
f <- "data/raw/Z scores Aim 1 AC.xlsx"
L1b <- importDataToList(f, levels)
L1b[["file"]]
```

```
## [1] "data/raw/Z scores Aim 1 AC.xlsx"
```

```r
L1b[["dim"]]
```

```
## [1] 899   4
```

```r
L1b[["names"]]
```

```
## [1] "condition"  "id"         "metabolite" "z_value"
```

```r
L1b[["head"]]
```

```
## # A tibble: 6 x 4
##   condition    id metabolite      z_value
##   <fct>     <dbl> <fct>             <dbl>
## 1 KOEX       1010 acetylcarnitine   -1.84
## 2 KOEX       1012 acetylcarnitine   -2.15
## 3 KOEX       1014 acetylcarnitine   -2.04
## 4 KOEX       1017 acetylcarnitine   -2.43
## 5 KOEX       1018 acetylcarnitine   -2.49
## 6 KOEX       1019 acetylcarnitine   -1.83
```

```r
D1b <- L1b[["data"]]
```

Import Aim 2 Amino Acids and TCA Cycle intermediates (Neg).


```r
f <- "data/raw/Z scores Aim 2 Neg-aa.xlsx"
levels <- c("wtex", "koex", "koc7", "koc8")
L2a <- importDataToList(f, levels)
L2a[["file"]]
```

```
## [1] "data/raw/Z scores Aim 2 Neg-aa.xlsx"
```

```r
L2a[["dim"]]
```

```
## [1] 1035    4
```

```r
L2a[["names"]]
```

```
## [1] "condition"  "id"         "metabolite" "z_value"
```

```r
L2a[["head"]]
```

```
## # A tibble: 6 x 4
##   condition    id metabolite z_value
##   <fct>     <dbl> <fct>        <dbl>
## 1 KOC7       1120 FUMARIC      -3.23
## 2 KOC7       1126 FUMARIC      -3.22
## 3 KOC7       1127 FUMARIC      -3.49
## 4 KOC7       1128 FUMARIC      -2.31
## 5 KOC7       1142 FUMARIC      -1.73
## 6 KOC7       1144 FUMARIC      -1.71
```

```r
D2a <- L2a[["data"]]
```

Import Aim 2 Acylcarnitines.


```r
f <- "data/raw/Z scores Aim 2 AC.xlsx"
L2b <- importDataToList(f, levels)
L2b[["file"]]
```

```
## [1] "data/raw/Z scores Aim 2 AC.xlsx"
```

```r
L2b[["dim"]]
```

```
## [1] 948   4
```

```r
L2b[["names"]]
```

```
## [1] "condition"  "id"         "metabolite" "z_value"
```

```r
L2b[["head"]]
```

```
## # A tibble: 6 x 4
##   condition    id metabolite      z_value
##   <fct>     <dbl> <fct>             <dbl>
## 1 KOC7       1120 acetylcarnitine   -2.57
## 2 KOC7       1126 acetylcarnitine   -2.43
## 3 KOC7       1127 acetylcarnitine   -2.58
## 4 KOC7       1128 acetylcarnitine   -2.11
## 5 KOC7       1142 acetylcarnitine   -2.67
## 6 KOC7       1143 acetylcarnitine   -2.55
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
The data was analyzed using R version 3.5.1 (2018-07-02) and the `nlme` package version 3.1.137.

Estimate model.
Specify the correlation structure using `cs`.
Use `corCompSymm`, *compound symmetry structure corresponding to a constant correlation*.


```r
estimateModel <- function (data) {
  require(magrittr)
  require(nlme)
  fixed <- formula(z_value ~
                     condition +
                     metabolite +
                     condition * metabolite)
  random <- formula(~ 1 | id)
  ctrl <- lmeControl(opt = "optim",
                     maxIter = 500, msMaxIter = 500,
                     tolerance = 1e-6, niterEM = 25, msMaxEval = 200, msTol = 1e-7)
  cs <-
    corCompSymm(form = random, fixed = FALSE) %>%
    Initialize(data = data)
  Dim(cs)
  lme(fixed,
      data = data,
      random = random,
      correlation = cs,
      control = ctrl)
}
```

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
##   R Core Team (2018). R: A language and environment for
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
##     year = {2018},
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


## Aim 1

* WT-rest (ref)
* WT-exercise
* KO-rest
* KO-exercise

### Amino Acids and TCA Cycle intermediates (Neg)


```r
t0 <- Sys.time()
M <- estimateModel(data = D1a)
Sys.time() - t0
```

```
## Time difference of 9.595944 secs
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



|                     | numDF| denDF|    F-value|   p-value|
|:--------------------|-----:|-----:|----------:|---------:|
|(Intercept)          |     1|   850|  0.6111944| 0.4345562|
|condition            |     3|    38|  0.3369410| 0.7986948|
|metabolite           |    23|   850| 10.5263614| 0.0000000|
|condition:metabolite |    69|   850|  3.9513981| 0.0000000|

```r
M %>% tidy(effects = "fixed") %>% kable()
```



|term                                     |   estimate| std.error|  statistic|   p.value|
|:----------------------------------------|----------:|---------:|----------:|---------:|
|(Intercept)                              |  0.0000000| 0.4757741|  0.0000000| 1.0000000|
|conditionWTEX                            | -0.3742641| 0.6564367| -0.5701451| 0.5719352|
|conditionKOREST                          |  0.4617075| 0.6728462|  0.6862007| 0.4967524|
|conditionKOEX                            |  0.0009073| 0.6539455|  0.0013875| 0.9989002|
|metabolitearginine                       |  0.0000000| 0.6402653|  0.0000000| 1.0000000|
|metaboliteasparagine                     |  0.0000000| 0.6402653|  0.0000000| 1.0000000|
|metaboliteAspartic Acid                  |  0.0000000| 0.6402653|  0.0000000| 1.0000000|
|metaboliteCITRIC                         |  0.0134814| 0.6584386|  0.0204748| 0.9836694|
|metaboliteFUMARIC                        |  0.0134814| 0.6584386|  0.0204748| 0.9836694|
|metaboliteglutamic                       |  0.0000000| 0.6402653|  0.0000000| 1.0000000|
|metaboliteglutamine                      | -0.0395393| 0.6583777| -0.0600557| 0.9521254|
|metaboliteglycine                        |  0.0000000| 0.6402653|  0.0000000| 1.0000000|
|metabolitehistidine                      |  0.0000000| 0.6402653|  0.0000000| 1.0000000|
|metaboliteisoleucine                     |  0.0000000| 0.6402653|  0.0000000| 1.0000000|
|metaboliteLACTIC                         |  0.0134814| 0.6584386|  0.0204748| 0.9836694|
|metaboliteleucine                        |  0.0000000| 0.6402653|  0.0000000| 1.0000000|
|metabolitelysine                         | -0.0395393| 0.6583777| -0.0600557| 0.9521254|
|metaboliteMALIC                          |  0.0134814| 0.6584386|  0.0204748| 0.9836694|
|metabolitemethionine                     | -0.0103602| 0.6583435| -0.0157368| 0.9874481|
|metaboliteMETHYLSUCCINIC                 | -0.0293151| 0.6803888| -0.0430858| 0.9656432|
|metabolitePhenylalanine                  |  0.0562397| 0.6583435|  0.0854260| 0.9319428|
|metaboliteserine                         |  0.0000000| 0.6402653|  0.0000000| 1.0000000|
|metaboliteSUCCINIC-2                     |  0.0134814| 0.6584386|  0.0204748| 0.9836694|
|metabolitethreonine                      |  0.0000000| 0.6402653|  0.0000000| 1.0000000|
|metaboliteTryptophan                     |  0.0000000| 0.6402653|  0.0000000| 1.0000000|
|metaboliteTyrosine                       |  0.0000000| 0.6402653|  0.0000000| 1.0000000|
|metabolitevaline                         |  0.0000000| 0.6402653|  0.0000000| 1.0000000|
|conditionWTEX:metabolitearginine         | -0.8311461| 0.8846537| -0.9395159| 0.3477329|
|conditionKOREST:metabolitearginine       | -0.3461279| 0.9054719| -0.3822625| 0.7023622|
|conditionKOEX:metabolitearginine         | -0.6910406| 0.8774366| -0.7875675| 0.4311693|
|conditionWTEX:metaboliteasparagine       |  0.0302293| 0.8846537|  0.0341707| 0.9727490|
|conditionKOREST:metaboliteasparagine     |  0.0009625| 0.9054719|  0.0010630| 0.9991521|
|conditionKOEX:metaboliteasparagine       |  0.5898843| 0.8779711|  0.6718721| 0.5018477|
|conditionWTEX:metaboliteAspartic Acid    | -0.2574452| 0.8954426| -0.2875061| 0.7737949|
|conditionKOREST:metaboliteAspartic Acid  | -0.3657526| 0.9054719| -0.4039358| 0.6863615|
|conditionKOEX:metaboliteAspartic Acid    |  0.3114769| 0.8749585|  0.3559904| 0.7219361|
|conditionWTEX:metaboliteCITRIC           |  1.4597571| 0.8978941|  1.6257565| 0.1043721|
|conditionKOREST:metaboliteCITRIC         | -1.2138172| 0.9184122| -1.3216475| 0.1866412|
|conditionKOEX:metaboliteCITRIC           | -1.1682481| 0.9017407| -1.2955477| 0.1954831|
|conditionWTEX:metaboliteFUMARIC          |  0.4118623| 0.8978941|  0.4586981| 0.6465682|
|conditionKOREST:metaboliteFUMARIC        | -1.2787945| 0.9311286| -1.3733812| 0.1699960|
|conditionKOEX:metaboliteFUMARIC          | -0.6363067| 0.8952350| -0.7107706| 0.4774214|
|conditionWTEX:metaboliteglutamic         |  0.3650516| 0.8846537|  0.4126491| 0.6799678|
|conditionKOREST:metaboliteglutamic       | -0.9216272| 0.9054719| -1.0178418| 0.3090427|
|conditionKOEX:metaboliteglutamic         | -0.9176003| 0.8772186| -1.0460338| 0.2958428|
|conditionWTEX:metaboliteglutamine        | -1.5275195| 0.8978494| -1.7013092| 0.0892507|
|conditionKOREST:metaboliteglutamine      | -0.4834446| 0.9183686| -0.5264167| 0.5987360|
|conditionKOEX:metaboliteglutamine        | -1.6795128| 0.8900555| -1.8869754| 0.0595042|
|conditionWTEX:metaboliteglycine          | -0.2594300| 0.8846537| -0.2932560| 0.7693981|
|conditionKOREST:metaboliteglycine        | -0.5030744| 0.9183675| -0.5477920| 0.5839785|
|conditionKOEX:metaboliteglycine          | -0.2438609| 0.8782050| -0.2776811| 0.7813247|
|conditionWTEX:metabolitehistidine        | -1.0349972| 0.8846537| -1.1699462| 0.2423506|
|conditionKOREST:metabolitehistidine      | -0.7441185| 0.9054719| -0.8218019| 0.4114201|
|conditionKOEX:metabolitehistidine        | -1.3127998| 0.8760541| -1.4985374| 0.1343649|
|conditionWTEX:metaboliteisoleucine       |  0.3032184| 0.8845989|  0.3427750| 0.7318525|
|conditionKOREST:metaboliteisoleucine     | -1.0609136| 0.9054719| -1.1716692| 0.2416582|
|conditionKOEX:metaboliteisoleucine       | -0.0870379| 0.8779974| -0.0991323| 0.9210566|
|conditionWTEX:metaboliteLACTIC           |  0.4259937| 0.8978941|  0.4744365| 0.6353105|
|conditionKOREST:metaboliteLACTIC         | -1.2519772| 0.9184122| -1.3631974| 0.1731813|
|conditionKOEX:metaboliteLACTIC           | -0.9272591| 0.8911827| -1.0404815| 0.2984121|
|conditionWTEX:metaboliteleucine          | -0.3598590| 0.8845989| -0.4068047| 0.6842538|
|conditionKOREST:metaboliteleucine        | -1.3278343| 0.9054719| -1.4664555| 0.1428941|
|conditionKOEX:metaboliteleucine          |  1.7845573| 0.8886678|  2.0081266| 0.0449458|
|conditionWTEX:metabolitelysine           | -1.7980628| 0.8978494| -2.0026329| 0.0455339|
|conditionKOREST:metabolitelysine         | -0.8627883| 0.9183686| -0.9394794| 0.3477516|
|conditionKOEX:metabolitelysine           | -1.8125358| 0.8898720| -2.0368501| 0.0419745|
|conditionWTEX:metaboliteMALIC            |  2.1166206| 0.8978941|  2.3573166| 0.0186334|
|conditionKOREST:metaboliteMALIC          | -1.1635017| 0.9184122| -1.2668622| 0.2055516|
|conditionKOEX:metaboliteMALIC            | -0.4005399| 0.9004875| -0.4448033| 0.6565751|
|conditionWTEX:metabolitemethionine       | -0.6878087| 0.8978244| -0.7660837| 0.4438392|
|conditionKOREST:metabolitemethionine     |  0.7677198| 0.9183441|  0.8359827| 0.4033996|
|conditionKOEX:metabolitemethionine       |  0.4034013| 0.8900149|  0.4532523| 0.6504827|
|conditionWTEX:metaboliteMETHYLSUCCINIC   | -0.2137941| 0.9245168| -0.2312495| 0.8171766|
|conditionKOREST:metaboliteMETHYLSUCCINIC | -1.4269004| 0.9342743| -1.5272820| 0.1270629|
|conditionKOEX:metaboliteMETHYLSUCCINIC   | -0.7541148| 0.9067013| -0.8317124| 0.4058049|
|conditionWTEX:metabolitePhenylalanine    |  0.7834502| 0.9083767|  0.8624728| 0.3886706|
|conditionKOREST:metabolitePhenylalanine  |  0.0257076| 0.9183441|  0.0279934| 0.9776740|
|conditionKOEX:metabolitePhenylalanine    |  1.2425627| 0.9022858|  1.3771276| 0.1688353|
|conditionWTEX:metaboliteserine           | -0.4865217| 0.8846537| -0.5499572| 0.5824932|
|conditionKOREST:metaboliteserine         | -0.1523948| 0.9183675| -0.1659410| 0.8682428|
|conditionKOEX:metaboliteserine           | -0.0110144| 0.8769417| -0.0125600| 0.9899818|
|conditionWTEX:metaboliteSUCCINIC-2       |  0.0163023| 0.8978941|  0.0181562| 0.9855185|
|conditionKOREST:metaboliteSUCCINIC-2     | -1.1103437| 0.9184122| -1.2089818| 0.2270060|
|conditionKOEX:metaboliteSUCCINIC-2       | -1.9068324| 0.8956432| -2.1290090| 0.0335404|
|conditionWTEX:metabolitethreonine        | -1.1543771| 0.8954426| -1.2891693| 0.1976899|
|conditionKOREST:metabolitethreonine      | -0.6311500| 0.9054719| -0.6970398| 0.4859684|
|conditionKOEX:metabolitethreonine        |  0.2066789| 0.8773770|  0.2355645| 0.8138274|
|conditionWTEX:metaboliteTryptophan       |  3.2330070| 0.8846537|  3.6545455| 0.0002734|
|conditionKOREST:metaboliteTryptophan     |  0.3144952| 0.9054719|  0.3473274| 0.7284313|
|conditionKOEX:metaboliteTryptophan       |  3.3330768| 0.8781535|  3.7955513| 0.0001578|
|conditionWTEX:metaboliteTyrosine         |  2.6516304| 0.8845989|  2.9975513| 0.0028008|
|conditionKOREST:metaboliteTyrosine       | -0.1932230| 0.9183429| -0.2104040| 0.8334028|
|conditionKOEX:metaboliteTyrosine         |  4.0757324| 0.8893361|  4.5828931| 0.0000053|
|conditionWTEX:metabolitevaline           |  0.8928731| 0.8845989|  1.0093537| 0.3130923|
|conditionKOREST:metabolitevaline         | -0.8085738| 0.9054719| -0.8929861| 0.3721174|
|conditionKOEX:metabolitevaline           |  2.3513989| 0.8776519|  2.6791930| 0.0075226|

```r
M %>% tidy(effects = "fixed") %>% write.csv(file = "data/processed/lmeFixedCoefAim1a.csv", row.names = FALSE)
M1a <- M
```

### Acylcarnitines


```r
t0 <- Sys.time()
M <- estimateModel(data = D1b)
Sys.time() - t0
```

```
## Time difference of 0.6244099 secs
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



|                     | numDF| denDF|   F-value|   p-value|
|:--------------------|-----:|-----:|---------:|---------:|
|(Intercept)          |     1|   773| 28.400453| 0.0000001|
|condition            |     3|    38|  5.345096| 0.0035915|
|metabolite           |    21|   773| 87.087957| 0.0000000|
|condition:metabolite |    63|   773| 33.515134| 0.0000000|

```r
M %>% tidy(effects = "fixed") %>% kable()
```



|term                                              |   estimate| std.error|  statistic|   p.value|
|:-------------------------------------------------|----------:|---------:|----------:|---------:|
|(Intercept)                                       |  0.0535035| 0.6502946|  0.0822759| 0.9344486|
|conditionWTEX                                     |  2.3493818| 0.8565498|  2.7428431| 0.0092402|
|conditionKOREST                                   |  0.1943739| 0.9506563|  0.2044629| 0.8390830|
|conditionKOEX                                     |  0.1025831| 0.8575237|  0.1196271| 0.9054088|
|metabolite3HMG                                    | -0.0535035| 0.8420109| -0.0635426| 0.9493509|
|metaboliteacetylcarnitine                         | -0.0535035| 0.8420109| -0.0635426| 0.9493509|
|metabolitebutyrylcarnitine                        | -0.0535035| 0.8420109| -0.0635426| 0.9493509|
|metaboliteC10:1 total                             | -0.0146475| 0.8618324| -0.0169957| 0.9864444|
|metaboliteC12                                     | -0.0535035| 0.8420109| -0.0635426| 0.9493509|
|metaboliteC14:1 total                             | -0.0535035| 0.8420109| -0.0635426| 0.9493509|
|metaboliteC14:2 total                             | -0.0535035| 0.8420109| -0.0635426| 0.9493509|
|metaboliteC15:1 total                             | -0.0535035| 0.8420109| -0.0635426| 0.9493509|
|metaboliteC16:1 total                             | -0.0830624| 0.8635372| -0.0961885| 0.9233958|
|metaboliteC17:1 total                             | -0.0535035| 0.8420109| -0.0635426| 0.9493509|
|metaboliteC18:1 total                             | -0.0830624| 0.8635372| -0.0961885| 0.9233958|
|metaboliteC18:2 total                             | -0.0830624| 0.8635372| -0.0961885| 0.9233958|
|metaboliteC19:1 total                             | -0.0830624| 0.8635372| -0.0961885| 0.9233958|
|metabolitecarnitine                               | -0.0535035| 0.8420109| -0.0635426| 0.9493509|
|metaboliteethylmalonylcarnitine                   | -0.0535035| 0.8420109| -0.0635426| 0.9493509|
|metaboliteisobutyrylcarnitine                     | -0.0448009| 0.8618324| -0.0519833| 0.9585554|
|metabolitemethylsuccinylcarnitine                 | -0.0535035| 0.8420109| -0.0635426| 0.9493509|
|metaboliten-decanoylcarnitine                     | -0.0535035| 0.8420109| -0.0635426| 0.9493509|
|metaboliten-hexanoylcarnitine                     | -0.0535035| 0.8420109| -0.0635426| 0.9493509|
|metaboliten-octanoylcarnitine                     | -0.0535035| 0.8420109| -0.0635426| 0.9493509|
|metabolitepropionylcarnitine                      | -0.0535035| 0.8420109| -0.0635426| 0.9493509|
|conditionWTEX:metabolite3HMG                      | -2.7317822| 1.1313121| -2.4147025| 0.0159793|
|conditionKOREST:metabolite3HMG                    | -0.8785031| 1.2286967| -0.7149878| 0.4748325|
|conditionKOEX:metabolite3HMG                      | -0.9703957| 1.1313654| -0.8577208| 0.3913125|
|conditionWTEX:metaboliteacetylcarnitine           | -2.5752214| 1.1313121| -2.2763138| 0.0230999|
|conditionKOREST:metaboliteacetylcarnitine         | -1.6021606| 1.2148851| -1.3187754| 0.1876348|
|conditionKOEX:metaboliteacetylcarnitine           | -2.3445605| 1.1313734| -2.0723135| 0.0385672|
|conditionWTEX:metabolitebutyrylcarnitine          | -2.9898929| 1.1313121| -2.6428541| 0.0083875|
|conditionKOREST:metabolitebutyrylcarnitine        | -1.7395579| 1.2148851| -1.4318703| 0.1525852|
|conditionKOEX:metabolitebutyrylcarnitine          | -3.1366307| 1.1313224| -2.7725349| 0.0056960|
|conditionWTEX:metaboliteC10:1 total               | -1.9843514| 1.1589086| -1.7122587| 0.0872501|
|conditionKOREST:metaboliteC10:1 total             | -2.2856069| 1.2287061| -1.8601739| 0.0632405|
|conditionKOEX:metaboliteC10:1 total               | -1.8885233| 1.1467628| -1.6468299| 0.0999995|
|conditionWTEX:metaboliteC12                       | -2.4999123| 1.1313121| -2.2097459| 0.0274152|
|conditionKOREST:metaboliteC12                     | -0.5874357| 1.2286769| -0.4781043| 0.6327111|
|conditionKOEX:metaboliteC12                       | -0.4150223| 1.1313263| -0.3668458| 0.7138344|
|conditionWTEX:metaboliteC14:1 total               | -2.6195379| 1.1313121| -2.3154865| 0.0208473|
|conditionKOREST:metaboliteC14:1 total             | -0.7352843| 1.2286769| -0.5984359| 0.5497244|
|conditionKOEX:metaboliteC14:1 total               | -0.5991443| 1.1313729| -0.5295727| 0.5965603|
|conditionWTEX:metaboliteC14:2 total               | -2.8482016| 1.1313121| -2.5176091| 0.0120165|
|conditionKOREST:metaboliteC14:2 total             | -0.6077469| 1.2148851| -0.5002505| 0.6170411|
|conditionKOEX:metaboliteC14:2 total               | -1.0557320| 1.1313249| -0.9331821| 0.3510174|
|conditionWTEX:metaboliteC15:1 total               | -2.7379024| 1.1313121| -2.4201123| 0.0157453|
|conditionKOREST:metaboliteC15:1 total             | -0.1911423| 1.2148851| -0.1573336| 0.8750230|
|conditionKOEX:metaboliteC15:1 total               | -0.6493711| 1.1318551| -0.5737229| 0.5663222|
|conditionWTEX:metaboliteC16:1 total               | -3.2521798| 1.1474238| -2.8343318| 0.0047119|
|conditionKOREST:metaboliteC16:1 total             |  4.8826656| 1.2299025|  3.9699615| 0.0000786|
|conditionKOEX:metaboliteC16:1 total               |  5.4429400| 1.1480122|  4.7411865| 0.0000025|
|conditionWTEX:metaboliteC17:1 total               | -2.9138252| 1.1313121| -2.5756156| 0.0101909|
|conditionKOREST:metaboliteC17:1 total             |  1.1834289| 1.2148851|  0.9741077| 0.3303078|
|conditionKOEX:metaboliteC17:1 total               |  0.6643010| 1.1313786|  0.5871607| 0.5572673|
|conditionWTEX:metaboliteC18:1 total               | -2.8945636| 1.1474238| -2.5226631| 0.0118466|
|conditionKOREST:metaboliteC18:1 total             |  7.6409522| 1.2299025|  6.2126488| 0.0000000|
|conditionKOEX:metaboliteC18:1 total               |  8.5666217| 1.1477120|  7.4640865| 0.0000000|
|conditionWTEX:metaboliteC18:2 total               | -3.6540268| 1.1474238| -3.1845486| 0.0015079|
|conditionKOREST:metaboliteC18:2 total             |  4.1186451| 1.2299025|  3.3487574| 0.0008511|
|conditionKOEX:metaboliteC18:2 total               |  2.3147734| 1.1474923|  2.0172452| 0.0440151|
|conditionWTEX:metaboliteC19:1 total               | -2.4174635| 1.1474238| -2.1068619| 0.0354512|
|conditionKOREST:metaboliteC19:1 total             | 21.3123445| 1.2299025| 17.3284831| 0.0000000|
|conditionKOEX:metaboliteC19:1 total               | 20.4516089| 1.1475001| 17.8227518| 0.0000000|
|conditionWTEX:metabolitecarnitine                 | -0.9540640| 1.1313121| -0.8433252| 0.3993075|
|conditionKOREST:metabolitecarnitine               | -1.7670066| 1.2148851| -1.4544639| 0.1462236|
|conditionKOEX:metabolitecarnitine                 | -2.2484981| 1.1450274| -1.9637068| 0.0499221|
|conditionWTEX:metaboliteethylmalonylcarnitine     |  0.3773822| 1.1313121|  0.3335792| 0.7387876|
|conditionKOREST:metaboliteethylmalonylcarnitine   | -2.4686276| 1.2148851| -2.0319844| 0.0424967|
|conditionKOEX:metaboliteethylmalonylcarnitine     | -2.4290124| 1.1313393| -2.1470239| 0.0321018|
|conditionWTEX:metaboliteisobutyrylcarnitine       | -1.3217549| 1.1589086| -1.1405170| 0.2544242|
|conditionKOREST:metaboliteisobutyrylcarnitine     |  0.2081126| 1.2805701|  0.1625156| 0.8709423|
|conditionKOEX:metaboliteisobutyrylcarnitine       |  1.8658907| 1.1461619|  1.6279469| 0.1039436|
|conditionWTEX:metabolitemethylsuccinylcarnitine   |  1.2796354| 1.1313121|  1.1311073| 0.2583607|
|conditionKOREST:metabolitemethylsuccinylcarnitine | -2.3155150| 1.2148851| -1.9059539| 0.0570268|
|conditionKOEX:metabolitemethylsuccinylcarnitine   | -2.2064553| 1.1458319| -1.9256362| 0.0545163|
|conditionWTEX:metaboliten-decanoylcarnitine       | -3.1098466| 1.1313121| -2.7488847| 0.0061194|
|conditionKOREST:metaboliten-decanoylcarnitine     | -0.9443264| 1.2148851| -0.7772968| 0.4372215|
|conditionKOEX:metaboliten-decanoylcarnitine       | -1.1731607| 1.1448904| -1.0246926| 0.3058288|
|conditionWTEX:metaboliten-hexanoylcarnitine       | -2.2713708| 1.1313121| -2.0077314| 0.0450194|
|conditionKOREST:metaboliten-hexanoylcarnitine     | -1.3210105| 1.2148851| -1.0873542| 0.2772192|
|conditionKOEX:metaboliten-hexanoylcarnitine       | -2.2367430| 1.1313436| -1.9770678| 0.0483889|
|conditionWTEX:metaboliten-octanoylcarnitine       | -1.7134011| 1.1313121| -1.5145256| 0.1303014|
|conditionKOREST:metaboliten-octanoylcarnitine     | -1.4925508| 1.2300129| -1.2134431| 0.2253311|
|conditionKOEX:metaboliten-octanoylcarnitine       | -1.7021729| 1.1442886| -1.4875380| 0.1372806|
|conditionWTEX:metabolitepropionylcarnitine        | -0.5605228| 1.1313121| -0.4954625| 0.6204145|
|conditionKOREST:metabolitepropionylcarnitine      | -2.1220764| 1.2148851| -1.7467300| 0.0810813|
|conditionKOEX:metabolitepropionylcarnitine        | -2.3785860| 1.1444573| -2.0783529| 0.0380062|

```r
M %>% tidy(effects = "fixed") %>% write.csv(file = "data/processed/lmeFixedCoefAim1b.csv", row.names = FALSE)
M1b <- M
```


## Aim 2

* WT-exer (ref)
* KO-C7
* KO-C8
* KO-exer

### Amino Acids and TCA Cycle intermediates (Neg)


```r
t0 <- Sys.time()
M <- estimateModel(data = D2a)
Sys.time() - t0
```

```
## Time difference of 1.78132 secs
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



|                     | numDF| denDF|  F-value| p-value|
|:--------------------|-----:|-----:|--------:|-------:|
|(Intercept)          |     1|   899| 42.19271|       0|
|condition            |     3|    40| 24.34463|       0|
|metabolite           |    23|   899| 67.07559|       0|
|condition:metabolite |    69|   899| 52.36200|       0|

```r
M %>% tidy(effects = "fixed") %>% kable()
```



|term                                   |   estimate| std.error|  statistic|   p.value|
|:--------------------------------------|----------:|---------:|----------:|---------:|
|(Intercept)                            |  0.0039017| 0.6667134|  0.0058522| 0.9953320|
|conditionKOEX                          |  0.6157507| 0.9170512|  0.6714463| 0.5057957|
|conditionKOC7                          | -0.1373235| 0.9432649| -0.1455832| 0.8849816|
|conditionKOC8                          |  1.3409355| 0.9654053|  1.3889871| 0.1725220|
|metabolitearginine                     | -0.0111821| 0.9030235| -0.0123830| 0.9901228|
|metaboliteasparagine                   | -0.0040280| 0.9030235| -0.0044606| 0.9964420|
|metaboliteAspartic Acid                |  0.0280093| 0.9258227|  0.0302534| 0.9758717|
|metaboliteCITRIC                       | -0.0083106| 0.9030235| -0.0092031| 0.9926591|
|metaboliteFUMARIC                      | -0.0039017| 0.9051566| -0.0043105| 0.9965617|
|metaboliteglutamic                     | -0.0044970| 0.9030235| -0.0049799| 0.9960277|
|metaboliteglutamine                    | -0.0039777| 0.9030235| -0.0044048| 0.9964864|
|metaboliteglycine                      | -0.0090560| 0.9030235| -0.0100285| 0.9920007|
|metabolitehistidine                    | -0.0108380| 0.9030235| -0.0120019| 0.9904268|
|metaboliteisoleucine                   | -0.0132796| 0.9030235| -0.0147058| 0.9882702|
|metaboliteLACTIC                       | -0.0039017| 0.9051566| -0.0043105| 0.9965617|
|metaboliteleucine                      | -0.0147016| 0.9030235| -0.0162804| 0.9870143|
|metabolitelysine                       |  0.0016179| 0.9030235|  0.0017916| 0.9985709|
|metaboliteMALIC                        | -0.0039017| 0.9051566| -0.0043105| 0.9965617|
|metabolitemethionine                   |  0.0042381| 0.9030235|  0.0046933| 0.9962564|
|metaboliteMETHYLSUCCINIC               | -0.0080403| 0.9275740| -0.0086681| 0.9930859|
|metabolitePhenylalanine                | -0.0070516| 0.9258227| -0.0076166| 0.9939246|
|metaboliteserine                       | -0.0072949| 0.9030235| -0.0080783| 0.9935563|
|metaboliteSUCCINIC-2                   | -0.0039017| 0.9051566| -0.0043105| 0.9965617|
|metabolitethreonine                    |  0.0307575| 0.9258227|  0.0332218| 0.9735051|
|metaboliteTryptophan                   | -0.0067022| 0.9030235| -0.0074219| 0.9940799|
|metaboliteTyrosine                     | -0.0039017| 0.9051566| -0.0043105| 0.9965617|
|metabolitevaline                       | -0.0078125| 0.9030235| -0.0086515| 0.9930991|
|conditionKOEX:metabolitearginine       |  0.1510645| 1.2292556|  0.1228911| 0.9022208|
|conditionKOC7:metabolitearginine       |  1.1667572| 1.2788647|  0.9123382| 0.3618353|
|conditionKOC8:metabolitearginine       | -0.8567436| 1.2952813| -0.6614343| 0.5085033|
|conditionKOEX:metaboliteasparagine     |  0.3126286| 1.2294393|  0.2542855| 0.7993331|
|conditionKOC7:metaboliteasparagine     |  1.1734040| 1.2788647|  0.9175357| 0.3591081|
|conditionKOC8:metaboliteasparagine     |  0.6942373| 1.2952813|  0.5359742| 0.5921090|
|conditionKOEX:metaboliteAspartic Acid  |  0.8179432| 1.2395063|  0.6598943| 0.5094907|
|conditionKOC7:metaboliteAspartic Acid  |  3.1067722| 1.2950642|  2.3989329| 0.0166457|
|conditionKOC8:metaboliteAspartic Acid  |  2.4400578| 1.3112780|  1.8608241| 0.0630953|
|conditionKOEX:metaboliteCITRIC         | -2.2192054| 1.2431849| -1.7850969| 0.0745827|
|conditionKOC7:metaboliteCITRIC         | -1.1411850| 1.2953303| -0.8809992| 0.3785537|
|conditionKOC8:metaboliteCITRIC         | -2.2388154| 1.2952813| -1.7284395| 0.0842528|
|conditionKOEX:metaboliteFUMARIC        | -1.6963514| 1.2403226| -1.3676695| 0.1717573|
|conditionKOC7:metaboliteFUMARIC        | -2.6028011| 1.3165641| -1.9769650| 0.0483508|
|conditionKOC8:metaboliteFUMARIC        | -2.6375614| 1.2967694| -2.0339479| 0.0422500|
|conditionKOEX:metaboliteglutamic       | -1.4732017| 1.2185274| -1.2090017| 0.2269801|
|conditionKOC7:metaboliteglutamic       |  0.0379106| 1.2788647|  0.0296440| 0.9763576|
|conditionKOC8:metaboliteglutamic       | -1.0841855| 1.2952813| -0.8370270| 0.4027998|
|conditionKOEX:metaboliteglutamine      | -0.4424009| 1.2303233| -0.3595810| 0.7192449|
|conditionKOC7:metaboliteglutamine      | -0.5883785| 1.2788647| -0.4600788| 0.6455709|
|conditionKOC8:metaboliteglutamine      | -0.4039481| 1.2952813| -0.3118613| 0.7552182|
|conditionKOEX:metaboliteglycine        | -0.1391318| 1.2321663| -0.1129164| 0.9101220|
|conditionKOC7:metaboliteglycine        |  0.2707250| 1.2788647|  0.2116917| 0.8323956|
|conditionKOC8:metaboliteglycine        | -0.0331649| 1.3105340| -0.0253064| 0.9798162|
|conditionKOEX:metabolitehistidine      | -0.5693372| 1.2372282| -0.4601716| 0.6455044|
|conditionKOC7:metabolitehistidine      |  0.6151717| 1.2788647|  0.4810295| 0.6306126|
|conditionKOC8:metabolitehistidine      | -1.7352620| 1.2952813| -1.3396796| 0.1806881|
|conditionKOEX:metaboliteisoleucine     | -0.5968069| 1.2328301| -0.4840950| 0.6284363|
|conditionKOC7:metaboliteisoleucine     | -0.6993975| 1.2788647| -0.5468893| 0.5845905|
|conditionKOC8:metaboliteisoleucine     | -2.3994282| 1.2952813| -1.8524378| 0.0642906|
|conditionKOEX:metaboliteLACTIC         | -2.4744865| 1.2197312| -2.0287146| 0.0427817|
|conditionKOC7:metaboliteLACTIC         | -2.7553836| 1.2968183| -2.1247261| 0.0338818|
|conditionKOC8:metaboliteLACTIC         | -5.3043789| 1.2967694| -4.0904565| 0.0000469|
|conditionKOEX:metaboliteleucine        |  0.4329293| 1.2483439|  0.3468029| 0.7288205|
|conditionKOC7:metaboliteleucine        | -0.4371929| 1.2788647| -0.3418601| 0.7325361|
|conditionKOC8:metaboliteleucine        | -1.7756742| 1.2952813| -1.3708792| 0.1707548|
|conditionKOEX:metabolitelysine         | -0.3950677| 1.2392184| -0.3188039| 0.7499493|
|conditionKOC7:metabolitelysine         | -0.0251930| 1.2788647| -0.0196995| 0.9842875|
|conditionKOC8:metabolitelysine         | -0.5465025| 1.2952813| -0.4219180| 0.6731858|
|conditionKOEX:metaboliteMALIC          | -2.5223717| 1.2450768| -2.0258765| 0.0430724|
|conditionKOC7:metaboliteMALIC          | -1.6855378| 1.2968183| -1.2997486| 0.1940203|
|conditionKOC8:metaboliteMALIC          | -3.8260606| 1.2967694| -2.9504556| 0.0032554|
|conditionKOEX:metabolitemethionine     | -0.1485534| 1.2246510| -0.1213027| 0.9034784|
|conditionKOC7:metabolitemethionine     |  1.3735703| 1.2788647|  1.0740544| 0.2830865|
|conditionKOC8:metabolitemethionine     |  0.6129785| 1.2952813|  0.4732396| 0.6361571|
|conditionKOEX:metaboliteMETHYLSUCCINIC | -1.1821646| 1.2509799| -0.9449909| 0.3449174|
|conditionKOC7:metaboliteMETHYLSUCCINIC | -0.9951659| 1.3320759| -0.7470790| 0.4552112|
|conditionKOC8:metaboliteMETHYLSUCCINIC | 42.2293977| 1.3125151| 32.1744086| 0.0000000|
|conditionKOEX:metabolitePhenylalanine  | -0.0994187| 1.2623800| -0.0787550| 0.9372450|
|conditionKOC7:metabolitePhenylalanine  | -0.3700238| 1.2950642| -0.2857185| 0.7751595|
|conditionKOC8:metabolitePhenylalanine  | -2.1422092| 1.3112780| -1.6336804| 0.1026761|
|conditionKOEX:metaboliteserine         |  0.4207857| 1.2293425|  0.3422852| 0.7322163|
|conditionKOC7:metaboliteserine         |  1.9370902| 1.2788647|  1.5146951| 0.1302011|
|conditionKOC8:metaboliteserine         |  3.2761774| 1.3105340|  2.4998797| 0.0126008|
|conditionKOEX:metaboliteSUCCINIC-2     | -1.8164653| 1.2422803| -1.4622025| 0.1440352|
|conditionKOC7:metaboliteSUCCINIC-2     |  0.0398480| 1.2968183|  0.0307275| 0.9754937|
|conditionKOC8:metaboliteSUCCINIC-2     | -1.6263740| 1.2967694| -1.2541737| 0.2101049|
|conditionKOEX:metabolitethreonine      |  3.0136288| 1.2386226|  2.4330485| 0.0151663|
|conditionKOC7:metabolitethreonine      | 10.9387180| 1.2950642|  8.4464674| 0.0000000|
|conditionKOC8:metabolitethreonine      |  7.8694065| 1.3112780|  6.0013256| 0.0000000|
|conditionKOEX:metaboliteTryptophan     | -0.2361494| 1.2256110| -0.1926789| 0.8472539|
|conditionKOC7:metaboliteTryptophan     | -0.1074139| 1.2788647| -0.0839916| 0.9330818|
|conditionKOC8:metaboliteTryptophan     | -3.4267670| 1.3115636| -2.6127341| 0.0091321|
|conditionKOEX:metaboliteTyrosine       | -0.0091659| 1.2463057| -0.0073545| 0.9941337|
|conditionKOC7:metaboliteTyrosine       | -0.1171141| 1.2803719| -0.0914688| 0.9271405|
|conditionKOC8:metaboliteTyrosine       | -1.6966254| 1.2967694| -1.3083478| 0.1910898|
|conditionKOEX:metabolitevaline         |  0.4272935| 1.2294024|  0.3475620| 0.7282505|
|conditionKOC7:metabolitevaline         | -0.5447894| 1.2788647| -0.4259945| 0.6702138|
|conditionKOC8:metabolitevaline         |  0.3923559| 1.2952813|  0.3029118| 0.7620272|

```r
M %>% tidy(effects = "fixed") %>% write.csv(file = "data/processed/lmeFixedCoefAim2a.csv", row.names = FALSE)
M2a <- M
```

### Acylcarnitines


```r
t0 <- Sys.time()
M <- estimateModel(data = D2b)
Sys.time() - t0
```

```
## Time difference of 0.8821671 secs
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



|                     | numDF| denDF|  F-value| p-value|
|:--------------------|-----:|-----:|--------:|-------:|
|(Intercept)          |     1|   820| 168.6622|       0|
|condition            |     3|    40| 108.6826|       0|
|metabolite           |    21|   820| 403.8874|       0|
|condition:metabolite |    63|   820| 165.8649|       0|

```r
M %>% tidy(effects = "fixed") %>% kable()
```



|term                                            |   estimate| std.error|  statistic|   p.value|
|:-----------------------------------------------|----------:|---------:|----------:|---------:|
|(Intercept)                                     |  0.0000000| 0.7719295|  0.0000000| 1.0000000|
|conditionKOEX                                   | -1.0241094| 1.0906711| -0.9389717| 0.3533809|
|conditionKOC7                                   |  0.0175332| 1.1167773|  0.0156998| 0.9875520|
|conditionKOC8                                   | -1.2156145| 1.0916731| -1.1135334| 0.2721262|
|metabolite3HMG                                  |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metaboliteacetylcarnitine                       |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metabolitebutyrylcarnitine                      |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metaboliteC10:1 total                           |  0.0015254| 1.0622811|  0.0014360| 0.9988546|
|metaboliteC12                                   |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metaboliteC14:1 total                           |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metaboliteC14:2 total                           |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metaboliteC15:1 total                           |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metaboliteC16:1 total                           |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metaboliteC17:1 total                           |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metaboliteC18:1 total                           |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metaboliteC18:2 total                           |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metaboliteC19:1 total                           |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metabolitecarnitine                             |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metaboliteethylmalonylcarnitine                 |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metaboliteisobutyrylcarnitine                   |  0.0245746| 1.0622811|  0.0231338| 0.9815492|
|metabolitemethylsuccinylcarnitine               |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metaboliten-decanoylcarnitine                   |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metaboliten-hexanoylcarnitine                   |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metaboliten-octanoylcarnitine                   |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|metabolitepropionylcarnitine                    |  0.0000000| 1.0358591|  0.0000000| 1.0000000|
|conditionKOEX:metabolite3HMG                    |  0.2722061| 1.4614212|  0.1862612| 0.8522860|
|conditionKOC7:metabolite3HMG                    | -0.8632232| 1.5024152| -0.5745570| 0.5657485|
|conditionKOC8:metabolite3HMG                    |  0.7305901| 1.4649260|  0.4987215| 0.6181093|
|conditionKOEX:metaboliteacetylcarnitine         | -2.1680276| 1.4623284| -1.4825859| 0.1385687|
|conditionKOC7:metaboliteacetylcarnitine         | -2.4879347| 1.4837283| -1.6768129| 0.0939602|
|conditionKOC8:metaboliteacetylcarnitine         | -0.8712095| 1.4649260| -0.5947123| 0.5521998|
|conditionKOEX:metabolitebutyrylcarnitine        | -1.4711570| 1.4633385| -1.0053429| 0.3150283|
|conditionKOC7:metabolitebutyrylcarnitine        | -1.7130344| 1.5025598| -1.1400774| 0.2545870|
|conditionKOC8:metabolitebutyrylcarnitine        |  0.7150015| 1.4649260|  0.4880803| 0.6256235|
|conditionKOEX:metaboliteC10:1 total             | -0.0046739| 1.4815763| -0.0031547| 0.9974837|
|conditionKOC7:metaboliteC10:1 total             | -0.9274611| 1.5022938| -0.6173633| 0.5371665|
|conditionKOC8:metaboliteC10:1 total             |  0.4905384| 1.5023583|  0.3265123| 0.7441201|
|conditionKOEX:metaboliteC12                     |  0.8616196| 1.4633027|  0.5888184| 0.5561453|
|conditionKOC7:metaboliteC12                     |  2.6557825| 1.5025598|  1.7675054| 0.0775155|
|conditionKOC8:metaboliteC12                     |  1.4133870| 1.4649260|  0.9648180| 0.3349204|
|conditionKOEX:metaboliteC14:1 total             |  0.8501425| 1.4602261|  0.5821992| 0.5605925|
|conditionKOC7:metaboliteC14:1 total             |  3.2722637| 1.4837283|  2.2054332| 0.0277005|
|conditionKOC8:metaboliteC14:1 total             |  0.9057273| 1.4837936|  0.6104132| 0.5417572|
|conditionKOEX:metaboliteC14:2 total             |  0.4068890| 1.4607765|  0.2785430| 0.7806658|
|conditionKOC7:metaboliteC14:2 total             |  1.8982474| 1.4837283|  1.2793767| 0.2011263|
|conditionKOC8:metaboliteC14:2 total             |  0.7198717| 1.4649260|  0.4914048| 0.6232717|
|conditionKOEX:metaboliteC15:1 total             |  0.8152323| 1.4653891|  0.5563248| 0.5781406|
|conditionKOC7:metaboliteC15:1 total             |  9.2087187| 1.4837283|  6.2064724| 0.0000000|
|conditionKOC8:metaboliteC15:1 total             |  0.9030802| 1.4649260|  0.6164681| 0.5377567|
|conditionKOEX:metaboliteC16:1 total             |  5.4602425| 1.4661478|  3.7242100| 0.0002093|
|conditionKOC7:metaboliteC16:1 total             |  5.9847363| 1.4837283|  4.0335797| 0.0000601|
|conditionKOC8:metaboliteC16:1 total             |  6.0671337| 1.4649260|  4.1415973| 0.0000381|
|conditionKOEX:metaboliteC17:1 total             |  4.1843188| 1.4600554|  2.8658630| 0.0042652|
|conditionKOC7:metaboliteC17:1 total             | 26.9662901| 1.4837283| 18.1746822| 0.0000000|
|conditionKOC8:metaboliteC17:1 total             |  3.2429032| 1.4649260|  2.2136976| 0.0271241|
|conditionKOEX:metaboliteC18:1 total             |  6.0524906| 1.4607084|  4.1435311| 0.0000378|
|conditionKOC7:metaboliteC18:1 total             |  2.3119974| 1.4837283|  1.5582350| 0.1195634|
|conditionKOC8:metaboliteC18:1 total             |  6.8165550| 1.4649260|  4.6531736| 0.0000038|
|conditionKOEX:metaboliteC18:2 total             |  6.9787991| 1.4603540|  4.7788407| 0.0000021|
|conditionKOC7:metaboliteC18:2 total             |  6.5680548| 1.4837283|  4.4267235| 0.0000109|
|conditionKOC8:metaboliteC18:2 total             |  6.2882413| 1.4649260|  4.2925316| 0.0000198|
|conditionKOEX:metaboliteC19:1 total             | 20.6873791| 1.4609942| 14.1597951| 0.0000000|
|conditionKOC7:metaboliteC19:1 total             | 95.3148684| 1.4837283| 64.2401098| 0.0000000|
|conditionKOC8:metaboliteC19:1 total             | 16.4486682| 1.4649260| 11.2283269| 0.0000000|
|conditionKOEX:metabolitecarnitine               | -2.6884543| 1.4787881| -1.8180119| 0.0694271|
|conditionKOC7:metabolitecarnitine               | -3.2607741| 1.4837283| -2.1976895| 0.0282501|
|conditionKOC8:metabolitecarnitine               | -2.0852062| 1.4649260| -1.4234208| 0.1549946|
|conditionKOEX:metaboliteethylmalonylcarnitine   | -4.0178874| 1.4629153| -2.7464936| 0.0061556|
|conditionKOC7:metaboliteethylmalonylcarnitine   | -4.8400656| 1.5025598| -3.2212133| 0.0013268|
|conditionKOC8:metaboliteethylmalonylcarnitine   | -3.7320303| 1.4837275| -2.5153071| 0.0120830|
|conditionKOEX:metaboliteisobutyrylcarnitine     |  1.4103948| 1.4819445|  0.9517191| 0.3415199|
|conditionKOC7:metaboliteisobutyrylcarnitine     |  0.6011983| 1.5208955|  0.3952923| 0.6927300|
|conditionKOC8:metaboliteisobutyrylcarnitine     |  1.6615857| 1.4837268|  1.1198731| 0.2630958|
|conditionKOEX:metabolitemethylsuccinylcarnitine | -3.4524129| 1.4786767| -2.3347989| 0.0197944|
|conditionKOC7:metabolitemethylsuccinylcarnitine | -4.3332826| 1.5025598| -2.8839335| 0.0040302|
|conditionKOC8:metabolitemethylsuccinylcarnitine | -3.0804273| 1.4837936| -2.0760483| 0.0382005|
|conditionKOEX:metaboliten-decanoylcarnitine     |  0.6425131| 1.4835969|  0.4330780| 0.6650720|
|conditionKOC7:metaboliten-decanoylcarnitine     |  2.7123558| 1.4837283|  1.8280677| 0.0679026|
|conditionKOC8:metaboliten-decanoylcarnitine     |  4.4792402| 1.4837936|  3.0187757| 0.0026166|
|conditionKOEX:metaboliten-hexanoylcarnitine     | -0.4795217| 1.4621964| -0.3279461| 0.7430361|
|conditionKOC7:metaboliten-hexanoylcarnitine     | -1.2785130| 1.4837283| -0.8616894| 0.3891104|
|conditionKOC8:metaboliten-hexanoylcarnitine     |  0.6568021| 1.4649260|  0.4483517| 0.6540178|
|conditionKOEX:metaboliten-octanoylcarnitine     | -0.5328835| 1.4842702| -0.3590206| 0.7196721|
|conditionKOC7:metaboliten-octanoylcarnitine     | -1.1292802| 1.4837283| -0.7611099| 0.4468103|
|conditionKOC8:metaboliten-octanoylcarnitine     |  1.6013774| 1.4649260|  1.0931456| 0.2746509|
|conditionKOEX:metabolitepropionylcarnitine      | -2.1389622| 1.4807886| -1.4444751| 0.1489874|
|conditionKOC7:metabolitepropionylcarnitine      | -2.3577460| 1.5025598| -1.5691528| 0.1169981|
|conditionKOC8:metabolitepropionylcarnitine      | -1.9888791| 1.4649260| -1.3576652| 0.1749434|

```r
M %>% tidy(effects = "fixed") %>% write.csv(file = "data/processed/lmeFixedCoefAim2b.csv", row.names = FALSE)
M2b <- M
```


## Save

Save `lme` objects for interactive work.


```r
save(M1a, M1b, M2a, M2b, file = "data/processed/lmeObjects.RData")
```
