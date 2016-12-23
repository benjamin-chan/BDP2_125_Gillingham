# Metabolomics of very long-chain aclCoA dehydrogenase knockout mice

Attribute | Value
---|---
Principal investigator | Melanie Gillingham, gillingm@ohsu.edu, 503-494-1682
Main contact | Melanie Gillingham, gillingm@ohsu.edu, 503-494-1682
Statistician | Benjamin Chan, chanb@ohsu.edu, 503-494-5491


## Objective

To compare wild type versus very long-chain acylCoA dehydrogenase knockout mice on the following experimental conditions

A. Rest versus exhaustive exercise  
B. MCT chow versus experimental chow with triheptanoin  

The design of both experiments is a two-factor balanced design with multiple measurements (300 metabolites per mouse). 
The primary challenge of the analysis will be to reduce the data (using principle components analysis or other data reduction techniques) or to account for the correlation structure of the multiple measurements.


## Deliverables

1. Analysis of experiment A
  * Standardize measurements
  * Perform data reduction and model reduced data
  * Model non-reduced data accounting for within-mouse correlation structure
  * Compare models
  * Summarize
2. Analysis of experiment B
  * Same as above
3. Text describing the statistical methods used, suitable for a manuscript


## Interpretation of analysis

* Aim 1, activity: exercise vs rest
  * None of the log-transformed metabolite values in either the wildtype or
    knockout mice were statistically significant at $\alpha$-level 0.05
    (controlling for the false discovery rate)
* Aim 2, chow: yellow C8 vs white C7
  * **3-hydroxybutyric** values were significantly **higher** in wildtype mice feed yellow C8 compared to wildtype mice fed white C7 ($\beta$ = 0.81, p = 0.007).
    The large estimated effect is likely due to an outlying value among the wilftype mice fed white C7.
  * **LC odd AC total** values were significantly **lower** in wildtype mice feed yellow C8 compared to wildtype mice fed white C7 ($\beta$ = -0.70, p = 0.040)
  * **Methylsuccinic** values were significantly **higher** in wildtype mice feed yellow C8 compared to wildtype mice fed white C7 ($\beta$ = 0.74, p = 0.040)


## Outputs

* Results: [Excel XLSX](data/processed/contrasts.xlsx)
  * Formatted versions of `data/processed/contrastsAim1.csv` and `data/processed/contrastsAim2.csv`
* Aim 1
  * Figure, metabolite values: [PNG](figures/plotDataAim1.png)
  * Contrast estimates: [CSV](data/processed/contrastsAim1.csv)
* Aim 2
  * Figure, metabolite values: [PNG](figures/plotDataAim2.png)
  * Contrast estimates: [CSV](data/processed/contrastsAim2.csv)
* Complete analysis: [HTML](docs/index.html) or [Markdown](docs/index.md)


## Methods

A mixed linear effects model was estimated for each aim.
Metabolomic values were transformed to the log-10 scale.
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
