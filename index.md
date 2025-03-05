---
layout: page
title: Exact matching and matching-adjusted indirect comparisons
tagline: 
description: Exact matching and matching-adjusted indirect comparisons
---

<!-- badges: start -->

[![CRAN status](https://badges.cranchecks.info/flavor/release/maicChecks.svg)](https://cran.r-project.org/web/checks/check_results_maicChecks.html) [![downloads](https://cranlogs.r-pkg.org/badges/maicChecks)](https://www.rdocumentation.org/trends)<a href="https://github.com/clyau/maicChecks"><img src="MaicChecks.png" align="right"/></a>

<!-- badges: end -->


## Installation

|Package        | Type                 | Source      | Command                          |
|:--------------|:---------------------|:------------|:---------------------------------|
|`maicChecks`   | Release (0.2.0)      | CRAN        | `install.packages("maicChecks")` |

## Version history

**Version 0.1.2** was released on CRAN on 4 January, 2022. It implements the following methods:

-   Matching-adjusted indirect comparison (MAIC) for comparing study with patient level data to study with only aggregated data ([Signorovitch (2010)](#reference))    
-   Exact matching by maximizing effective sample size (ESS) for comparing study with patient level data to study with only aggregated data ([Glimm and Yau (2022)](#reference))     
-   Various checks to assess feasibility of conducting MAIC as introduced in ([Glimm and Yau (2022)](#reference))

**Version 0.2.0** was released on CRAN on 3 March, 2025. It implements the following methods:

-   Exact matching by maximizing effective sample size (ESS) for comparing two studies both with patient level data ([Glimm and Yau (2025)](#reference))

## Overview

The comparison of different medical treatments from observational studies or across different clinical studies is often biased by confounding factors such as systematic differences in patient demographics or in the inclusion criteria for the trials. The confounding must be adjusted before indirect comparisons can be conducted. The adjustment is usually accomplished by matching the baseline covariates so that patients from one or both studies are each assigned a weight, which are taken into account when comparing clinical outcomes.

From a data availability prespective, two situations arise:

1.  Individual patient-level data (IPD) are available for both studies, i.e., an IPD vs IPD comparison.
2.  IPD are available for one study, and only aggregated data (AD), i.e., summary statistics, are available for the other, i.e., and IPD vs AD comparison.

In both situations, one study can be considered a "target" population, and the other is then matched onto it. However, in the first situation when IPD are available in both, the two study populations can also be matched onto a common pooled population. 

`maicChecks` is an R package that offers two different but related methods for matching baseline covariates for these two situations. The methods are:

1.  Exact matching: used for IPD vs IPD and IPD vs AD 
2.  Matching-adjusted indirect comparison: used for IPD vs AD

## Exact matching

Exact matching should be used when IPD are available from both studies. It is an alternative to propensity score matching. The method ensures that after matching, the weighted means of the baseline covariates between the two studies are exactly the same. Details on methodology can be found in [Glimm & Yau (2025)](#reference).

### Methods

#### Linear check

The function `maicChecks::exmLP.2ipd()` checks if there is an overlap between the two IPD. If yes, matching of baseline covariates can be performed. In practice, there is almost always an overlap between the two IPD. 

#### Exact matching

The function `maicChecks::exmWt.2ipd()` matches the baseline covariates from the two IPD, and assigns a weight to each patient in the two studies. The algorithm treats matching as a constrained optimization problem Constrained optimization is a purely algebraic technique and solves the convex optimization problem in a finite number of steps. In other words, the method does not require numerical approximation.

Although the weighted means are the same for the two studies after matching, it can happen that in one of the covariates used in matching it is not between the two observed means. To avoid this, an additional constraint can be added to force the weighted means to be always between the observed means. Naturally, with an additional constraint, the likelihood of a non-existing solution increases.

(Note that `maicChecks::exmLP.2ipd()` which checks for feasibility of a match uses linear programming. Hence, it can happen that the check thinks it's fine to go ahead with the matching, but `maicChecks::exmWt.2ipd()` yields no solution.)

### Usage and examples

Summary statistics, observed and weight, of dataset `sim110` are presented in Table 4 of [Glimm & Yau (2025)](#reference). It is included in the package.

```r
require(maicChecks)
head(sim110)
colnames(sim110)
summary(sim110)
table(sim110$study)
## check if there is an overlap between IPD A and IPD B
ipd1 <- sim110[sim110$study == 'IPD A',]
ipd2 <- sim110[sim110$study == 'IPD B',]
exmLP.2ipd(ipd1 = ipd1, ipd2 = ipd2, 
           vars_to_match = paste0('X', 1:5), 
           cat_vars_to_01 = paste0('X', 1:3), 
           mean.constrained = TRUE)
```

The check returns 0, indicating a solution should exists. Note that by default the additional constraint is set to false (`mean.constrained = FALSE`). In this example, it is added to the check (`mean.constrained = TRUE`)

To perform the match with the additional constraint:

```r
x <- exmWt.2ipd(ipd1 = ipd1, ipd2 = ipd2, 
                vars_to_match = paste0('X', 1:5), 
                cat_vars_to_01 = paste0('X', 1:3), 
                mean.constrained = TRUE)
names(x)
View(x[["ipd1"]])
View(x[["ipd2"]])
View(x[["wtd.summ"]])
```

The result `x` is a list of three objects: `ipd1`, `ipd2`, and `wtd.summ`. The first two are identical dataframes containing the following:

-   The weights for each individual in the two studies    
-   The original input dataframes including all the variables used or not used in the matching    
-   The 0-1 indicators of the categorical variables in `cat_vars_to_01`  

The third object `wtd.summ` contains the effective sample sizes (ESS) for the two studies, and the weighted means of the variables used in matching.

## MAIC


### Methods

-   **Convex Hull Check**: Checks if the AD lies within the convex hull of the IPD; if yes, then it is guaranteed that a unique solution for MAIC weights can be found. This method uses linear programming to determine if the AD is within the convex hull of the IPD, ensuring numerical compatibility for MAIC.
-   **Principal Component Analysis (PCA)**: Provides a visual assessment of the AD's position relative to the IPD in a multi-dimensional space. PCA is used to visualize the AD's position relative to the IPD, providing a graphical representation of data overlap.
-   **Mahalanobis Distance and Hotelling's T² Test**: Tests whether matching IPD to AD is necessary by assessing the similarity of their distributions. These statistical tests assess the similarity between IPD and AD, determining if matching is necessary.

### Usage and examples

``` r
require(maicChecks)
# eAD[1,] is the scenario A in the reference paper,
# i.e. when AD is within IPD convex hull
# eAD[3,] is the scenario C in the reference paper,
# i.e. when AD is outside IPD convex hull

# Perform Convex Hull check
maicLP(eIPD, eAD[1,2:3])
maicLP(eIPD, eAD[3,2:3])

# Visualize data using PCA
a1 <- maicPCA(eIPD, eAD[1,2:3])
a1 ## the dot plots of PC's for IPD and AD
a3 <- maicPCA(eIPD, eAD[3,2:3])
a3 ## the dot plots of PC's for IPD and AD

# Conduct Mahalanobis Distance test
md <- maicMD(eIPD, eAD[1,2:3])
md ## a dot-plot of IPD Mahalanobis distances along with AD in the same metric.

# Conduct Hotelling's T² test
maicT2Test(eIPD, eAD[1,2:3])

# Estimate the MAIC weights
m1 <- maicWt(eIPD, eAD[1,2:3])
```

## Reference

-   Glimm E and Yau L. (2025). "Exact matching as an alternative to propensity score matching." [arXiv:2503.02850v1](https://doi.org/10.48550/arXiv.2503.02850).
-   Glimm E and Yau L. (2022). "Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons." [*Pharmaceutical Statistics*. 21(5):974-987](https://onlinelibrary.wiley.com/doi/full/10.1002/pst.2210).    
-   Signorovitch JE, Wu EQ, Andrew P, et al. (2010). "Comparative effectiveness without head-to-head trials: a method for matching-adjusted indirect comparisons applied to psoriasis treatment with adalimumab or etanercept." *PharmacoEconomics*. 28(10):935-945.

## Package authors

-   Lillian Yau
-   Ekkehard Glimm
-   Xinlei Deng
