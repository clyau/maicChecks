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
| `maicChecks`&nbsp;&nbsp;&nbsp;| Release (version 0.2.0)&nbsp;&nbsp;&nbsp;| CRAN&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;| `install.packages("maicChecks")` |

## Version history

### Version 0.2.0

**_CURRENT:_ Version 0.2.0** was released on CRAN on 3 March, 2025. The following method is added:

-   [Exact matching](#exact-matching) by maximizing effective sample size (ESS) for comparing two studies both with patient level data ([Glimm and Yau (2025)](#reference))

### Version 0.1.2

Initial **version 0.1.2** was released on CRAN on 4 January, 2022. It implements the following methods:

-   Matching-adjusted indirect comparison (MAIC) for comparing study with patient level data to study with only aggregated data ([Signorovitch (2010)](#reference))    
-   MAIC by maximizing effective sample size (ESS) for comparing study with patient level data to study with only aggregated data ([Glimm and Yau (2022)](#reference))     
-   Various checks to assess feasibility of conducting the above two methods ([Glimm and Yau (2022)](#reference))

All functionalities in the initial version are retained in the current version 0.2.0.

## Overview

The comparison of different medical treatments from observational studies or across different clinical studies is often biased by confounding factors such as systematic differences in patient demographics or in the inclusion criteria for the trials. The confounding must be adjusted before indirect comparisons can be conducted. The adjustment is usually accomplished by matching the baseline covariates so that patients from one or both studies are each assigned a weight, which are taken into account when comparing clinical outcomes.

From a data availability prespective, two situations arise:

1.  Individual patient-level data (IPD) are available for both studies, i.e., an IPD vs IPD comparison.
2.  IPD are available for one study, and only aggregated data (AD), i.e., summary statistics, are available for the other, i.e., and IPD vs AD comparison.

In both situations, one study can be considered a "target" population, and the other is then matched onto it. However, in the first situation when IPD are available in both, the two study populations can also be matched onto a common pooled population. 

`maicChecks` is an R package that offers two different but related methods for matching baseline covariates for these two situations. The methods are:

1.  [Exact matching](#exact-matching): used for IPD vs IPD   
2.  [MAIC](#matching-adjusted-indirect-comparison-maic): used for IPD vs AD

## Exact matching

Exact matching can be used when IPD are available from both studies. It is an alternative to propensity score matching. The method ensures that after matching, the weighted means of the baseline covariates between the two studies are exactly the same. Details on methodology can be found in [Glimm & Yau (2025)](#reference).

### Method: Linear check

The function `maicChecks::exmLP.2ipd()` checks if there is an overlap between the two IPD. If yes, matching of baseline covariates can be performed. In practice, there is almost always an overlap between the two IPD. 

#### Usage and example

Included in the package is data set `sim110`. The summary statistics, observed and weight, of the dataset are presented in Table 4 of [Glimm & Yau (2025)](#reference). 

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

### Method: Exact matching for IPD vs IPD

The function `maicChecks::exmWt.2ipd()` matches the baseline covariates from the two IPD, and assigns a weight to each patient in the two studies. The algorithm treats matching as a constrained optimization problem Constrained optimization is a purely algebraic technique and solves the convex optimization problem in a finite number of steps. In other words, the method does not require numerical approximation.

Although the weighted means are the same for the two studies after matching, it can happen that in one of the covariates used in matching it is not between the two observed means. To avoid this, an additional constraint can be added to force the weighted means to be always between the observed means. Naturally, with an additional constraint, the likelihood of a non-existing solution increases.

(Note that `maicChecks::exmLP.2ipd()` which checks for feasibility of a match uses linear programming. Hence, it can happen that the check thinks it's fine to go ahead with the matching, but `maicChecks::exmWt.2ipd()` yields no solution.)

#### Usage and example

The following code perform the matching with the additional constraint. To Perform the matching without it, leave `mean.constrained = FALSE` as is the default.

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

The weights for each individual in the output can be used directly for desired analysis comparing the outcomes of the two studies.

#### Note on matching one IPD onto the other

In the above example, the two IPD's are matched onto a common (pooled) population. However, if one IPD represents the target population and should be matched onto, then the summary statistics of the baseline covariates can be first calculated for this IPD. Subsequently, the MAIC approach described in the next section can be applied. 

## Matching-Adjusted Indirect Comparison (MAIC)

This method is used when IPD is available for one study but only AD is available for the other.

The package focuses on the following two areas of methods related to MAIC:

1.   Checking feasibility of conducting MAIC
2.   Implementing MAIC:    
   (a)   as originally proposed in [Signorovitch (2010)](#reference)    
   (b)   by maximizing ESS as proposed in [Glimm and Yau (2022)](#reference)    

### Method: checking feasibility of conducting MAIC

Movitation and methods for checking whether MAIC are described in [Glimm and Yau (2022)](#reference). The three methods implemented are summarized below:

-   **Convex Hull Check**: Checks if the AD lies within the convex hull of the IPD; if yes, then it is guaranteed that a unique solution for MAIC weights can be found. This method uses linear programming to determine if the AD is within the convex hull of the IPD, ensuring numerical compatibility for MAIC.
-   **Principal Component Analysis (PCA)**: Provides a visual assessment of the AD's position relative to the IPD in a multi-dimensional space. PCA is used to visualize the AD's position relative to the IPD, providing a graphical representation of data overlap.
-   **Mahalanobis Distance and Hotelling's T² Test**: Tests whether matching IPD to AD is necessary by assessing the similarity of their distributions. These statistical tests assess the similarity between IPD and AD, determining if matching is necessary.

#### Usage and examples

The simulated datasets presented in Section 2 of [Glimm and Yau (2022)](#reference) are used here to illustrate syntax. These datasets are also included in the package.

A couple simple examples:

``` r
require(maicChecks)
# eAD[1,] is the scenario A in the reference paper,
# i.e. when AD is within IPD convex hull
# eAD[3,] is the scenario C in the reference paper,
# i.e. when AD is outside IPD convex hull

print(eAD)
head(eIPD) ## the IPD dataset
dim(eIPD)

# Perform Convex Hull check
maicLP(eIPD, eAD[1,2:3])
maicLP(eIPD, eAD[3,2:3])

# Visualize data using PCA
a1 <- maicPCA(eIPD, eAD[1,2:3])
print(a1) ## the dot plots of PC's for IPD and AD
a3 <- maicPCA(eIPD, eAD[3,2:3])
print(a3) ## the dot plots of PC's for IPD and AD

# Conduct Mahalanobis Distance test
md <- maicMD(eIPD, eAD[1,2:3])
md ## a dot-plot of IPD Mahalanobis distances along with AD in the same metric.

# Conduct Hotelling's T² test
maicT2Test(eIPD, eAD[1,2:3])
```

Two points to note:  
1.   It is important that all variables are in the same orders in IPD and in AD. The functions assume this is the case but do not check it.
2.   If there are categorical variables (e.g., region, disease status, or median or quantiles of a continuous variable) to be used in matching, they need to be first convert to 0-1 indicator variables.

**_Converting categorical variables before checking or matching:_**

As an example, two additional variables are added in `eIPD` and `eAD[1,]` to be also used in matching: (1) patients' baseline disease risk category (low, intermediate, high), and (2) another continuous variable (y3) but instead of the mean, only median is available in AD.

The new AD data based on `eAD[1,]`:

``` r
eAD.1x <- data.frame(c(eAD[1,2:3], 
           ds.low = 0.24, ds.int = 0.58, ds.hi = 0.18, 
           y3 = 0.1))
print(eAD.1x)
```

Here, `y1` and `y2` are, as before, the means of two continous covariates; `ds.low`, `ds.int`, and `ds.hi` are proportions of patients in disease status categories low, intermediate, or high, respectively; and `y3` is the **median** of a continuous variable `y3`.

Before `eAD.1x` can be used for checking/matching, `y3` need to be converted to the proportion of patients above (or below) the median. In this case, the proportion is obviously 50%.

``` r
eAD.1x <- eAD.1x %>% 
  mutate(y3.med = 0.50) %>%
  select(-y3, -ds.hi) ## remove y3 and ds.hi
```

Two new variables `ds` and `y3` are also added to the `eIPD` dataset:

``` r
eIPD.x <- data.frame(ds = sample(c('low', 'intermediate', 'high'),
                                 size = nrow(eIPD), 
                                 replace = TRUE),
                     y3 = rnorm(n = nrow(eIPD), mean = 0, sd = 1.2)
                     ) %>%
  cbind(eIPD, .)
##
head(eIPD.x)
```

Before checking/matching can be performed, indicator variables must be created for `ds` and `y3` median:

``` r
## indicators are created for ds low and ds intermediate 
eIPD.x <- eIPD.x %>% 
  mutate(ds.low = ifelse(ds == 'low', 1, 0),
         ds.int = ifelse(ds == 'intermediate', 1, 0)) %>%
  ## y3.med is the indicator for whether a patient in eIPD whose y3 values ...
  ## ... are below 0.1, the median in AD study.
  mutate(y3.med = ifelse(y3 <= 0.1, 1, 0)) %>%  
  select(-y3, -ds) ## remove y3 and ds
## make sure the variables are in the same order
head(eIPD.x)
print(eAD.1x)
```

The same syntax is used for checking for the feasibility.

``` r
# Perform Convex Hull check
maicLP(eIPD.x, eAD.1x)

# Visualize data using PCA
a1.x <- maicPCA(eIPD.x, eAD.1x)
print(a1.x)

# Conduct Mahalanobis Distance test
md.x <- maicMD(eIPD.x, eAD.1x)
md.x

# Conduct Hotelling's T² test
maicT2Test(eIPD.x, eAD.1x)
```

Other than `maicMD()`, other functions (including matching, see below) do not require the indicator variables to correspond to a full-rank design matrix, i.e. with `k-1` indicators for a categorical variable with `k` levels. In other words, if `maicMD()` is not needed, it is fine to have `k` indicators for a `k`-level categorical variable. 

### Method: MAIC as proposed by [Signorovitch (2010)](#reference). 

Details of the method can be found in the reference publication. The following snytax can be used to perform the matching.

``` r
# Estimate the MAIC weights
m1 <- maicWt(eIPD, eAD[1,2:3]) ## default max.it = 25
```

The output `m1` is a list contains results inherited from `optim()` function. Although in most cases the default `max.it = 25` is sufficient when the checks have passed, **_it is still important to first check whether convergence is achieved:_** if `m1$optim.out$convergence = 1` then `optim()` has converged; if `m1$optim.out$convergence = 0` then no, and `max.it` needs to be increased.

The rest are related to the matching:

-    `maic.wt`: a vector with weights each patient in IPD study receives after matching
-    `maic.wt.rs`: re-scaled weights so that the sum is the total number of patients in IPD. It is recommended to use this weight.
-    `ipd.ess`: ESS for the IPD study
-    `ipd.wtsumm`: weighted means of the matching variables. These should be identical to the values of AD.

For the datasets with categorical variables, the matching are done the same way:

``` r
## make sure the variables are in the same order
head(eIPD.x)
print(eAD.1x)

# Estimate the MAIC weights
m1.x <- maicWt(eIPD.x, eAD.1x)
```

### Method: MAIC by maximizing ESS as introduced in [Glimm and Yau (2022)](#reference)

Details of the method can be found in the reference publication. The following snytax can be used to perform the matching.

``` r
# Estimate the MAIC weights
me1 <- maxessWt(eIPD, eAD[1,2:3])

# The example with categorical variables
me1.x <- maxessWt(eIPD.x, eAD.1x)
```

The outputs `me1` and `me1.x` contain the following:

-    `maxess.wt`: (re-scaled) weights each patient in IPD study receives after matching. The total add up to the number of patients in IPD.
-    `ipd.ess`: ESS for the IPD study
-    `ipd.wtsumm`: weighted means of the matching variables. These should be identical to the values of AD.

## Reference

-   Glimm E and Yau L. (2025). "Exact matching as an alternative to propensity score matching." [*Statistics in Biopharmaceutical Research*](https://doi.org/10.1080/19466315.2025.2507378).
-   Glimm E and Yau L. (2022). "Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons." [*Pharmaceutical Statistics*. 21(5):974-987](https://onlinelibrary.wiley.com/doi/full/10.1002/pst.2210).    
-   Signorovitch JE, Wu EQ, Andrew P, et al. (2010). "Comparative effectiveness without head-to-head trials: a method for matching-adjusted indirect comparisons applied to psoriasis treatment with adalimumab or etanercept." *PharmacoEconomics*. 28(10):935-945.

## Package authors

-   Lillian Yau
-   Ekkehard Glimm
-   Xinlei Deng
