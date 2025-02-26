# maicChecks <a href="https://github.com/clyau/maicChecks"><img src="https://github.com/user-attachments/assets/16780205-98c6-4bbd-af23-dfec94010d54" align="right" height="140"/></a>

<!-- badges: start -->

[![CRAN status](https://badges.cranchecks.info/flavor/release/maicChecks.svg)](https://cran.r-project.org/web/checks/check_results_maicChecks.html) [![downloads](https://cranlogs.r-pkg.org/badges/maicChecks)](https://www.rdocumentation.org/trends)

<!-- badges: end -->

`maicChecks`: This package provides tools for data analysts to conduct additional statistical checks and ensure the numerical feasibility of implementing MAIC, a method widely used in clinical trials and observational studies for comparing treatment effects.

## Installation

| Type                | Source | Command                                       |
|------------------|------------------|------------------------------------|
| Release (0.1.2)     | CRAN   | `install.packages("maicChecks")`              |
| Development (0.2.0) | GitHub | `remotes::install_github("clyau/maicChecks")` |

## Overview

The `maicChecks` package offers a comprehensive set of tools to assess the appropriateness of MAIC for a given set of baseline characteristics. The package includes methods for checking the overlap between individual patient data (IPD) and aggregated data (AD), visualizing data distributions, and performing statistical tests to ensure valid comparisons.

## Methods

-   **Convex Hull Check**: Ensures that the AD lies within the convex hull of the IPD, guaranteeing a unique solution for MAIC weights. This method uses linear programming to determine if the AD is within the convex hull of the IPD, ensuring numerical compatibility for MAIC.
-   **Principal Component Analysis (PCA)**: Provides a visual assessment of the AD's position relative to the IPD in a multi-dimensional space. PCA is used to visualize the AD's position relative to the IPD, providing a graphical representation of data overlap.
-   **Mahalanobis Distance and Hotelling's T² Test**: Tests whether matching IPD to AD is necessary by assessing the similarity of their distributions. These statistical tests assess the similarity between IPD and AD, determining if matching is necessary.
-   **Exact Matching**: An alternative to propensity score matching, offering exact matching of covariate means between treatment groups. This method provides exact matching of covariate means between treatment groups, ensuring balanced comparisons.

## Usage/Example
``` r

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
-   [Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons](https://onlinelibrary.wiley.com/doi/full/10.1002/pst.2210)

## Package authors

-   Lillian Yau
-   Ekkehard Glimm
-   Xinlei Deng

