% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maxessWt.R
\name{maxessWt}
\alias{maxessWt}
\title{Maximum ESS Weights}
\usage{
maxessWt(ipd, ad)
}
\arguments{
\item{ipd}{a dataframe with n row and p column, where n is number of subjects and p is the number of variables used in matching.}

\item{ad}{a dataframe with 1 row and p column. The matching variables should be in the same order as that in \code{ipd}. The function does not check this.}
}
\value{
\item{maxess.wt }{maximum ESS weights. Scaled to sum up to the total IPD sample size, i.e. nrow(ipd)}
\item{ipd.ess }{effective sample size. It is no smaller than the ESS given by the MAIC weights.}
\item{ipd.wtsumm}{weighted summary statistics of the matching variables after matching. they should be identical to the input AD when AD is within the IPD convex hull.}
}
\description{
Estimates an alternative set of weights which maximizes effective sample size (ESS) for a given set of variates used in the matching. Should only be used after it is ascertained that AD is indeed within the convex hull of IPD.
}
\details{
The weights maximize the ESS subject to the set of baseline covariates used in the matching.
}
\examples{
## eAD[1,] is scenario A in the reference manuscript
m0 <- maxessWt(eIPD, eAD[1,2:3])
}
\references{
Glimm & Yau (2021). "Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons", Pharmaceutical Statistics, 21(5):974-987. doi:10.1002/pst.2210.
}
