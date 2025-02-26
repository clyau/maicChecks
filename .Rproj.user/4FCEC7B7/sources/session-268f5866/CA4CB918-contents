#' @title Hotelling's T-square test to check whether maic is needed
#'
#' @param ipd a dataframe with n row and p column, where n is number of subjects and p is the number of variables used in matching.
#' @param ad a dataframe with 1 row and p column. The matching variables should be in the same order as that in \code{ipd}. The function does not check this.
#' @param n.ad default is Inf assuming \code{ad} is a fixed (known) quantity with infinit accuracy. In most MAIC applications \code{ad} is the sample statistics and \code{n.ad} is known.
#'
#' @details When \code{n.ad} is not Inf, the covariance matrix is adjusted by the factor n.ad/(n.ipd + n.ad)), where n.ipd is nrow(ipd), the sample size of \code{ipd}.
#' @return
#' \item{T.sq.f }{the value of the T^2 test statistic}
#' \item{p.val }{the p-value corresponding to the test statistic. When the p-value is small, matching is necessary.}
#' @references Glimm & Yau (2021). "Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons", Pharmaceutical Statistics, 21(5):974-987. doi:10.1002/pst.2210.
#' @export maicT2Test
#'
#' @examples
#' ## eAD[1,] is the scenario A in the reference paper,
#' ## i.e. when AD is perfectly within IPD
#' maicT2Test(eIPD, eAD[1,2:3])
maicT2Test <- function(ipd, ad, n.ad = Inf){
  ##
  ipd.bar <- colMeans(ipd)
  n.ipd <- nrow(ipd)
  p.var <- ncol(ipd)
  ## cov.matrix of ipdx
  sig.hat <- (1 / (n.ipd - 1)) * t(ipd) %*%
    (diag(n.ipd) - matrix(1/n.ipd,
                          nrow = n.ipd,
                          ncol = n.ipd)) %*%
    as.matrix(ipd)
  ## hotelling's t-squared
  if(n.ad == Inf){ # when ad is fixed
    T.squared <- n.ipd * as.matrix(ipd.bar - ad) %*%
      solve(sig.hat) %*%
      as.matrix(t(ipd.bar - ad))
  } else {
    T.squared <- ((n.ipd * n.ad) / (n.ipd + n.ad)) *
      as.matrix(ipd.bar - ad) %*%
      solve(sig.hat) %*%
      as.matrix(t(ipd.bar - ad))
  }
  ## hotelling's t-squared adjusting the factor to have an F-distn
  T.squared.f <- T.squared * (n.ipd - p.var) / (p.var * (n.ipd - 1))
  ##
  cat('T.sq.f = ', T.squared.f, '\n')
  p.val <- 1 - pf(T.squared.f, p.var, (n.ipd - p.var))
  cat('p.val = ', p.val, '\n')
}
