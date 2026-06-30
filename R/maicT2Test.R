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
#' @references Glimm E and Yau L. (2022). 'Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons.' \emph{Pharmaceutical Statistics}, 21(5):974-987. \doi{10.1002/pst.2210}.
#' @export maicT2Test
#'
#' @examples
#' ## eIPD contains response columns; subset to matching columns y1, y2.
#' ## eAD[1,] is the scenario A in the reference paper,
#' ## i.e. when AD is perfectly within IPD
#' maicT2Test(eIPD[, c('y1', 'y2')], eAD[1, c('y1', 'y2')])
maicT2Test <- function(ipd, ad, n.ad = Inf) {
  ipd     <- as.matrix(ipd)
  ad      <- as.numeric(unlist(ad))
  n.ipd   <- nrow(ipd)
  p.var   <- ncol(ipd)
  if (length(ad) != p.var)
    stop('`ad` must have the same number of variables as `ipd`.', call. = FALSE)

  ipd.bar <- colMeans(ipd)
  S       <- stats::cov(ipd)
  diff    <- ipd.bar - ad

  ## scaling factor: n.ipd when AD is fixed; (n.ipd*n.ad)/(n.ipd+n.ad) otherwise
  scale <- if (is.finite(n.ad)) (n.ipd * n.ad) / (n.ipd + n.ad) else n.ipd
  T.squared   <- scale * as.numeric(crossprod(diff, solve(S, diff)))
  T.squared.f <- T.squared * (n.ipd - p.var) / (p.var * (n.ipd - 1))
  p.val       <- 1 - stats::pf(T.squared.f, p.var, n.ipd - p.var)

  cat('T.sq.f = ', T.squared.f, '\n')
  cat('p.val  = ', p.val, '\n')

  invisible(list(T.sq.f = T.squared.f,
                 p.val  = p.val))
}
