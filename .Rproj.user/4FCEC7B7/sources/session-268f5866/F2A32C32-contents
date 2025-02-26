#' @title Maximum ESS Weights
#' @description Estimates an alternative set of weights which maximizes effective sample size (ESS) for a given set of variates used in the matching. Should only be used after it is ascertained that AD is indeed within the convex hull of IPD.
#'
#' @param ipd a dataframe with n row and p column, where n is number of subjects and p is the number of variables used in matching.
#' @param ad a dataframe with 1 row and p column. The matching variables should be in the same order as that in \code{ipd}. The function does not check this.
#'
#' @details The weights maximize the ESS subject to the set of baseline covariates used in the matching.
#' @return
#' \item{maxess.wt }{maximum ESS weights. Scaled to sum up to the total IPD sample size, i.e. nrow(ipd)}
#' \item{ipd.ess }{effective sample size. It is no smaller than the ESS given by the MAIC weights.}
#' \item{ipd.wtsumm}{weighted summary statistics of the matching variables after matching. they should be identical to the input AD when AD is within the IPD convex hull.}
#'
#' @references Glimm & Yau (2021). "Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons", Pharmaceutical Statistics, 21(5):974-987. doi:10.1002/pst.2210.
#' @export maxessWt
#'
#' @examples
#' ## eAD[1,] is scenario A in the reference manuscript
#' m0 <- maxessWt(eIPD, eAD[1,2:3])
maxessWt <- function(ipd, ad) {
  ##
  ## assume ipd is a dataframe with n row and p coln
  ## ... n = number of subjects, p = number of matching variables
  ## assume ad is a dataframe with 1 row and p coln
  ##
  ipd.n <- nrow(ipd)
  ones  <- rep(1, ipd.n)
  ipd0  <- data.frame(cbind(ipd, ones))
  ad0   <- matrix(c(ad, 1, rep(0, nrow(ipd0))), nrow = 1)
  p     <- length(ad) # number of variables used in matching
  #
  ## ipd0 serve as the constraint = Amat
  Amat  <- as.matrix(ipd0) # this includes a row of 1's
  ## optimization object: a quadratic function
  Dmat  <- diag(ipd.n)
  Amat0 <- as.matrix(data.frame(cbind(Amat, Dmat)))
  dvec  <- rep(0, ipd.n) # not in the lp-solve
  ## the right hand side
  bvec <- ad0 # include the 1 already
  #
  x1 <- quadprog::solve.QP(Dmat = Dmat,
                           dvec = dvec,
                           Amat = Amat0,
                           bvec = bvec,
                           meq  = p)
  #
  # weights scaled to total number of patients in ipd
  ipd.wts.me <- x1[["solution"]] * ipd.n
  # ess of new weight
  ipd.ess.me <- round(sum(ipd.wts.me)^2 / sum(ipd.wts.me^2), 1)
  ipd.wtsumm.me <- colMeans(ipd * ipd.wts.me)
  ##
  return(list(maxess.wt = ipd.wts.me,
              ipd.ess = ipd.ess.me,
              ipd.wtsumm = ipd.wtsumm.me))
}
