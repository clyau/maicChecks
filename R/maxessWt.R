#' @title Maximum ESS Weights
#' @description Estimates an alternative set of weights which maximizes effective sample size (ESS) for a given set of variates used in the matching. The function first calls \code{\link{maicLP}} as a feasibility gate; if \code{ad} is not within the convex hull of \code{ipd} the function stops with an error message.
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
#' @references Glimm E and Yau L. (2022). 'Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons.' \emph{Pharmaceutical Statistics}, 21(5):974-987. \doi{10.1002/pst.2210}.
#' @export maxessWt
#'
#' @examples
#' ## eIPD contains response columns (r.cont, r.bin) in addition to the
#' ## matching columns y1, y2. Subset to the matching columns explicitly.
#' ## eAD[1,] is scenario A in the reference manuscript
#' m0 <- maxessWt(eIPD[, c('y1', 'y2')], eAD[1, c('y1', 'y2')])
maxessWt <- function(ipd, ad) {
  ##
  ## assume ipd is a dataframe with n row and p coln
  ## ... n = number of subjects, p = number of matching variables
  ## assume ad is a dataframe with 1 row and p coln
  ##
  ## ---- feasibility gate -----------------------------------------
  lp.chk <- maicLP(ipd, ad)$lp.check
  if (lp.chk != 0)
    stop('AD is outside the convex hull of IPD (maicLP status = ', lp.chk,
         '). max-ESS weights cannot be computed; matching is not feasible.',
         call. = FALSE)

  ipd   <- as.matrix(ipd)
  ad    <- as.numeric(unlist(ad))
  ipd.n <- nrow(ipd)
  p     <- length(ad)

  ## constraints: t(A) %*% w >= b, first p+1 are equality (meq = p+1 in concept;
  ## the original code uses meq = p with a redundant sum>=1 constraint -- preserved)
  Amat <- cbind(ipd, rep(1, ipd.n), diag(ipd.n))
  bvec <- c(ad, 1, rep(0, ipd.n))
  Dmat <- diag(ipd.n)
  dvec <- rep(0, ipd.n)

  x1 <- quadprog::solve.QP(Dmat = Dmat,
                           dvec = dvec,
                           Amat = Amat,
                           bvec = bvec,
                           meq  = p)

  w.sol <- x1[['solution']]
  tol   <- max(1e-12, 1e-10 * max(1, max(abs(w.sol))))

  if (any(w.sol < -tol))
    stop('`maxessWt` produced materially negative weights; check feasibility / numerics.',
         call. = FALSE)

  if (any(w.sol < 0)) {
    warning('`maxessWt` produced tiny negative weights from numerical noise; set to 0.',
            call. = FALSE)
    w.sol[w.sol < 0] <- 0
    if (sum(w.sol) <= 0)
      stop('weights must sum to a positive number after clipping.', call. = FALSE)
    w.sol <- w.sol / sum(w.sol)
  }

  ipd.wts.me    <- w.sol * ipd.n            ## scaled to sum to n
  ipd.ess.me    <- round(sum(ipd.wts.me)^2 / sum(ipd.wts.me^2), 1)
  ipd.wtsumm.me <- colMeans(ipd * ipd.wts.me)

  return(list(maxess.wt  = ipd.wts.me,
              ipd.ess    = ipd.ess.me,
              ipd.wtsumm = ipd.wtsumm.me))
}
