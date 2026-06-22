#' @title Estimates the MAIC weights
#' @description Estimates the MAIC weights for each individual in the IPD. The function first calls \code{\link{maicLP}} as a feasibility gate; if \code{ad} is not within the convex hull of \code{ipd} the function stops with an informative error rather than letting \code{optim()} fail silently.
#'
#' @param ipd a dataframe with n row and p column, where n is number of subjects and p is the number of variables used in matching.
#' @param ad a dataframe with 1 row and p coln. The matching variables should be in the same order as that in \code{ipd}. The function does not check this.
#' @param max.it maximum iteration passed to optim(). if \code{ad} is within \code{ipd} convex hull, then the default 25 iterations of optim() should be enough.
#'
#' @return The main code are taken from Philippo (2016). It returns the following:
#' \item{optim.out}{results of optim()}
#' \item{maic.wt}{MAIC un-scaled weights for each subject in the IPD set}
#' \item{maic.wt.rs}{re-scaled weights which add up to the original total sample size, i.e. nrow(ipd)}
#' \item{ipd.ess}{effective sample size}
#' \item{ipd.wtsumm}{weighted summary statistics of the matching variables after matching. they should be identical to the input AD when AD is within the IPD convex hull.}
#'
#' @references Phillippo DM, Ades AE, Dias S, et al. (2016). Methods for population-adjusted indirect comparisons in submissions to NICE. NICE Decision Support Unit Technical Support Document 18.
#' @export maicWt
#'
#' @examples
#' ## eIPD contains response columns (r.cont, r.bin) in addition to the
#' ## matching columns y1, y2. Subset to the matching columns explicitly.
#' ## eAD[1,] is scenario A in the reference manuscript
#' m1 <- maicWt(eIPD[, c('y1', 'y2')], eAD[1, c('y1', 'y2')])
maicWt <- function(ipd, ad, max.it = 25) {
  ##
  ## assume ipd is a dataframe with n row and p coln
  ## ... n = number of subjects, p = number of matching variables
  ## assume ad is a dataframe with 1 row and p coln
  ##
  ## ---- feasibility gate -----------------------------------------
  lp.chk <- maicLP(ipd, ad)$lp.check
  if (lp.chk != 0)
    stop('AD is outside the convex hull of IPD (maicLP status = ', lp.chk,
         '). MAIC weights cannot be computed; matching is not feasible.',
         call. = FALSE)

  objfn  <- function(a1, X) sum(exp(X %*% a1))
  gradfn <- function(a1, X) colSums(X * as.vector(exp(X %*% a1)))

  ipd     <- as.data.frame(ipd)
  p       <- ncol(ipd)
  ipd.n   <- nrow(ipd)
  X.EM.1  <- sweep(as.matrix(ipd), 2, as.numeric(unlist(ad)), '-')

  ## Estimate alpha_2 (See Philippo 2016)
  op2 <- stats::optim(par     = rep(0, p),
                      fn      = objfn,
                      gr      = gradfn,
                      X       = X.EM.1,
                      method  = 'BFGS',
                      control = list(maxit = max.it))

  a2    <- op2$par
  wt    <- exp(X.EM.1 %*% a2)               ## un-scaled weights
  wt.rs <- (wt / sum(wt)) * ipd.n           ## rescaled to sum to n

  ipd.ess    <- round(sum(wt.rs)^2 / sum(wt.rs^2), 1)
  ipd.wtsumm <- colMeans(ipd * as.vector(wt.rs))

  return(list(optim.out  = op2,
              maic.wt    = wt,
              maic.wt.rs = wt.rs,
              ipd.ess    = ipd.ess,
              ipd.wtsumm = ipd.wtsumm))
}
