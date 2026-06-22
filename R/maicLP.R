#' @title Checks if AD is within the convex hull of IPD using lp-solve
#' @description Checks if AD is within the convex hull of IPD using lp-solve
#'
#' @param ipd a dataframe with n row and p column, where n is number of subjects and p is the number of variables used in matching.
#' @param ad a dataframe with 1 row and p column. The matching variables should be in the same order as that in \code{ipd}. The function does not check this.
#'
#' @return
#' \item{lp.check }{0 = AD is inside IPD, and MAIC can be conducted; 2 = otherwise}
#'
#' @references Glimm E and Yau L. (2022). 'Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons.' \emph{Pharmaceutical Statistics}, 21(5):974-987. \doi{10.1002/pst.2210}.
#'
#' @export maicLP
#'
#' @examples
#' ## eIPD now contains response columns (r.cont, r.bin) in addition to the
#' ## matching columns y1, y2. Subset to the matching columns explicitly.
#'
#' ## eAD[1,] is the scenario A in the reference paper,
#' ## i.e. when AD is within IPD convex hull
#' maicLP(eIPD[, c('y1', 'y2')], eAD[1, c('y1', 'y2')])
#'
#' ## eAD[3,] is the scenario C in the reference paper,
#' ## i.e. when AD is outside IPD convex hull
#' maicLP(eIPD[, c('y1', 'y2')], eAD[3, c('y1', 'y2')])
maicLP <- function(ipd, ad) {
  ##
  ## assume ipd is a dataframe with n row and p coln
  ## ... n = number of subjects, p = number of matching variables
  ## assume ad is a dataframe / numeric vector with p elements
  ##
  ipd <- as.matrix(ipd)
  ad  <- as.numeric(unlist(ad))
  n   <- nrow(ipd)
  p   <- ncol(ipd)
  if (length(ad) != p)
    stop('`ad` must have the same number of variables as `ipd` (got ',
         length(ad), ' vs ', p, ').', call. = FALSE)

  ## stack a row of 1's to enforce sum(w) = 1
  f.con    <- rbind(t(ipd), rep(1, n))
  f.obj    <- rep(0.5, n)               ## dummy objective
  f.rhs    <- c(ad, 1)
  f.dir    <- rep('=', p + 1)

  lp.check <- lpSolve::lp(direction    = 'max',
                          objective.in = f.obj,
                          const.mat    = f.con,
                          const.dir    = f.dir,
                          const.rhs    = f.rhs)$status
  return(list(lp.check = lp.check))
}
