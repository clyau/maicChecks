#' @title Checks if AD is within the convex hull of IPD using lp-solve
#' @description Checks if AD is within the convex hull of IPD using lp-solve
#'
#' @param ipd a dataframe with n row and p column, where n is number of subjects and p is the number of variables used in matching.
#' @param ad a dataframe with 1 row and p column. The matching variables should be in the same order as that in \code{ipd}. The function does not check this.
#'
#' @return
#' \item{lp.check }{0 = AD is inside IPD, and MAIC can be conducted; 2 = otherwise}
#'
#' @references Glimm & Yau (2021). "Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons", Pharmaceutical Statistics, 21(5):974-987. doi:10.1002/pst.2210.
#'
#' @export maicLP
#'
#' @examples
#' ## eAD[1,] is the scenario A in the reference paper,
#' ## i.e. when AD is within IPD convex hull
#' maicLP(eIPD, eAD[1,2:3])
#'
#' ## eAD[3,] is the scenario C in the reference paper,
#' ## i.e. when AD is outside IPD convex hull
#' maicLP(eIPD, eAD[3,2:3])
maicLP <- function(ipd, ad) {
  ##
  ## assume ipd is a dataframe with n row and p coln
  ## ... n = number of subjects, p = number of matching variables
  ## assume ad is a dataframe with 1 row and p coln
  ##
  ## constrain the sum of the weights to 1
  ones  <- rep(1, nrow(ipd))
  ipd   <- data.frame(cbind(ipd, ones))
  p     <- ncol(ad) ## p = number of variables to match
  ad    <- data.frame(c(ad, 1))
  ## the ipd serve as the constraint
  f.con <- as.matrix(t(ipd))
  ## a dummy object to be optimized
  f.obj <- rep(0.5, ncol(f.con))
  ## the right hand side is ad
  f.rhs <- ad
  ## direction of constraint
  f.dir <- rep("=", p+1)
  ## solve
  lp.check <- lpSolve::lp (direction = "max", 
                           objective.in = f.obj, 
                           const.mat = f.con, 
                           const.dir = f.dir, 
                           const.rhs = f.rhs)$status
  ##
  return(list(lp.check = lp.check))
}
