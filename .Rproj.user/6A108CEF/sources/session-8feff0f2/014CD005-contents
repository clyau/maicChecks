#' @title Checks whether two IPD datasets can be matched with lpSolve::lp
#'
#' @param ipd1 a dataframe with n1 row and p column, where n1 is number of subjects of the first IPD, and p is the number of variables used in standardization.
#' @param ipd2 a dataframe with n2 row and p column, where n2 is number of subjects of the second IPD, and p is the number of variables used in standardization.
#' @param catigorical.var a list of variable names for the categorical variables in the data
#' @param mean.constrained whether to restrict the weighted means to be within the ranges of observed means. Default is FALSE. When it is TRUE, there is a higher chance of not having a solution.
#'
#' @return \item{lp.check}{0 = OS can be conducted; 2 = OS cannot be conducted}
#' @export exmLP.2ipd
#'
#' @author Lillian Yau
## exmLP(ipd1, ipd2) ## this would be the example, but ipd1 and ipd2 are not in the package yet
exmLP.2ipd <- function (ipd1, ipd2, mean.constrained = FALSE, catigorical.var = NULL)
{
  ipd <- as.data.frame(rbind(-1 * ipd1, ipd2))
  oneszeros <- c(rep(1, nrow(ipd1)), rep(0, nrow(ipd2)))
  zerosones <- c(rep(0, nrow(ipd1)), rep(1, nrow(ipd2)))
  ipd <- as.data.frame(cbind(ipd, oneszeros, zerosones))
  ##
  p <- ncol(ipd1) ## p = number of variables to be matched
  f.con <- as.matrix(t(ipd))
  f.obj <- rep(0.5, ncol(f.con)) ## this is a dummy
  f.rhs <- as.data.frame(t(c(rep(0, p ), 1, 1)))
  f.dir <- rep("=", p+2)
  ##  
  if(mean.constrained){
    ## means of each variable for ipd1 and ipd2
    ipd1.bar <- colMeans(ipd1)
    ipd2.bar <- colMeans(ipd2)
    x <- as.data.frame(rbind(ipd1.bar, ipd2.bar))
    ##
    bar.min <- apply(x, 2, min)
    bar.max <- apply(x, 2, max) 
    ##
    f.rhs <- cbind(f.rhs, 2 * t(bar.min), 2 * t(bar.max))
    f.dir <- c(f.dir, rep(">=", p), rep("<=", p))
    f.con <- data.frame(rbind(f.con, 
                              cbind(t(ipd1), t(ipd2)), 
                              cbind(t(ipd1), t(ipd2))))
    rm(ipd1.bar, ipd2.bar, x, bar.min, bar.max)
  } 
  lp.check <- lpSolve::lp(direction = "max", 
                          objective.in = f.obj, 
                          const.mat = f.con, 
                          const.dir = f.dir, 
                          const.rhs = f.rhs,
                          transpose.constraints = TRUE)$status
  return(list(lp.check = lp.check))
}