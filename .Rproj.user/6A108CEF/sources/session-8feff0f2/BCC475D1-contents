#' @title Exact matching for two IPD's
#'
#' @param ipd1 a dataframe with n row and p column, where n is number of subjects and p is the number of variables used in matching.
#' @param ipd2 the other IPD with the same number of columns
#' @param catigorical.var a list of variable names for the categorical variables in the data
#' @param mean.constrained whether to restrict the weighted means to be within the ranges of observed means. Default is FALSE. When it is TRUE, there is a higher chance of not having a solution.
#' 
#' @details If dummy variables are already created for the categorical variables in the data set, and are present in \code{ipd1} and \code{ipd2}, then \code{categorical.var} should be left as NULL.
#' 
#' @return
# \item{osess.wt }{weights of optimal standardization by maximizing ESS. Scaled to .......}
# \item{ipd.ess }{effective sample size. It is no smaller than the ESS given by the MAIC weights.}
# \item{ipd.wtsumm}{weighted summary statistics of the matching variables after matching. they should be identical to the input AD when AD is within the IPD convex hull.}
#' \item{wt.1 }{weights of the exact matching by maximizing ESS for IPD 1}
#' \item{ess.1 }{effective sample size for IPD 1}
#' \item{wtd.means.1 }{weighted means of the matching variables for IPD 1. Same as wtd.means.2}
#' \item{wt.2 }{weights of the exact matching by maximizing ESS for IPD 2}
#' \item{ess.2 }{effective sample size for IPD 2}
#' \item{wtd.means.2 }{weghted means of the matching bariables for IPD 2. Same as wtd.means.1}
#' 
#' @export exmWt.2ipd
#'
exmWt.2ipd <- function (ipd1, ipd2, mean.constrained = FALSE, catigorical.var = NULL) 
{
  ipd <- as.data.frame(rbind(-1 * ipd1, ipd2))
  oneszeros <- c(rep(1, nrow(ipd1)), rep(0, nrow(ipd2)))
  zerosones <- c(rep(0, nrow(ipd1)), rep(1, nrow(ipd2)))
  ipd <- as.data.frame(cbind(ipd, oneszeros, zerosones))
  rm(oneszeros, zerosones)
  ipd.n <- nrow(ipd)
  p <- ncol(ipd1) ## number of matching variables
  ##
  ## per default no constraint on the mean
  ## then bvec and Amat0 do not have the two last colns (see below when constrain is true)
  ##
  bvec <- ad0 <- matrix(c(rep(0, p), 
                          1, 
                          1, 
                          rep(0, nrow(ipd))), 
                        nrow = 1)
  Amat <- as.matrix(ipd)
  Dmat <- diag(ipd.n)
  Amat0 <- as.matrix(data.frame(cbind(Amat, Dmat)))
  ## dvec is not really needed for this purpose
  dvec <- rep(0, ipd.n)
  ##
  if(mean.constrained){
    ## means of each variable for ipd1 and ipd2
    ipd1.bar <- colMeans(ipd1)
    ipd2.bar <- colMeans(ipd2)
    x <- as.data.frame(rbind(ipd1.bar, ipd2.bar))
    ##
    bar.min <- apply(x, 2, min)
    bar.max <- apply(x, 2, max) 
    bvec <- matrix(c(ad0, 2 * bar.min, 2 * bar.max * (-1)), nrow = 1)
    rm(ipd1.bar, ipd2.bar, x, bar.min, bar.max)
    ## ipd blocks
    x0 <- as.data.frame(rbind(ipd1, ipd2))
    x1 <- as.data.frame(rbind(-1 * ipd1, -1 * ipd2))
    Amat0 <- as.matrix(data.frame(cbind(Amat, Dmat, x0, x1)))
    rm(x0, x1)
  } 
  ## Dmat = Q in the (soon to be) paper
  ## dvec = b in the (soon to be) paper, which are all 0's
  ## Amat = Amat0 = A in (soon to be) paper
  ## bvec = c in the (soon to be) paper, plus the two additional restrictions (if mean.constrained == TRUE) to ensure the weighted means are within observed means
  ## meq = signs of equality/inequality, which we have p+2 equals, 
  ## ... and there is one '>=' which in the syntax do not need to be specified
  wts <- quadprog::solve.QP(Dmat = Dmat, 
                            dvec = dvec, 
                            Amat = Amat0, 
                            bvec = bvec, 
                            meq = p+2)
  ## all weights
  ipd.wts.me <- wts[["solution"]]
  ipd.ess.me <- round(sum(ipd.wts.me)^2/sum(ipd.wts.me^2), 1)
  ipd.wt.mean.me <- colMeans(ipd[, 1:p] * ipd.wts.me)
  ## weights etc. for ipd1
  ipd.1.wts <- ipd.wts.me[1 : nrow(ipd1)]
  ipd.1.wts <- ipd.1.wts * nrow(ipd1)
  ipd.1.ess <- round(sum(ipd.1.wts)^2/sum(ipd.1.wts^2), 1)
  ipd.1.wt.mean <- colMeans(ipd1 * ipd.1.wts)
  ## weights etc. for ipd2
  ipd.2.wts <- ipd.wts.me[(1 + nrow(ipd1)) : (nrow(ipd1) + nrow(ipd2))]
  ipd.2.wts <- ipd.2.wts * nrow(ipd2) 
  ipd.2.ess <- round(sum(ipd.2.wts)^2/sum(ipd.2.wts^2), 1)
  ipd.2.wt.mean <- colMeans(ipd2 * ipd.2.wts)
  ##
##<<<<<<< HEAD:R/osWt.2ipd.R
  return(list(ipd.all.wts = ipd.wts.me, 
              ipd.all.ess = ipd.ess.me, 
              ipd.all.wts.mean = ipd.wt.mean.me,
              ipd.1.wts = ipd.1.wts,
              ipd.1.ess = ipd.1.ess,
              ipd.1.wts.mean = ipd.1.wt.mean,
              ipd.2.wts = ipd.2.wts,
              ipd.2.ess = ipd.2.ess,
              ipd.2.wts.mean = ipd.2.wt.mean
# =======
#   return(list(## wt.all = ipd.wts.me, ## ipd.all.wt, weights from the 2 ipd combined
#               ## ess.all = ipd.ess.me, ## ipd.all.ess, ess from the 2 ipd combined
#               ## wtd.means = ipd.wt.mean.me, ## ipd.all.wt.mean, weighted means from the 2 ipd combined
#               wt.1 = ipd.1.wts, ## ipd.1.wts, weights for ipd 1
#               ess.1 = ipd.1.ess, ## ipd.1.ess , ess from ipd 1
#               wtd.means.1 = ipd.1.wt.mean, ## ipd.1.wt.mean, weighted means from ipd 1 
#               wt.2 = ipd.2.wts, ## ipd.2.wts, weights for ipd 2 
#               ess.2 = ipd.2.ess, ## ipd.2.ess, ess for ipd 2
#               wtd.means.2 = ipd.2.wt.mean ## ipd.2.wt.mean, weighted means for ipd 2
# >>>>>>> ee91a79fcc65741c26fab3c9af4d94368bb3048b:R/exmWt.2ipd.R
  ))
}