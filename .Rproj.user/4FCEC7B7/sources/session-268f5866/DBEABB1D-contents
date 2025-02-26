#' @title Exact matching for two IPD's
#'
#' @param ipd1 a dataframe with n row and p column, where n is number of subjects and p is the number of variables used in matching.
#' @param ipd2 the other IPD with the same number of columns
#' @param vars_to_match variables used for matching. if NULL, use all variables.
#' @param cat_vars_to_01 a list of variable names for the categorical variables that need to be converted to indicator variables.
#' @param mean.constrained whether to restrict the weighted means to be within the ranges of observed means. Default is FALSE. When it is TRUE, there is a higher chance of not having a solution.
#' 
#' @details If dummy variables are already created for the categorical variables in the data set, and are present in \code{ipd1} and \code{ipd2}, then \code{cat_vars_to_01} should be left as NULL.
#' 
#' @return
#' \item{ipd1 }{re-scaled weights of the exact matching by maximizing ESS for IPD 1, and the input IPD 1 data with categorical variables converted to 0-1 indicators}
#' \item{ipd2 }{re-scaled weights of the exact matching by maximizing ESS for IPD 2, and the input IPD 2 data with categorical variables converted to 0-1 indicators}
#' \item{wtd.summ }{ESS for IPD 1, ESS for IPD 2, and weighted means of the matching variables}
#' 
#' @export exmWt.2ipd
#' 
#' @author Lillian Yau
#' 
exmWt.2ipd <- function (ipd1, ipd2, vars_to_match = NULL, cat_vars_to_01 = NULL, mean.constrained = FALSE) 
{
  ## check vars_to_match
  vars_to_match <- .check_data(ipd1, 
                              ipd2, 
                              v.ars_to_match = vars_to_match,
                              c.at_vars_to_01 = cat_vars_to_01
  )
  
  ## keep variables not used
  ipd1_not_used <- ipd1[!(colnames(ipd1) %in% vars_to_match)]
  ipd2_not_used <- ipd2[!(colnames(ipd1) %in% vars_to_match)]
  
  ## extract only vars_to_match from both ipd's
  ipd1 <- ipd1[vars_to_match]
  ipd2 <- ipd2[vars_to_match]
  
  ## save original input data and create new ones
  if(!is.null(cat_vars_to_01)){
    ipd1.o <- ipd1 ## save original input data for later use
    ipd2.o <- ipd2
    ipd1 <- .cat201_minus1(ipd1.o, v.cat = cat_vars_to_01)
    ipd2 <- .cat201_minus1(ipd2.o, v.cat = cat_vars_to_01)
  }
  
  ##
  ## derivation for weights calculation starts here ::::
  ##
  ipd <- as.data.frame(rbind(-1 * ipd1, ipd2))
  oneszeros <- c(rep(1, nrow(ipd1)), rep(0, nrow(ipd2)))
  zerosones <- c(rep(0, nrow(ipd1)), rep(1, nrow(ipd2)))
  ipd <- as.data.frame(cbind(ipd, oneszeros, zerosones))
  rm(oneszeros, zerosones)
  ipd.n <- nrow(ipd)
  p <- ncol(ipd1)
  bvec <- matrix(c(rep(0, p), 
                   1, 
                   1, 
                   rep(0, nrow(ipd))), 
                 nrow = 1)
  ad0 <- bvec
  Amat <- as.matrix(ipd)
  Dmat <- diag(ipd.n)
  Amat0 <- as.matrix(data.frame(cbind(Amat, Dmat)))
  dvec <- rep(0, ipd.n)
  
  if (mean.constrained) { ## if constrained means == TRUE
    if(!is.null(cat_vars_to_01)){
      ipd1 <- .cat201(ipd1.o, v.cat = cat_vars_to_01)
      ipd2 <- .cat201(ipd2.o, v.cat = cat_vars_to_01)
    }
    ipd1.bar <- colMeans(ipd1) 
    ipd2.bar <- colMeans(ipd2)
    x <- as.data.frame(rbind(ipd1.bar, ipd2.bar))
    bar.min <- apply(x, 2, min)
    bar.max <- apply(x, 2, max)
    bvec <- matrix(c(ad0, 2 * bar.min, 2 * bar.max * (-1)), 
                   nrow = 1)
    x0 <- as.data.frame(rbind(ipd1, ipd2))
    x1 <- as.data.frame(rbind(-1 * ipd1, -1 * ipd2))
    Amat0 <- as.matrix(data.frame(cbind(Amat, Dmat, x0, x1)))
  }
  ##
  wts <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat0, 
                            bvec = bvec, meq = p + 2)
  
  ## save results
  ## weights
  ipd.wts.me <- wts[["solution"]]
  ipd.1.wts.u <- ipd.wts.me[1:nrow(ipd1)]
  ipd.2.wts.u <- ipd.wts.me[(1 + nrow(ipd1)):(nrow(ipd1) + nrow(ipd2))]
  
  ## re-scale weights
  ipd.1.wts <- ipd.1.wts.u * nrow(ipd1)
  ipd.2.wts <- ipd.2.wts.u * nrow(ipd2)
  
  ## weighted covariates means, only need one since they are the same
  ipd.wtd.mean <- colMeans(ipd1 * ipd.1.wts)
  ##ipd2.wtd.mean <- colMeans(ipd2 * ipd.2.wts)
  
  ## ess
  ipd1.ess <- round(sum(ipd.1.wts)^2/sum(ipd.1.wts^2), 1)
  ipd2.ess <- round(sum(ipd.2.wts)^2/sum(ipd.2.wts^2), 1)
  
  ## combine re-scaled weights with ipd1 data
  ipd1 <- data.frame(exm.wts = ipd.1.wts, ipd1)
  ## combine with un-used ipd1 variables
  if(!is.null(ipd1_not_used))
    ipd1 <- data.frame(cbind(ipd1_not_used, ipd1))
  
  ## combine re-scaled weights with ipd2 data
  ipd2 <- data.frame(exm.wts = ipd.2.wts, ipd2)
  ## combined with un-used ipd2 variables
  if(!is.null(ipd2_not_used))
    ipd2 <- data.frame(cbind(ipd2_not_used, ipd2))
  
  ## summary (weighted means) and ess
  wtd.summ <- rbind(c(ipd1.ess = ipd1.ess, 
                      ipd2.ess = ipd2.ess, 
                      ipd.wtd.mean))

  return(list(ipd1 = ipd1,
              ipd2 = ipd2,
              wtd.summ = wtd.summ))
}