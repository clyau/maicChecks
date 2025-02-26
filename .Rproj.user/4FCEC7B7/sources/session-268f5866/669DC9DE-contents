#' @title Estimates the MAIC weights
#' @description Estimates the MAIC weights for each individual in the IPD. Should only be used after it is ascertained that AD is indeed within the convex hull of IPD.
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
#' ## eAD[1,] is scenario A in the reference manuscript
#' m1 <- maicWt(eIPD, eAD[1,2:3])
maicWt <- function(ipd, ad, max.it = 25) {
  ##
  ## assume ipd is a dataframe with n row and p coln
  ## ... n = number of subjects, p = number of matching variables
  ## assume ad is a dataframe with 1 row and p coln
  ##
  objfn <- function(a1, X){
    sum(exp(X %*% a1))
  }
  gradfn <- function(a1, X){
    colSums(sweep(X, 1, exp(X %*% a1), "*"))
  }
  ## to take care the case when only 1 variable is matched
  ipd <- as.data.frame(ipd)
  x <- ncol(ipd)
  ipd.n <- nrow(ipd)
  #
  X.EM.1 <- sweep(as.matrix(ipd), 2, as.matrix(ad), '-') # , give.warning = FALSE
  #
  # Estimate $\alpha_2$ (See Philippo 2016)
  #
  op2 <- stats::optim(par = rep(0, x),
               fn = objfn,
               gr = gradfn,
               X = X.EM.1,
               method = "BFGS",
               control = list(maxit = max.it)
  )
  a2 <- op2$par
  wt <- exp(X.EM.1 %*% a2) # weights for each subject in IPD
  wt.rs <- (wt / sum(wt)) * ipd.n # rescaled weights
  #
  ipd.ess <- round(sum(wt.rs)^2 / sum(wt.rs^2), 1)
  ipd.wtsumm <- colMeans(ipd * wt.rs)
  #
  return(list(optim.out = op2,
              maic.wt = wt,
              maic.wt.rs = wt.rs,
              ipd.ess = ipd.ess,
              ipd.wtsumm = ipd.wtsumm))
}
