#' @title Exact matching for two IPD's
#'
#' @description
#' Computes exact-matching weights for two individual-patient datasets
#' (IPDs). Two modes are supported:
#' \itemize{
#'   \item \strong{Default (\code{target.ipd = NULL})} -- the method described
#'     in Glimm & Yau (2026): both IPDs are reweighted simultaneously so that
#'     their weighted covariate means match each other and the joint
#'     effective sample size (ESS) is maximised. The optimisation is solved
#'     via \code{\link[quadprog]{solve.QP}}.
#'   \item \strong{One-sided MAIC (\code{target.ipd = 'ipd1'} or
#'     \code{'ipd2'})} -- the chosen IPD plays the role of the aggregate
#'     reference: its column means (after dummy expansion) are used as the
#'     AD, and the other IPD is weighted with either \code{\link{maicWt}}
#'     (default) or \code{\link{maxessWt}}. The target IPD is returned with
#'     uniform weights of 1.
#' }
#' In all modes the function first calls \code{\link{exmLP.2ipd}} as a
#' feasibility gate. If matching is infeasible, a \code{message()} is
#' emitted and \code{NA} weights are returned together with
#' \code{lp.check = 2}, instead of letting the QP / optimiser fail.
#'
#' @param ipd1 a dataframe with n1 row and p column.
#' @param ipd2 a dataframe with n2 row and p column.
#' @param vars_to_match variables used for matching. If \code{NULL} (default), all variables shared by \code{ipd1} and \code{ipd2} are used.
#' @param cat_vars_to_01 a vector of variable names for the categorical variables that need to be converted to indicator variables. If the dummies are already present in the data, leave as \code{NULL}.
#' @param mean.constrained whether to restrict the weighted means to be within the ranges of observed means. Only used when \code{target.ipd = NULL}. Default is \code{FALSE}. When \code{TRUE}, the QP is more likely to have no solution.
#' @param target.ipd one of \code{NULL} (default), \code{'ipd1'}, or \code{'ipd2'}. When \code{NULL}, the default method (Glimm & Yau, 2026) is used. When \code{'ipd1'}, the means of \code{ipd1} are used as the AD and \code{ipd2} is weighted via \code{method}; vice versa for \code{'ipd2'}.
#' @param method one of \code{'maicWt'} (default) or \code{'maxessWt'}. Ignored when \code{target.ipd = NULL}.
#'
#' @details If dummy variables are already created for the categorical variables in the data set, and are present in \code{ipd1} and \code{ipd2}, then \code{cat_vars_to_01} should be left as \code{NULL}.
#'
#' @return A list with the following slots:
#' \item{ipd1}{re-scaled exact-matching weights for IPD 1 (column \code{exm.wts}) together with the IPD 1 data, categorical variables converted to 0/1 indicators. \code{NA} if the matching is infeasible.}
#' \item{ipd2}{same for IPD 2.}
#' \item{wtd.summ}{a 1-row matrix with ESS for IPD 1, ESS for IPD 2, and the weighted means of the matching variables.}
#' \item{lp.check}{the \code{\link{exmLP.2ipd}} (or \code{\link{maicLP}}) feasibility status: \code{0} = feasible, \code{2} = infeasible.}
#' \item{target.ipd}{the value of \code{target.ipd} that was used.}
#' \item{method}{the weighting method used: \code{'qp'} for the default Glimm & Yau (2026) method, otherwise \code{'maicWt'} or \code{'maxessWt'}.}
#'
#' @export exmWt.2ipd
#'
#' @author Lillian Yau
#'
#' @references Glimm E and Yau L. (2026). 'Exact matching as an alternative to propensity score matching.' \emph{Statistics in Biopharmaceutical Research}, 18(1):106-116. \doi{10.1080/19466315.2025.2507378}.
#'
#' @examples
#' \dontrun{
#' ipd1 <- sim110[sim110$study == 'IPD A', ]
#' ipd2 <- sim110[sim110$study == 'IPD B', ]
#'
#' ## (1) Default: both IPDs reweighted simultaneously (Glimm & Yau, 2026)
#' x0 <- exmWt.2ipd(ipd1, ipd2,
#'                  vars_to_match    = paste0('X', 1:5),
#'                  cat_vars_to_01   = paste0('X', 1:3),
#'                  mean.constrained = FALSE)
#'
#' ## (2) One-sided MAIC: target = IPD A, weight IPD B via maicWt
#' x1 <- exmWt.2ipd(ipd1, ipd2,
#'                  vars_to_match  = paste0('X', 1:5),
#'                  cat_vars_to_01 = paste0('X', 1:3),
#'                  target.ipd     = 'ipd1',
#'                  method         = 'maicWt')
#'
#' ## (3) One-sided MAIC: target = IPD B, weight IPD A via maxessWt
#' x2 <- exmWt.2ipd(ipd1, ipd2,
#'                  vars_to_match  = paste0('X', 1:5),
#'                  cat_vars_to_01 = paste0('X', 1:3),
#'                  target.ipd     = 'ipd2',
#'                  method         = 'maxessWt')
#' }
exmWt.2ipd <- function(ipd1, ipd2,
                       vars_to_match    = NULL,
                       cat_vars_to_01   = NULL,
                       mean.constrained = FALSE,
                       target.ipd       = NULL,
                       method           = c('maicWt', 'maxessWt')) {

  method <- match.arg(method)
  if (!is.null(target.ipd) && !(target.ipd %in% c('ipd1', 'ipd2')))
    stop("`target.ipd` must be NULL, 'ipd1', or 'ipd2'.", call. = FALSE)

  ## ---- validate inputs ----------------------------------------------------
  vars_to_match <- .check_data(ipd1, ipd2,
                               v.ars_to_match  = vars_to_match,
                               c.at_vars_to_01 = cat_vars_to_01)

  ## keep variables not used (so we can return the original IPDs with weights)
  ipd1_not_used <- ipd1[!(colnames(ipd1) %in% vars_to_match)]
  ipd2_not_used <- ipd2[!(colnames(ipd2) %in% vars_to_match)]

  ipd1 <- data.frame(ipd1[vars_to_match])
  ipd2 <- data.frame(ipd2[vars_to_match])

  ## expand categorical variables to k-1 dummies
  if (!is.null(cat_vars_to_01)) {
    ipd1.o <- ipd1
    ipd2.o <- ipd2
    ipd1 <- .cat201_minus1(ipd1.o, v.cat = cat_vars_to_01)
    ipd2 <- .cat201_minus1(ipd2.o, v.cat = cat_vars_to_01)
  }

  n1 <- nrow(ipd1)
  n2 <- nrow(ipd2)

  ## ---- helper: package up an infeasible-return list -----------------------
  .infeasible_return <- function(status) {
    message('matching not feasible')
    list(ipd1       = NA,
         ipd2       = NA,
         wtd.summ   = NA,
         lp.check   = status,
         target.ipd = target.ipd,
         method     = if (is.null(target.ipd)) 'qp' else method)
  }

  ## =========================================================================
  ##  Branch 1: one-sided MAIC (target.ipd is 'ipd1' or 'ipd2')
  ## =========================================================================
  if (!is.null(target.ipd)) {

    if (target.ipd == 'ipd1') {
      target.mean <- colMeans(ipd1)
      other       <- ipd2
    } else {
      target.mean <- colMeans(ipd2)
      other       <- ipd1
    }

    ## feasibility gate (maicLP on the other IPD vs target means)
    lp.chk <- maicLP(other, target.mean)$lp.check
    if (lp.chk != 0) return(.infeasible_return(lp.chk))

    ## delegate to the chosen MAIC weighting routine
    other.wts <- switch(method,
                        maicWt   = maicWt  (other, target.mean)$maic.wt.rs,
                        maxessWt = maxessWt(other, target.mean)$maxess.wt)
    other.wts <- as.numeric(other.wts)

    ## assemble weights: target gets uniform 1's, other gets the MAIC weights
    if (target.ipd == 'ipd1') {
      ipd.1.wts <- rep(1, n1)
      ipd.2.wts <- other.wts
    } else {
      ipd.1.wts <- other.wts
      ipd.2.wts <- rep(1, n2)
    }

    ipd1.ess <- round(sum(ipd.1.wts)^2 / sum(ipd.1.wts^2), 1)
    ipd2.ess <- round(sum(ipd.2.wts)^2 / sum(ipd.2.wts^2), 1)
    ## weighted means (will equal target.mean when feasible)
    ipd.wtd.mean <- if (target.ipd == 'ipd1') colMeans(ipd2 * ipd.2.wts)
                    else                       colMeans(ipd1 * ipd.1.wts)

  } else {
    ## =======================================================================
    ##  Branch 2: legacy symmetric QP (target.ipd = NULL)
    ## =======================================================================

    ## pre-feasibility check
    lp.chk <- exmLP.2ipd(ipd1            = if (!is.null(cat_vars_to_01)) ipd1.o else ipd1,
                         ipd2            = if (!is.null(cat_vars_to_01)) ipd2.o else ipd2,
                         vars_to_match   = vars_to_match,
                         cat_vars_to_01  = cat_vars_to_01,
                         mean.constrained = mean.constrained)$lp.check
    if (lp.chk != 0) return(.infeasible_return(lp.chk))

    ## ---- build QP -------------------------------------------------------
    ipd <- as.data.frame(rbind(-1 * ipd1, ipd2))
    ipd <- as.data.frame(cbind(ipd,
                               oneszeros = c(rep(1, n1), rep(0, n2)),
                               zerosones = c(rep(0, n1), rep(1, n2))))
    p     <- ncol(ipd1)
    ipd.n <- nrow(ipd)
    bvec  <- matrix(c(rep(0, p), 1, 1, rep(0, ipd.n)), nrow = 1)
    ad0   <- bvec
    Amat  <- as.matrix(ipd)
    Dmat  <- diag(ipd.n)
    Amat0 <- as.matrix(data.frame(cbind(Amat, Dmat)))
    dvec  <- rep(0, ipd.n)

    if (isTRUE(mean.constrained)) {
      if (!is.null(cat_vars_to_01)) {
        ipd1 <- .cat201(ipd1.o, v.cat = cat_vars_to_01)
        ipd2 <- .cat201(ipd2.o, v.cat = cat_vars_to_01)
      }
      ipd1.bar <- colMeans(ipd1)
      ipd2.bar <- colMeans(ipd2)
      x        <- as.data.frame(rbind(ipd1.bar, ipd2.bar))
      bar.min  <- apply(x, 2, min)
      bar.max  <- apply(x, 2, max)
      bvec     <- matrix(c(ad0, 2 * bar.min, 2 * bar.max * (-1)), nrow = 1)
      x0       <- as.data.frame(rbind(ipd1, ipd2))
      x1       <- as.data.frame(rbind(-1 * ipd1, -1 * ipd2))
      Amat0    <- as.matrix(data.frame(cbind(Amat, Dmat, x0, x1)))
    }

    wts <- quadprog::solve.QP(Dmat = Dmat,
                              dvec = dvec,
                              Amat = Amat0,
                              bvec = bvec,
                              meq  = p + 2)

    ipd.wts.me  <- wts[['solution']]
    ipd.1.wts.u <- ipd.wts.me[seq_len(n1)]
    ipd.2.wts.u <- ipd.wts.me[n1 + seq_len(n2)]

    ## re-scale weights to sum to the per-IPD sample size
    ipd.1.wts <- ipd.1.wts.u * n1
    ipd.2.wts <- ipd.2.wts.u * n2

    ## the symmetric QP equates the two weighted means
    ipd.wtd.mean <- colMeans(ipd1 * ipd.1.wts)

    ipd1.ess <- round(sum(ipd.1.wts)^2 / sum(ipd.1.wts^2), 1)
    ipd2.ess <- round(sum(ipd.2.wts)^2 / sum(ipd.2.wts^2), 1)
  }

  ## ---- assemble return objects -------------------------------------------
  ipd1.out <- data.frame(exm.wts = ipd.1.wts, ipd1)
  if (ncol(ipd1_not_used) > 0)
    ipd1.out <- data.frame(cbind(ipd1_not_used, ipd1.out))

  ipd2.out <- data.frame(exm.wts = ipd.2.wts, ipd2)
  if (ncol(ipd2_not_used) > 0)
    ipd2.out <- data.frame(cbind(ipd2_not_used, ipd2.out))

  wtd.summ <- rbind(c(ipd1.ess = ipd1.ess,
                      ipd2.ess = ipd2.ess,
                      ipd.wtd.mean))

  return(list(ipd1       = ipd1.out,
              ipd2       = ipd2.out,
              wtd.summ   = wtd.summ,
              lp.check   = lp.chk,
              target.ipd = target.ipd,
              method     = if (is.null(target.ipd)) 'qp' else method))
}
