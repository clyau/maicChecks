#' @title Weighted treatment-effect difference with Wald confidence interval
#'
#' @description
#' Computes the weighted mean response in each of two arms (IPD vs IPD or
#' IPD vs AD), their difference, and a Wald (normal-based) confidence
#' interval for the difference using the variance estimator described in
#' Section 5 of Glimm & Yau (2026). Two modes are supported:
#' \itemize{
#'   \item \strong{IPD vs IPD} -- supply per-subject treatment-effect
#'     responses and weights for both arms
#'     (\code{ipd1.te, w1, ipd2.te, w2}). Weights are typically obtained
#'     from \code{\link{exmWt.2ipd}}.
#'   \item \strong{IPD vs AD} -- supply per-subject responses and weights
#'     for one arm (\code{ipd1.te, w1}) together with an aggregate summary
#'     for the other (\code{ad.mean}, optionally \code{ad.sd} and
#'     \code{ad.n}). IPD weights are typically obtained from
#'     \code{\link{maicWt}} or \code{\link{maxessWt}}.
#' }
#' The function is intended for \strong{binary (0/1) or continuous}
#' responses only. It is \emph{not} suitable for time-to-event / survival
#' endpoints. This is not checked by the function -- the user is responsible
#' for supplying an appropriate response type.
#'
#' @details
#' The weighted point estimate in arm \eqn{k} is
#' \deqn{\hat{\mu}_k = \frac{\sum_i w_i y_i}{\sum_i w_i},}
#' and following Glimm & Yau (2026, Section 5) its conditional variance is
#' estimated as
#' \deqn{\widehat{\mathrm{var}}(\hat{\mu}_k \mid X^{(k)}) =
#'       \frac{\sum_i w_i^2}{(\sum_i w_i)^2}\, s_k^2 = \frac{s_k^2}{\mathrm{ESS}_k},}
#' where \eqn{s_k^2 = \frac{1}{n_k}\sum_i (y_i - \bar{y}_k)^2} uses the
#' \emph{unweighted} sample mean \eqn{\bar{y}_k}, and
#' \eqn{\mathrm{ESS}_k = (\sum_i w_i)^2 / \sum_i w_i^2} is the effective
#' sample size. Under the assumptions stated in Section 5 (no hidden
#' confounders, within-study exchangeability, and
#' \eqn{\mathrm{var}(Y) \geq \mathrm{var}(Y \mid X)}), this is a
#' \emph{conservative} estimator of the conditional variance.
#'
#' For the \strong{IPD vs AD} mode, if both \code{ad.sd} and \code{ad.n} are
#' supplied the AD variance is \eqn{\hat{\sigma}_{ad}^2 / n_{ad}}; otherwise
#' the AD mean is treated as a fixed constant
#' (\eqn{\widehat{\mathrm{var}} = 0}). For binary AD endpoints the user may
#' supply \code{ad.sd = sqrt(p * (1 - p))}.
#'
#' Assuming the two arms are independent (e.g. they come from different
#' studies),
#' \deqn{\mathrm{SE}(\hat{\Delta}) = \sqrt{\widehat{\mathrm{var}}(\hat{\mu}_1)
#'        + \widehat{\mathrm{var}}(\hat{\mu}_2)},}
#' and the \eqn{(1-\alpha)} Wald confidence interval for
#' \eqn{\hat{\Delta} = \hat{\mu}_1 - \hat{\mu}_2} is
#' \eqn{\hat{\Delta} \pm z_{1-\alpha/2} \cdot \mathrm{SE}(\hat{\Delta})}.
#'
#' @param ipd1.te numeric vector of per-subject treatment-effect responses for IPD 1 (0/1 or continuous). The suffix \code{.te} stands for \emph{treatment effect}.
#' @param w1 numeric vector of weights for IPD 1 (e.g. \code{maic.wt.rs} from \code{\link{maicWt}}, \code{maxess.wt} from \code{\link{maxessWt}}, or \code{exm.wts} from \code{\link{exmWt.2ipd}}). Must have the same length as \code{ipd1.te}.
#' @param ipd2.te optional numeric vector of per-subject treatment-effect responses for IPD 2.
#' @param w2 optional numeric vector of weights for IPD 2. Must have the same length as \code{ipd2.te}.
#' @param ad.mean optional scalar, the aggregate-data (AD) mean response. Used only when \code{ipd2.te} and \code{w2} are not supplied.
#' @param ad.sd optional scalar, the AD standard deviation of the response. If \code{NULL} the AD mean is treated as a constant.
#' @param ad.n optional scalar, the AD sample size.
#' @param conf.level confidence level for the Wald CI. Default \code{0.95}.
#'
#' @return A list with the following slots:
#' \item{wt.y1}{weighted mean response in arm 1.}
#' \item{wt.y2}{weighted mean response in arm 2 (or \code{ad.mean}).}
#' \item{wt.diff}{the difference \code{wt.y1 - wt.y2}.}
#' \item{se1}{standard error of \code{wt.y1}.}
#' \item{se2}{standard error of \code{wt.y2} (\code{0} when AD is treated as a constant).}
#' \item{se.diff}{standard error of \code{wt.diff}.}
#' \item{ci.lower}{lower Wald confidence limit for \code{wt.diff}.}
#' \item{ci.upper}{upper Wald confidence limit for \code{wt.diff}.}
#' \item{conf.level}{the confidence level used.}
#' \item{ess1}{effective sample size for arm 1.}
#' \item{ess2}{effective sample size for arm 2 (\code{ad.n} when supplied, \code{NA} otherwise).}
#'
#' @references
#' Glimm E and Yau L. (2026). 'Exact matching as an alternative to propensity score matching.' \emph{Statistics in Biopharmaceutical Research}, 18(1):106-116. \doi{10.1080/19466315.2025.2507378}.
#'
#' @export wtTrtDiff
#'
#' @author Lillian Yau
#'
#' @examples
#' ## ------------------------------------------------------------------
#' ## IPD vs AD: weight eIPD onto scenario A of eAD with maicWt, then
#' ## compare the IPD continuous / binary responses (eIPD$r.cont,
#' ## eIPD$r.bin) to the AD summaries stored in the same scenario row of
#' ## eAD. Note y1, y2 in eIPD / eAD are the *matching covariates*, while
#' ## r.cont, r.bin (and r.cont.mean/sd/n, r.bin.p/n in eAD) are the
#' ## response data.
#' ## ------------------------------------------------------------------
#' w.out <- maicWt(eIPD[, c('y1', 'y2')], eAD[1, c('y1', 'y2')])
#'
#' ## continuous response
#' wtTrtDiff(ipd1.te = eIPD$r.cont, w1 = w.out$maic.wt.rs,
#'           ad.mean = eAD$r.cont.mean[1],
#'           ad.sd   = eAD$r.cont.sd[1],
#'           ad.n    = eAD$r.cont.n[1])
#'
#' ## binary response, treating AD as a fixed constant
#' wtTrtDiff(ipd1.te = eIPD$r.bin, w1 = w.out$maic.wt.rs,
#'           ad.mean = eAD$r.bin.p[1])
#'
#' ## binary response, supplying ad.sd = sqrt(p * (1 - p))
#' wtTrtDiff(ipd1.te = eIPD$r.bin, w1 = w.out$maic.wt.rs,
#'           ad.mean = eAD$r.bin.p[1],
#'           ad.sd   = sqrt(eAD$r.bin.p[1] * (1 - eAD$r.bin.p[1])),
#'           ad.n    = eAD$r.bin.n[1])
#'
#' \dontrun{
#' ## ------------------------------------------------------------------
#' ## IPD vs IPD: symmetric exact-matching weights on sim110 IPD A vs B,
#' ## then compare the simulated continuous (Y) and binary (Y.bin)
#' ## responses.
#' ## ------------------------------------------------------------------
#' ipd1  <- sim110[sim110$study == 'IPD A', ]
#' ipd2  <- sim110[sim110$study == 'IPD B', ]
#' w.out <- exmWt.2ipd(ipd1, ipd2,
#'                     vars_to_match  = paste0('X', 1:5),
#'                     cat_vars_to_01 = paste0('X', 1:3))
#' ## continuous response
#' wtTrtDiff(ipd1.te = ipd1$Y,     w1 = w.out$ipd1$exm.wts,
#'           ipd2.te = ipd2$Y,     w2 = w.out$ipd2$exm.wts)
#' ## binary response
#' wtTrtDiff(ipd1.te = ipd1$Y.bin, w1 = w.out$ipd1$exm.wts,
#'           ipd2.te = ipd2$Y.bin, w2 = w.out$ipd2$exm.wts)
#' }
wtTrtDiff <- function(ipd1.te, w1,
                      ipd2.te = NULL, w2 = NULL,
                      ad.mean = NULL, ad.sd = NULL, ad.n = NULL,
                      conf.level = 0.95) {

  ## ---- validate IPD-1 inputs --------------------------------------------
  if (!is.numeric(ipd1.te) || !is.numeric(w1))
    stop('`ipd1.te` and `w1` must be numeric vectors.', call. = FALSE)
  if (length(ipd1.te) != length(w1))
    stop('`ipd1.te` and `w1` must have the same length.', call. = FALSE)
  if (length(ipd1.te) < 2L)
    stop('`ipd1.te` must have length >= 2.', call. = FALSE)
  if (anyNA(ipd1.te) || anyNA(w1))
    stop('`ipd1.te` and `w1` must not contain NA.', call. = FALSE)
  if (any(w1 < 0))
    stop('weights in `w1` must be non-negative.', call. = FALSE)
  if (sum(w1) <= 0)
    stop('weights in `w1` must sum to a positive number.', call. = FALSE)
  if (length(conf.level) != 1L || !is.numeric(conf.level) ||
      conf.level <= 0 || conf.level >= 1)
    stop('`conf.level` must be a single number in (0, 1).', call. = FALSE)

  ## ---- decide mode ------------------------------------------------------
  have.ipd2 <- !is.null(ipd2.te) || !is.null(w2)
  have.ad   <- !is.null(ad.mean)

  if (have.ipd2 && have.ad)
    stop('Specify either `ipd2.te`/`w2` (IPD vs IPD) or `ad.mean` (IPD vs AD), not both.',
         call. = FALSE)
  if (!have.ipd2 && !have.ad)
    stop('Either `ipd2.te` and `w2` (IPD vs IPD) or `ad.mean` (IPD vs AD) must be supplied.',
         call. = FALSE)

  ## ---- IPD 1: weighted mean and Section-5 variance ----------------------
  sumw1  <- sum(w1)
  sumw1s <- sum(w1^2)
  wt.y1  <- sum(w1 * ipd1.te) / sumw1
  ess1   <- sumw1^2 / sumw1s
  s1.sq  <- mean((ipd1.te - mean(ipd1.te))^2)   ## unweighted mean; divisor n1
  var1   <- (sumw1s / sumw1^2) * s1.sq          ## == s1.sq / ess1
  se1    <- sqrt(var1)

  ## ---- IPD 2 ------------------------------------------------------------
  if (have.ipd2) {
    if (is.null(ipd2.te) || is.null(w2))
      stop('Both `ipd2.te` and `w2` must be supplied for IPD-vs-IPD mode.',
           call. = FALSE)
    if (!is.numeric(ipd2.te) || !is.numeric(w2))
      stop('`ipd2.te` and `w2` must be numeric vectors.', call. = FALSE)
    if (length(ipd2.te) != length(w2))
      stop('`ipd2.te` and `w2` must have the same length.', call. = FALSE)
    if (length(ipd2.te) < 2L)
      stop('`ipd2.te` must have length >= 2.', call. = FALSE)
    if (anyNA(ipd2.te) || anyNA(w2))
      stop('`ipd2.te` and `w2` must not contain NA.', call. = FALSE)
    if (any(w2 < 0))
      stop('weights in `w2` must be non-negative.', call. = FALSE)
    if (sum(w2) <= 0)
      stop('weights in `w2` must sum to a positive number.', call. = FALSE)

    sumw2  <- sum(w2)
    sumw2s <- sum(w2^2)
    wt.y2  <- sum(w2 * ipd2.te) / sumw2
    ess2   <- sumw2^2 / sumw2s
    s2.sq  <- mean((ipd2.te - mean(ipd2.te))^2)
    var2   <- (sumw2s / sumw2^2) * s2.sq
    se2    <- sqrt(var2)
  } else {
    if (!is.numeric(ad.mean) || length(ad.mean) != 1L)
      stop('`ad.mean` must be a single number.', call. = FALSE)
    wt.y2 <- ad.mean

    have.ad.sd <- !is.null(ad.sd)
    have.ad.n  <- !is.null(ad.n)
    if (have.ad.sd != have.ad.n)
      stop('`ad.sd` and `ad.n` must be supplied together (or both omitted).',
           call. = FALSE)

    if (have.ad.sd) {
      if (!is.numeric(ad.sd) || length(ad.sd) != 1L || ad.sd < 0)
        stop('`ad.sd` must be a single non-negative number.', call. = FALSE)
      if (!is.numeric(ad.n) || length(ad.n) != 1L || ad.n <= 0)
        stop('`ad.n` must be a single positive number.', call. = FALSE)
      var2 <- ad.sd^2 / ad.n
      ess2 <- ad.n
    } else {
      var2 <- 0
      ess2 <- NA_real_
    }
    se2  <- sqrt(var2)
  }

  ## ---- difference and Wald CI -------------------------------------------
  wt.diff <- wt.y1 - wt.y2
  se.diff <- sqrt(var1 + var2)
  alpha   <- 1 - conf.level
  z       <- stats::qnorm(1 - alpha / 2)
  ci.lo   <- wt.diff - z * se.diff
  ci.hi   <- wt.diff + z * se.diff

  return(list(wt.y1      = wt.y1,
              wt.y2      = wt.y2,
              wt.diff    = wt.diff,
              se1        = se1,
              se2        = se2,
              se.diff    = se.diff,
              ci.lower   = ci.lo,
              ci.upper   = ci.hi,
              conf.level = conf.level,
              ess1       = ess1,
              ess2       = ess2))
}
