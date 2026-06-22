#' an IPD set
#'
#' An artificial individual-patient-data (IPD) set used in Glimm & Yau (2022).
#' Columns \code{y1} and \code{y2} are the \strong{matching covariates}
#' (despite the \code{y} prefix they play the role of
#' baseline covariates X in MAIC / exact matching). Columns \code{r.cont}
#' and \code{r.bin} are simulated \strong{response variables} added in
#' v0.3.0 for use with \code{\link{wtTrtDiff}}.
#'
#' @format
#' \describe{
#' \item{\code{y1}}{numeric, matching covariate 1.}
#' \item{\code{y2}}{numeric, matching covariate 2.}
#' \item{\code{r.cont}}{numeric, simulated \emph{continuous} response
#'   (\code{rnorm} with mean 1.0 and SD 0.5).}
#' \item{\code{r.bin}}{integer 0/1, simulated \emph{binary} response
#'   (\code{rbinom} with prob 0.6).}
#' }
#'
#' @docType data
#'
#' @usage data(eIPD)
#'
#' @keywords datasets
#'
#' @references Glimm E and Yau L. (2022). 'Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons.' \emph{Pharmaceutical Statistics}, 21(5):974-987. \doi{10.1002/pst.2210}.
#'
#'
#' @examples
#' data(eIPD)
"eIPD"
