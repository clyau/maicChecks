#' @title Checks if AD is within the convex hull of IPD using Mahalanobis distance
#' @description Should only be used when all matching variables are normally distributed
#'
#' @param ipd a dataframe with n row and p column, where n is number of subjects and p is the number of variables used in matching.
#' @param ad a dataframe with 1 row and p column. The matching variables should be in the same order as that in \code{ipd}. The function does not check this.
#' @param n.ad default is Inf assuming \code{ad} is a fixed (known) quantity with infinit accuracy. In most MAIC applications \code{ad} is only the sample statistics and n.ad is known.
#'
#' @details When AD does not have the largest Mahalanobis distance, in the original scale AD can still be outside of the IPD convex hull. On the other hand, when AD does have the largest Mahalanobis distance, in the original scale, AD is for sure outside the IPD convex hull.
#' @return Prints a message whether AD is furthest away from 0, i.e. IPD center in terms of Mahalanobis distance. Also returns ggplot object for plotting.
#' \item{md.dplot }{dot-plot of AD and IPD in Mahalanobis distance}
#' \item{md.check }{0 = AD has the largest Mahalanobis distance to the IPD center; 2 = otherwise}
#'
#' @references Glimm E and Yau L. (2022). 'Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons.' \emph{Pharmaceutical Statistics}, 21(5):974-987. \doi{10.1002/pst.2210}.
#' @export maicMD
#'
#' @examples
#' \dontrun{
#' ## eIPD contains response columns; subset to matching columns y1, y2.
#' ## eAD[1,] is the scenario A in the reference paper,
#' ## i.e. when AD is perfectly within IPD
#' md <- maicMD(eIPD[, c('y1', 'y2')], eAD[1, c('y1', 'y2')])
#' md ## a dot-plot of IPD Mahalanobis distances along with AD in the same metric.
#' }
#
# mahalonobis distance, manually
#
maicMD <- function(ipd, ad, n.ad = Inf) {
  ##
  ## assume ipd is a dataframe with n row and p coln
  ## ... n = number of subjects, p = number of matching variables
  ## assume ad is a 1-row dataframe / p-vector
  ##
  md <- y <- NULL ## silence R CMD check NOTE for aes()
  ipd <- as.matrix(ipd)
  ad  <- as.numeric(unlist(ad))
  if (length(ad) != ncol(ipd))
    stop('`ad` must have the same number of variables as `ipd`.', call. = FALSE)

  mu <- colMeans(ipd)
  S  <- stats::cov(ipd)

  md.ipd <- stats::mahalanobis(ipd, center = mu, cov = S)
  md.ad  <- stats::mahalanobis(matrix(ad, nrow = 1), center = mu, cov = S)
  if (is.finite(n.ad))
    md.ad <- (n.ad / (nrow(ipd) + n.ad)) * md.ad

  md.check <- if (md.ad > max(md.ipd)) 0L else 2L

  ##
  ## plot
  ##
  ipd.plot <- data.frame(md = md.ipd, y = 1)
  ad.plot  <- data.frame(md = as.numeric(md.ad), y = 1)

  dplot <- ggplot(data = ipd.plot, aes(x = md, y = y)) +
    theme_bw() +
    ylab('') +
    xlab('Mahalanobis distance') +
    geom_point(shape = 1, color = 'grey60', size = 2) +
    geom_point(data = ad.plot,
               aes(x = md, y = y),
               shape = 16,
               size = 2.5) +
    theme(panel.grid = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

  return(list(md.plot  = dplot,
              md.check = md.check))
}

