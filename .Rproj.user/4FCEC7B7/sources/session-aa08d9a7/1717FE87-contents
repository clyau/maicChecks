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
#' @references Glimm & Yau (2021). "Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons", Pharmaceutical Statistics, 21(5):974-987. doi:10.1002/pst.2210.
#' @export maicMD
#'
#' @examples
#' \dontrun{
#' ## eAD[1,] is the scenario A in the reference paper,
#' ## i.e. when AD is perfectly within IPD
#' md <- maicMD(eIPD, eAD[1,2:3])
#' md ## a dot-plot of IPD Mahalanobis distances along with AD in the same metric.
#' }
#
# mahalonobis distance, manually
#
maicMD <- function(ipd, ad, n.ad = Inf) {
  ##
  ## assume ipd is a dataframe with n row and p coln
  ## ... n = number of subjects, p = number of matching variables
  ## assume ad is a dataframe with 1 row and p coln
  ##
  md <- y <- NULL
  # ipd and ad must have identical coln's in identical order
  # ipd must be a dataframe
  if (is.data.frame(ipd) == FALSE) {
    ipd <- as.data.frame(ipd)
  }
  # ad must be a row vector
  if (is.data.frame(ad) == FALSE) {
    ad <- data.frame(ad)
    if(dim(ad)[1] > 1) {
      ad <- t(ad)
    }
  }
  # center ipd
  zipd <- as.matrix(sweep(ipd, 2, apply(ipd, 2, mean)))
  # covariance of ipd
  zpz1 <- (dim(ipd)[1] - 1) * solve(t(zipd)%*%(zipd))
  # m-dist. for ipd
  md.ipd <- diag(zipd %*% zpz1 %*% t(zipd))
  # center ad
  zad1 <- as.matrix(ad - apply(ipd, 2, mean))
  # m-dist of ad
  if(n.ad == Inf){
    md.ad <- zad1 %*% zpz1 %*% t(zad1)
  } else{
    md.ad <- n.ad / (dim(ipd)[1] + n.ad) *
      zad1 %*% zpz1 %*% t(zad1)
  }
  if(md.ad > max(md.ipd)) {
    md.check <- 0
  ##  cat('\nM-distance check: AD has the largest Mahalanobis distance to the IPD center.\nTherefore, in original scale AD is *outside* of the IPD convex hull.\n\n')
  } else {
    if(md.ad <= max(md.ipd)){
      md.check <- 2
   ##   cat('\nM-distance check: AD does not have the largest Mahalanobis distance to the IPD center.\n\n')
    }
  }
  #
  # plot
  #
  ipd.plot <- data.frame(md = md.ipd,
                         y = rep(1, length(md.ipd)))
  ad.plot <- data.frame(md = md.ad, y = 1)
  colnames(ad.plot) <- c('md', 'y')
  # coordinate ratio
  # coo.r <- max(md.ipd) / 1.5
  #
  dplot <- ggplot(data = ipd.plot,
                  aes(x = md, y = y),
  ) +
    theme_bw() +
    ylab("") +
    xlab('Mahalanobis distance') +
    geom_point(shape = 1, color = "grey60", size = 2) +
    geom_point(data = ad.plot,
               aes(x = md, y = y),
               shape = 16,
               size = 2.5
    ) +
    theme(panel.grid = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) #+
    #coord_fixed(coo.r)
  #
  return(list(##md.ipd = md.ipd,
              ##md.ad = md.ad,
              md.plot = dplot,
              md.check = md.check)
         )
}

