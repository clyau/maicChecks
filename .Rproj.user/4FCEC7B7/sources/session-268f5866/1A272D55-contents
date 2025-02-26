#' @title  Checks whether AD is outside IPD in PC coordinates
#' @description Checks whether AD is outside IPD in principal component (PC) coordinates
#'
#' @param ipd a dataframe with n row and p column, where n is number of subjects in IPD set and p is the number of variables used in matching.
#' @param ad a dataframe with 1 row and p column. The matching variables should be in the same order as that in \code{ipd}. The function does not check this.
#'
#' @details When AD is within the IPD PC ranges, AD can still be outside the IPD convex hull in the original scale. On the other hand, if AD is outside the IPD PC ranges, in the original scale AD is for sure outside the IPD convex hull.
#' @return Prints a message whether AD is inside or outside IPD PC coordinates. Also returns a ggplot object to be plotted.
#' \item{pc.dplot }{dot-plot of AD and IPD both in IPD's PC coordinates}
#' \item{pca.check }{0 = AD within the ranges of IPD's PC coordinates; 2 = otherwise}
#'
#' @export maicPCA
#'
#' @references Glimm & Yau (2021). "Geometric approaches to assessing the numerical feasibility for conducting matching-adjusted indirect comparisons", Pharmaceutical Statistics, 21(5):974-987. doi:10.1002/pst.2210.
#' @examples
#' \dontrun{
#' ## eAD[1,] is the scenario A in the reference paper,
#' ## i.e. when AD is perfectly within IPD
#' a1 <- maicPCA(eIPD, eAD[1,2:3])
#' a1 ## the dot plots of PC's for IPD and AD
#'
#' ## eAD[3,] is the scenario C in the reference paper,
#' ## i.e. when AD is outside IPD
#' a3 <- maicPCA(eIPD, eAD[3,2:3])
#' a3 ## the dot plots of PC's for IPD and AD
#' }
maicPCA <- function (ipd, ad) {
  ##
  ## assume ipd is a dataframe with n row and p coln
  ## ... n = number of subjects, p = number of matching variables
  ## assume ad is a dataframe with 1 row and p coln
  ##
  pc.o <- NULL
  ##
  ## standardize ipd and ad
  ##
  ipdm <- apply(ipd, 2, mean)
  ipdsd <- apply(ipd, 2, sd)
  zi1 <- sweep(ipd, 2, ipdm) # scaling done in prcomp()
  w <- (ad - ipdm) / ipdsd # standardize ad w.r.t. ipd
  ##
  ## pca
  ##
  pc <- stats::prcomp(zi1, retx = TRUE, scale. = TRUE)
  pc.ipd <- pc$x
  ## ad in ipd's pca coordinates
  pc.ad <- data.frame(t(pc$rotation) %*% t(w))
  ##
  ## check if ad is within all ipd's pc coordinates
  ##
  x <- data.frame(t(apply(pc.ipd, 2, range)))
  colnames(x) <- c('min', 'max')
  if (all(data.table::between(t(pc.ad), x$min, x$max))) {
    pc.check <- 0 ## pca.check <- 0 ... return asks for pc.check
  ##  cat("PCA check: AD within the ranges of IPD's PC coordinates.\n")
  }
  else {
    pc.check <- 2 ## pca.check <- 2 ... return asks for pc.check
  ##  cat("PCA check: AD outside the ranges of IPD's PC coordinates.\n")
  }
  ##
  ## create plot
  ##
  pc.ipd.long <- tidyr::gather(as.data.frame(pc.ipd), pc, x)
  pc.ipd.long$pc.o <- rep(1:dim(ipd)[2], each = dim(ipd)[1])
  ## ad in pca scale
  pc.ad$pc.o <- 1:nrow(pc.ad)
  colnames(pc.ad) <- c("w", "pc.o")
  ##
  ## plot
  ##
  pc.dplot <-
    ggplot(pc.ipd.long,
           aes(x, factor(pc.o))) +
    geom_point(shape = 1, color = "grey60", size = 2) +
    ## vertical line to go thru 0, i.e. ipd centers in pc coordinates
    geom_vline(xintercept = 0,
               color = "gray25",
               size = .5,
               alpha = .5,
               linetype = 'dashed') +
    geom_point(data = pc.ad,
               mapping = aes(w, factor(pc.o)),
               color = "black",
               size = 2.5,
               shape = 17) +
    theme_bw(base_size = 10) +
    scale_y_discrete(breaks = 1:max(pc.ipd.long$pc.o),
                     labels = paste0("PC", 1:max(pc.ipd.long$pc.o))) +
    ylab("") +
    geom_hline(yintercept = seq(1.5,
                                max(pc.ipd.long$pc.o),
                                by = 1),
               color = "gray",
               size = 0.5,
               alpha = 0.5
    ) +
    guides(size = guide_legend("# of obs.")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor.x = element_blank()) +
    xlab("IPD PC values")
  ##
  return(list(pc.dplot = pc.dplot,
              pc.check = pc.check)
         )
}
