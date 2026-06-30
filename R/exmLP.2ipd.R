#' @title Checks whether two IPD datasets can be matched with lpSolve::lp
#'
#' @param ipd1 a dataframe with n1 row and p column, where n1 is number of subjects of the first IPD, and p is the number of variables used in standardization.
#' @param ipd2 a dataframe with n2 row and p column, where n2 is number of subjects of the second IPD, and p is the number of variables used in standardization.
#' @param vars_to_match variables used for matching. if NULL, use all variables.
#' @param cat_vars_to_01 variable names for the categorical variables that need to be converted to indicator variables.
#' @param mean.constrained whether to restrict the weighted means to be within the ranges of observed means. Default is FALSE. When it is TRUE, there is a higher chance of not having a solution.
#'
#' @details If dummy variables are already created for the categorical variables in the data set, and are present in \code{ipd1} and \code{ipd2}, then \code{cat_vars_to_01} should be left as NULL.
#'
#' @return \item{lp.check}{0 = exact matching is feasible; 2 = otherwise (no solution; see \code{\link[lpSolve]{lp}} status codes).}
#'
#' @export exmLP.2ipd
#'
#' @author Lillian Yau
#'
#' @references Glimm E and Yau L. (2026). 'Exact matching as an alternative to propensity score matching.' \emph{Statistics in Biopharmaceutical Research}, 18(1):106-116. \doi{10.1080/19466315.2025.2507378}.
#'
#' @examples
#' \dontrun{
#' ipd1 <- sim110[sim110$study == 'IPD A', ]
#' ipd2 <- sim110[sim110$study == 'IPD B', ]
#' x <- exmLP.2ipd(ipd1, ipd2,
#'                 vars_to_match    = paste0('X', 1:5),
#'                 cat_vars_to_01   = paste0('X', 1:3),
#'                 mean.constrained = FALSE)
#' }
exmLP.2ipd <- function(ipd1, ipd2,
                       vars_to_match    = NULL,
                       cat_vars_to_01   = NULL,
                       mean.constrained = FALSE) {

  ## ---- validate inputs ----------------------------------------------------
  vars_to_match <- .check_data(ipd1, ipd2,
                               v.ars_to_match  = vars_to_match,
                               c.at_vars_to_01 = cat_vars_to_01)

  ipd1 <- data.frame(ipd1[vars_to_match])
  ipd2 <- data.frame(ipd2[vars_to_match])

  if (!is.null(cat_vars_to_01)) {
    ipd1.o <- ipd1    ## keep originals (with reference levels) for the
    ipd2.o <- ipd2    ## mean.constrained branch below
    ipd1 <- .cat201_minus1(ipd1.o, v.cat = cat_vars_to_01)
    ipd2 <- .cat201_minus1(ipd2.o, v.cat = cat_vars_to_01)
  }

  ## ---- build LP -----------------------------------------------------------
  ## stack -ipd1 on top of ipd2; append two indicator columns so the LP
  ## fixes the sums of the two weight blocks each to 1
  ipd <- as.data.frame(rbind(-1 * ipd1, ipd2))
  ipd <- as.data.frame(cbind(ipd,
                             oneszeros = c(rep(1, nrow(ipd1)), rep(0, nrow(ipd2))),
                             zerosones = c(rep(0, nrow(ipd1)), rep(1, nrow(ipd2)))))
  p     <- ncol(ipd1)
  f.con <- as.matrix(t(ipd))
  f.obj <- rep(0.5, ncol(f.con))
  f.rhs <- as.data.frame(t(c(rep(0, p), 1, 1)))
  f.dir <- rep('=', p + 2)

  if (isTRUE(mean.constrained)) {
    ## use FULL k-dummy expansion for non-binary categorical vars
    if (!is.null(cat_vars_to_01)) {
      ipd1 <- .cat201(ipd1.o, v.cat = cat_vars_to_01)
      ipd2 <- .cat201(ipd2.o, v.cat = cat_vars_to_01)
    }
    ipd1.bar <- colMeans(ipd1)
    ipd2.bar <- colMeans(ipd2)
    x        <- as.data.frame(rbind(ipd1.bar, ipd2.bar))
    bar.min  <- apply(x, 2, min)
    bar.max  <- apply(x, 2, max)
    f.rhs    <- cbind(f.rhs, 2 * t(bar.min), 2 * t(bar.max))
    f.dir    <- c(f.dir, rep('>=', ncol(ipd1)), rep('<=', ncol(ipd1)))
    f.con    <- data.frame(rbind(f.con,
                                 cbind(t(ipd1), t(ipd2)),
                                 cbind(t(ipd1), t(ipd2))))
  }

  lp.check <- lpSolve::lp(direction             = 'max',
                          objective.in          = f.obj,
                          const.mat             = f.con,
                          const.dir             = f.dir,
                          const.rhs             = f.rhs,
                          transpose.constraints = TRUE)$status

  return(list(lp.check = lp.check))
}


## =============================================================================
##  internal helpers (dot-prefixed, not exported)
## =============================================================================

#' @keywords internal
#' @noRd
.check_data <- function(ipd1, ipd2, v.ars_to_match, c.at_vars_to_01) {

  ## resolve vars_to_match
  if (is.null(v.ars_to_match)) {
    if (!setequal(colnames(ipd1), colnames(ipd2)))
      stop('`ipd1` and `ipd2` do not have the same variables.', call. = FALSE)
    v.ars_to_match <- colnames(ipd1)
  } else {
    if (!all(v.ars_to_match %in% colnames(ipd1)))
      stop('Some `vars_to_match` are not in `ipd1`.', call. = FALSE)
    if (!all(v.ars_to_match %in% colnames(ipd2)))
      stop('Some `vars_to_match` are not in `ipd2`.', call. = FALSE)
  }

  ## cat_vars_to_01 must be a subset of vars_to_match
  if (!all(c.at_vars_to_01 %in% v.ars_to_match))
    stop('Some categorical variables are not in `vars_to_match`.', call. = FALSE)

  ## sanity: there shouldn't be more character cols than declared
  d1 <- ipd1[v.ars_to_match]
  d2 <- ipd2[v.ars_to_match]
  n.char1 <- sum(vapply(d1, is.character, logical(1)))
  n.char2 <- sum(vapply(d2, is.character, logical(1)))
  if (n.char1 > length(c.at_vars_to_01))
    stop('There are more character-type variables in `vars_to_match` in `ipd1` ',
         'than specified in `cat_vars_to_01`.', call. = FALSE)
  if (n.char2 > length(c.at_vars_to_01))
    stop('There are more character-type variables in `vars_to_match` in `ipd2` ',
         'than specified in `cat_vars_to_01`.', call. = FALSE)

  v.ars_to_match
}


#' @keywords internal
#' @noRd
.cat201_minus1 <- function(df, v.cat) {
  ## k-1 indicator variables are created for ALL categorical variables
  ## (including binary) with k levels. Used by LP / QP feasibility paths.
  if (!all(v.cat %in% colnames(df)))
    stop('Some columns not present in the dataframe.', call. = FALSE)

  df[v.cat] <- lapply(df[v.cat], as.factor)

  dummy_list <- lapply(v.cat, function(var) {
    f  <- df[[var]]
    mm <- stats::model.matrix(~ f)[, -1, drop = FALSE]   ## drop intercept
    colnames(mm) <- paste(var, levels(f)[-1], sep = '.')
    as.data.frame(mm)
  })
  dummy_df <- do.call(cbind, dummy_list)
  v.dummy  <- names(dummy_df)

  v.num <- setdiff(names(df), v.cat)
  cbind(df, dummy_df)[c(v.num, v.dummy)]
}


#' @keywords internal
#' @noRd
.cat201 <- function(df, v.cat) {
  ## 1 dummy for binary, k dummies for non-binary categorical.
  ## Used by the mean.constrained = TRUE branch.
  if (!all(v.cat %in% colnames(df)))
    stop('Some columns not present in the dataframe.', call. = FALSE)

  df_not_used <- df[!(colnames(df) %in% v.cat)]
  df <- as.data.frame(lapply(df[v.cat], as.factor))
  names(df) <- v.cat

  binary_TF <- vapply(df[v.cat],
                      function(col) length(unique(col)) == 2L,
                      logical(1))

  ## non-binary: k dummies (no intercept)
  non.bin <- v.cat[!binary_TF]
  if (length(non.bin)) {
    dummy_nb <- lapply(non.bin, function(var) {
      f  <- df[[var]]
      mm <- stats::model.matrix(~ f - 1)
      colnames(mm) <- paste(var, levels(f), sep = '.')
      as.data.frame(mm)
    })
    dummy_nb <- do.call(cbind, dummy_nb)
  } else {
    dummy_nb <- NULL
  }

  ## binary: 1 dummy (drop reference)
  bin <- v.cat[binary_TF]
  if (length(bin)) {
    dummy_b <- lapply(bin, function(var) {
      f  <- df[[var]]
      mm <- stats::model.matrix(~ f)[, -1, drop = FALSE]
      colnames(mm) <- paste(var, levels(f)[-1], sep = '.')
      as.data.frame(mm)
    })
    dummy_b <- do.call(cbind, dummy_b)
  } else {
    dummy_b <- NULL
  }

  cbind(df_not_used, dummy_b, dummy_nb)
}
