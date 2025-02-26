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
#' @return \item{lp.check}{0 = OS can be conducted; 2 = OS cannot be conducted}
#' 
#' @export exmLP.2ipd
#'
#' @author Lillian Yau
#' 
## exmLP.2ipd(ipd1, ipd2) ## this would be the example, but ipd1 and ipd2 are not in the package yet

exmLP.2ipd <- function (ipd1, ipd2, vars_to_match = NULL, cat_vars_to_01 = NULL, mean.constrained = FALSE)
{
  ## check vars_to_match
  vars_to_match <- .check_data(ipd1, 
                              ipd2, 
                              v.ars_to_match = vars_to_match,
                              c.at_vars_to_01 = cat_vars_to_01
  )
  ## extract only vars_to_match from both ipd's
  ipd1 <- data.frame(ipd1[vars_to_match])
  ipd2 <- data.frame(ipd2[vars_to_match])
  
  if(!is.null(cat_vars_to_01)){

    ## save original input data for later use
    ipd1.o <- ipd1 
    ipd2.o <- ipd2
    
    ## convert categorical variables to indicator variables
    ipd1 <- .cat201_minus1(ipd1.o, v.cat = cat_vars_to_01)
    ipd2 <- .cat201_minus1(ipd2.o, v.cat = cat_vars_to_01)
  }
  ##
  ## derivation for lpCheck starts here ::::
  ##
  ipd <- as.data.frame(rbind(-1 * ipd1, ipd2))
  oneszeros <- c(rep(1, nrow(ipd1)), rep(0, nrow(ipd2)))
  zerosones <- c(rep(0, nrow(ipd1)), rep(1, nrow(ipd2)))
  ipd <- as.data.frame(cbind(ipd, oneszeros, zerosones))
  p <- ncol(ipd1)
  ## f.con is the A matrix's left 3 colns in the appendix of the paper
  f.con <- as.matrix(t(ipd))  
  f.obj <- rep(0.5, ncol(f.con))
  f.rhs <- as.data.frame(t(c(rep(0, p), 1, 1)))
  f.dir <- rep("=", p + 2)
  
  if (mean.constrained == TRUE) { 
    
    ## re-define ipd1o and ipd2o keeping the reference level for all categorical variables
    if(!is.null(cat_vars_to_01)){
      ipd1 <- .cat201(ipd1.o, v.cat = cat_vars_to_01)
      ipd2 <- .cat201(ipd2.o, v.cat = cat_vars_to_01)
    }
    ##
    ipd1.bar <- colMeans(ipd1)
    ipd2.bar <- colMeans(ipd2)
    x <- as.data.frame(rbind(ipd1.bar, ipd2.bar))
    bar.min <- apply(x, 2, min)
    bar.max <- apply(x, 2, max)
    f.rhs <- cbind(f.rhs, 2 * t(bar.min), 2 * t(bar.max))
    f.dir <- c(f.dir, rep(">=", ncol(ipd1)), rep("<=", ncol(ipd1)))
    f.con <- data.frame(rbind(f.con, 
                              cbind(t(ipd1), t(ipd2)), 
                              cbind(t(ipd1), t(ipd2))))
  }
  lp.check <- lpSolve::lp(direction = "max", objective.in = f.obj, 
                          const.mat = f.con, const.dir = f.dir, const.rhs = f.rhs, 
                          transpose.constraints = TRUE)$status
  return(list(lp.check = lp.check))
}


.check_data <- function(ipd1, 
                        ipd2, 
                        v.ars_to_match,
                        c.at_vars_to_01){
  ## extract variable names if input is NULL
  if(is.null(v.ars_to_match)){
    vars_to_match_ipd1 <- colnames(ipd1)
    
    vars_to_match_ipd2 <- colnames(ipd2)
    
    ## check if variables are the same in ipd1 and ipd2
    if(!(setequal(vars_to_match_ipd1, vars_to_match_ipd2)))
      stop("ipd1 and ipd2 do not have the same variables.")
    
    ## otherwise
    v.ars_to_match <- vars_to_match_ipd1
  } else{
    ## if v.ars_to_match is user provided:::
    # Check if 'vars_to_match' are in the both data sets
    if(!all(v.ars_to_match %in% colnames(ipd1))) {
      stop("Some `vars_to_match` are not in ipd1")
    }
    if(!all(v.ars_to_match %in% colnames(ipd2))) {
      stop("Some `vars_to_match` are not in ipd2")
    }
  }
  
  ## check if 'cat_vars_to_01' are in 'vars_to_match'
  if(!all(c.at_vars_to_01 %in% v.ars_to_match)) 
    stop("Some categorical variables are not in `vars_to_match`")
  
  ipd1 <- ipd1[v.ars_to_match]
  ipd2 <- ipd2[v.ars_to_match]
  
  ## :::::::::::::::::: check cat_vars_to_01 ::::::::::::::::::::::
  
  ## T/F for each coln in ipd1 for character variables
  is_char_ipd1 <- sapply(ipd1, is.character)
  
  if(sum(is_char_ipd1) > length(c.at_vars_to_01))
    stop('There are more character type variables in `vars_to_match` in ipd1 then specified in `cat_vars_to_01`.')
  
  is_char_ipd2 <- sapply(ipd2, is.character)
  
  if(sum(is_char_ipd2) > length(c.at_vars_to_01))
    stop('There are more character type variables in `vars_to_match` in ipd2 then specified in `cat_vars_to_01`.')
  
  return(v.ars_to_match)
}



.cat201_minus1 <- function(df, v.cat) {
  ## in cat201_minus1(), k-1 indicator variables are created for 
  ## ... all categorical variables (including binary) with k levels
  
  # Check if selected columns are in the dataframe
  if(!all(v.cat %in% colnames(df))) stop("Some columns not present in the dataframe")
  
  # Convert columns in v.cat to factors if they are not already
  df[v.cat] <- lapply(df[v.cat], factor)
  
  # Get dummy variables for each categorical variable
  dummy_vars <- lapply(v.cat, function(var) {
    model_matrix <- model.matrix(~ as.factor(df[[var]])) ## w/out reference level
    colnames(model_matrix) <- c('(Intercept)', 
                                paste(var, 
                                      ## label new indicators as X1.B, or X3.B ... X3.D
                                      ## ... i.e. without the lowest level
                                      levels(df[[var]])[2:length(levels(df[[var]]))], 
                                      sep = '.'))
    as.data.frame(model_matrix)
  })
  
  # Combine all the dummy variables with the original dataframe
  dummy_df <- do.call(cbind, dummy_vars)
  
  ## remove colns with names "Intercept"
  cols_to_remove <- grepl("(Intercept)", colnames(dummy_df))
  dummy_df <- dummy_df[, !cols_to_remove]
  
  v.dummy <- names(dummy_df)
  
  # Remove original categorical columns, but need to worry when only 1 numerical var remaining
  all.var <- names(df)       ## names of all variables in original dataframe
  x <- match(v.cat, all.var) ## positions of categorical variables in all variable name list
  v.num <- all.var[-x]       ## names of numerical variables
  
  # Combine the original dataframe with dummy variables
  ## keep only original numerical variables and (new) dummy variables
  df <- cbind(df, dummy_df)[c(v.num, v.dummy)]
  
  return(df)
}


.cat201 <- function(df, v.cat) {
  ## in cat201(), 1 indicator variable is created for binary variables
  ## ... k indicator variables are created for categorical variables with k levels
  
  # Check if selected columns are in the dataframe
  if(!all(v.cat %in% colnames(df))) stop("Some columns not present in the dataframe")
  
  df_not_used <- df[!(colnames(df) %in% v.cat)]
  
  # Convert columns in v.cat to factors if they are not already
  ## !!! df now contains only the binary/categorical variables
  df <- data.frame(lapply(df[v.cat], factor))
  
  # Identify (T/F) binary variables (variables with exactly two levels)
  binary_vars_TF <- sapply(df[v.cat], function(col) length(unique(col)) == 2)
  
  # Separate categorical (non-binary) variables
  categorical_df <- data.frame(df[v.cat[!binary_vars_TF]])
  
  ## create k dummy variables for each categorical variable with k levels
  dummy_categorical_vars <- lapply(v.cat[!binary_vars_TF], 
                                   function(var) {
                                     model_matrix <- model.matrix(~ as.factor(categorical_df[[var]]) - 1) ## all levels
                                     colnames(model_matrix) <- paste(var, levels(categorical_df[[var]]), sep = '.')
                                     as.data.frame(model_matrix)
                                   })
  dummy_df <- do.call(cbind, dummy_categorical_vars)
  
  ## ..................................................................... ##
  ## if there are binary variables:
  if(any(binary_vars_TF)){
    
    # Separate binary variables
    binary_df <- data.frame(df[v.cat[binary_vars_TF]])
    
    # create 1 dummy variable for each binary variable
    dummy_binary_vars <- lapply(v.cat[binary_vars_TF], 
                                function(var) {
                                  ##model_matrix <- model.matrix(~ as.factor(df[[var]]) - 1) ## w all levels
                                  model_matrix <- model.matrix(~ as.factor(binary_df[[var]])) ## w/out reference level
                                  colnames(model_matrix) <- c('(Intercept)', 
                                                              paste(var, 
                                                                    ## label new indicators as X1.B, or X3.B ... X3.D
                                                                    ## ... i.e. without the lowest level
                                                                    levels(binary_df[[var]])[2:length(levels(binary_df[[var]]))], 
                                                                    sep = '.'))
                                  as.data.frame(model_matrix)
                                })
    ## combine dummy binary indicators into a dataframe
    dummy_binary_df <- do.call(cbind, dummy_binary_vars)
    
    # Combine with the dummy variables for other categorical variables
    dummy_df <- data.frame(cbind(dummy_binary_df, dummy_df))
  }
  
  # Combine the not re-coded original data with new dummy variables
  df <- cbind(df_not_used, dummy_df)
  
  ## remove colns with names "Intercept"
  cols_to_remove <- grepl("(Intercept)", colnames(df))
  df <- df[, !cols_to_remove]
  
  return(df)
}