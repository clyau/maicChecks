cat201_minus1 <- function(df, v.cat) {
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