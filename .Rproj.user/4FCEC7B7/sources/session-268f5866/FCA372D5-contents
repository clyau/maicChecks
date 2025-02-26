cat201 <- function(df, v.cat) {
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