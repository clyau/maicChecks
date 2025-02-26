check_data <- function(ipd1, 
                       ipd2, 
                       v.ars_to_match = vars_to_match,
                       c.at_vars_to_01 = cat_vars_to_01
){
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
