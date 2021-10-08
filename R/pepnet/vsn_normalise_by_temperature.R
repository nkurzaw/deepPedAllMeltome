vsn_normalize_by_temperature <- function(e_set,
                                         grouping_col = "temperature",
                                         inv_trans_fun = inverse_glog2){
  
  unique_groups <- unique(pData(e_set)[[grouping_col]])
  e_set_norm <- e_set
  
  for (group in unique_groups){
    exprs(e_set_norm)[ ,pData(e_set_norm)[[grouping_col]] == group] <- vsn::justvsn(exprs(e_set_norm)[ ,pData(e_set_norm)[[grouping_col]] == group])
  }
  
  exprs(e_set_norm) <- inv_trans_fun(exprs(e_set_norm))
  
  return(e_set_norm)
}

inverse_glog2 <- function (x) {
  ids_below_one <- which(x < 1)
  x_trans <- sinh((log(2) * x) + log(2))
  x_trans[ids_below_one] <- 0
  return(x_trans)
}