#' Calculate correlations
#'
#' @param data \code{matrix} (columns: peptides, rows: samples)
#' @param method correlation method
#' @param feature_data \code{data.frame} with additional feature information
#' @param verbose flag to enable verbose reporting
#'
#' @return \code{data.frame} with correlations and additional parameters
#' 
#' @importFrom Hmisc rcorr
#' @importFrom fdrtool fdrtool
#' 
#' @export
#'
#' @examples
calculate_correlations <- function (data,
                                    method = c("pearson", "spearman"),
                                    feature_data,
                                    verbose = TRUE) {
  
  # calculate correlations
  sim_obj <- rcorr(data, type = method)
  
  sim_df <- sim_obj_to_df(sim_obj = sim_obj, feature_data = feature_data)
  
  if (nrow(sim_df) > 0) {
    # run fdrtool
    fdrtool_obj <- fdrtool(x = sim_df$p[!is.na(sim_df$p)],
                           statistic = "pvalue",
                           verbose = verbose,
                           plot = FALSE)
    
    # add info to dataframe
    sim_df$q <- NA
    sim_df[!is.na(sim_df$p), "q"] <- fdrtool_obj$qval
  } else {
    sim_df <- NULL
  }
  
  return(sim_df)
}