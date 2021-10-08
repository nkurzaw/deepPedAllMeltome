#' Calculate ratios to lowest temperature for TPP experiments
#'
#' @param e_set \code{ExpressionSet} with quantitative feature information
#' @param sample_col column name of the samples identifier
#'
#' @return \code{ExpressionSet} with the ratios
#' 
#' @importFrom magrittr "%>%"
#' @importFrom Biobase pData exprs
#' 
#' @export
#'
#' @examples
build_ratios_to_lowest_temperature <- function (e_set,
                                                sample_col = "sample_name") {
  
  # separate by sample
  sample_names <- e_set %>%
    pData() %>%
    .[[sample_col]] %>%
    unique()
  
  ratio_matrices <- sapply(X = sample_names,
                           FUN = function (current_sample_name) {
                             # subset
                             e_set_sample <- e_set[, pData(e_set)[[sample_col]] == current_sample_name]
                             
                             # lowest temperature column
                             denom_col <- e_set_sample %>%
                               pData() %>%
                               filter(temperature == min(temperature)) %>%
                               .$sample_id
                             
                             ratio_matrix <- e_set_sample %>%
                               exprs() %>%
                               sweep(MARGIN = 1,
                                     STATS = .[, denom_col],
                                     FUN = "/")
                             
                             return(ratio_matrix)
                           },
                           simplify = FALSE,
                           USE.NAMES = TRUE)
  
  combined_ratio_matrices <- do.call(cbind, ratio_matrices)
  
  # replace non-finite values by NA
  combined_ratio_matrices[!is.finite(combined_ratio_matrices)] <- NA
  
  # change row order to original ExpressionSet
  combined_ratio_matrices <- combined_ratio_matrices[, colnames(exprs(e_set))]
  
  # check, if the colnames are still in original order
  stopifnot(identical(colnames(exprs(e_set)), colnames(combined_ratio_matrices)))
  
  # assign matrix to ExpressionSet
  ratios_e_set <- e_set
  exprs(ratios_e_set) <- combined_ratio_matrices
  
  return(ratios_e_set)
}