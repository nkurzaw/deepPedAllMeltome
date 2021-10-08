#' Sample loading normalisation
#'
#' @param e_set \code{ExpressionSet}
#'
#' @return \code{ExpressionSet}
#' 
#' @importFrom magrittr "%>%"
#' @importFrom Biobase pData exprs
#' @importFrom matrixStats colSums2
#' 
#' @export
#'
#' @examples
normalise_sample_loading <- function (e_set) {
  
 # which sets are there?
  sets <- e_set %>%
    pData() %>%
    .$set %>%
    unique()
  
  col_sums <- sapply(X = sets,
                     FUN = function (current_set) {
                       e_set %>%
                         .[, pData(e_set)$set == current_set] %>%
                         exprs() %>%
                         colSums2(na.rm = TRUE)
                     },
                     simplify = FALSE)
  
  # determine the overall mean
  overall_mean <- mean(col_sums %>% unlist())
  
  # calculate normalisation factors
  norm_factors <- lapply(col_sums, function (x) overall_mean / x)
  
  # apply norm factors
  data_norm_list <- lapply(X = sets,
                           FUN = function (current_set) {
                               e_set %>%
                               .[, pData(e_set)$set == current_set] %>%
                               exprs() %>%
                               sweep(MARGIN = 2,
                                     STATS = norm_factors[[current_set]],
                                     FUN = "*")
                           })
  
  data_norm <- do.call(cbind, data_norm_list)[, colnames(e_set %>% exprs())]
  
  # check if peptide and position order stayed the same
  stopifnot(identical(fData(e_set)$peptide, row.names(data_norm)))
  stopifnot(identical(row.names(pData(e_set)), colnames(data_norm)))
  
  # write back to ExpressionSet
  exprs(e_set) <- data_norm
  
  return(e_set)
}