#' IRS normalisation
#'
#' @param e_set \code{ExpressionSet}
#' @param ref_channels columns of the reference channels
#'
#' @return \code{ExpressionSet}
#' 
#' @importFrom magrittr "%>%"
#' @importFrom Biobase pData exprs
#' 
#' @export
#'
#' @examples
normalise_irs <- function (e_set,
                           ref_channels) {
  
  # retrieve the available sets
  sets <- e_set %>%
    pData() %>%
    .$set %>%
    unique()
  
  # calculate IRS factors
  # TODO: change to vapply
  irs_factors <- sapply(X = sets,
                        FUN = function (current_set) {
                          e_set %>%
                            .[, pData(e_set)$position == ref_channels[[current_set]]] %>%
                            exprs() %>%
                            .[, 1]
                        },
                        simplify = TRUE) %>%
    as.data.frame() %>%
    mutate(average = apply(X = .,
                           MARGIN = 1,
                           FUN = function (x) exp(mean(log(x), na.rm = TRUE))))
  
  # compute the scaling factors
  scaling_factors <- apply(X = irs_factors,
                           MARGIN = 2,
                           FUN = function (col) irs_factors$average / col) %>%
    as.data.frame()
  
  # apply factors
  data_norm_list <- sapply(X = sets,
                           FUN = function (current_set) {
                             e_set %>%
                               .[, pData(e_set)$set == current_set] %>%
                               exprs() %>%
                               sweep(MARGIN = 1,
                                     STATS = scaling_factors[[current_set]],
                                     FUN = "*")
                           },
                           simplify = FALSE)
  
  data_norm <- do.call(cbind, data_norm_list)[, colnames(exprs(e_set))]
  
  # check if row and column names fit
  stopifnot(identical(row.names(exprs(e_set)), row.names(data_norm)))
  stopifnot(identical(colnames(exprs(e_set)), colnames(data_norm)))
  
  exprs(e_set) <- data_norm
  
  return(e_set)
}