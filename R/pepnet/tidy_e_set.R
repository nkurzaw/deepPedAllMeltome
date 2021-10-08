#' Transform \code{ExpressionSet} to tidy format
#'
#' @param e_set \code{ExpressionSet}
#'
#' @return \code{data.frame}
#' 
#' @importFrom magrittr "%>%"
#' @importFrom Biobase exprs fData
#' @import dplyr
#' 
#' @export
#'
#' @examples
tidy_e_set <- function (e_set) {
  
  df <- e_set %>%
    exprs() %>%
    as.data.frame() %>%
    cbind(., fData(e_set)) %>%
    gather(key = "sample_id", value = "value", seq_len(ncol(e_set))) %>%
    left_join(pData(e_set), by = "sample_id")
  
  return(df)
}