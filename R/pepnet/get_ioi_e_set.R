#' Subset an \code{ExpressionSet} by id
#'
#' @param ioi id of interest
#' @param e_set \code{ExpressionSet}
#' @param include_ambiguous_ids flag to include ambiguous ids (e.g. ACT1;ACT2)
#'
#' @return \code{ExpressionSet}
#' 
#' @importFrom Biobase fData
#' 
#' @export
#'
#' @examples
get_ioi_e_set <- function (ioi,
                           e_set,
                           include_ambiguous_ids = TRUE) {
  
  if (include_ambiguous_ids) {
    ioi_search_pattern <- paste0("(^|\\;\\s?)", ioi, "(\\;|$)")
    e_set_ioi <- e_set[grepl(ioi_search_pattern, fData(e_set)$id), ]
  } else {
    e_set_ioi <- e_set[fData(e_set)$id == ioi, ]
  }
  
  return(e_set_ioi)
}