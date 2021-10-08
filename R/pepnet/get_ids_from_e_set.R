#' Collect ids from an \code{ExpressionSet}
#'
#' @param e_set \code{ExpressionSet}
#'
#' @return named \code{character vector} with unique identifiers
#' 
#' @importFrom magrittr "%>%"
#' @importFrom Biobase fData
#' @importFrom stringr str_split
#' 
#' @export
#'
#' @examples
get_ids_from_e_set <- function (e_set) {
  
  gois <- e_set %>%
    fData() %>%
    .$id %>%
    str_split(";") %>%
    unlist() %>%
    .[. != ""] %>%
    .[. != " "] %>%
    .[!is.na(.)] %>%
    unique() %>%
    trimws() %>%
    set_names(.)
  
  return(gois)
}