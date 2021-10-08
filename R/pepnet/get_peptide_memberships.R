#' Retrieve peptide community memberships
#'
#' @param graphs \code{list} of graphs
#'
#' @return \code{data.frame}
#' 
#' @importFrom magrittr "%>%"
#' @importFrom igraph V
#' 
#' @export
#'
#' @examples
get_peptide_memberships <- function (graphs) {
  lapply(X = names(graphs),
         FUN = function (ioi) {
           g <- graphs[[ioi]]
           if (length(V(g)) == 0) return(NA)
           data.frame(ioi = ioi,
                      peptide = V(g)$name,
                      membership = V(g)$membership,
                      stringsAsFactors = FALSE)
         }) %>%
    .[!is.na(.)] %>%
    do.call(rbind, .)
}