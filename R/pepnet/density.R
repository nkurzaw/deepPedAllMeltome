#' Density
#'
#' @param g graph 
#'
#' @return density vector
#' 
#' @importFrom magrittr "%>%"
#' @importFrom igraph V E "%--%"
#' 
#' @export
#'
#' @examples
density <- function (g) {
  # for each community
  res <- vapply(X = V(g)$membership %>% unique() %>% as.character(),
                FUN = function (community) {
                  # verteces within the community
                  v_c <- V(g)[membership == community]
                  
                  # edges within the community
                  e_c <- E(g)[v_c %--% v_c][weight > 0]

                  return(sum(e_c$weight) / ((length(v_c) * (length(v_c) - 1)) / 2))
                },
                FUN.VALUE = numeric(1))
  
  return(res)
}