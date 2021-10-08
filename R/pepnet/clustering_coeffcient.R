#' Clustering coefficient
#'
#' @param g graph
#'
#' @return clustering coefficent vector
#' 
#' @importFrom magrittr "%>%"
#' @importFrom igraph V E induced_subgraph diameter
#' 
#' @export
#'
#' @examples
clustering_coefficient <- function (g) {
  res <- vapply(X = V(g)$membership %>% unique() %>% as.character(),
                FUN = function (community) {
                  # verteces within the community
                  v_c <- V(g)[membership == community]
                  
                  # subgraph
                  g_community <- induced_subgraph(graph = g, vids = as.numeric(v_c))
                  
                  return(transitivity(graph = g_community, type = "global"))
                },
                FUN.VALUE = numeric(1))
  
  return(res)
}