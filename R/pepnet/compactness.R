#' Compactness
#'
#' @param g graph 
#'
#' @return compactness vector
#' 
#' @importFrom magrittr "%>%"
#' @importFrom igraph V E induced_subgraph diameter "%--%"
#' 
#' @export
#'
#' @examples
compactness <- function (g) {
  res <- vapply(X = V(g)$membership %>% unique() %>% as.character(),
                FUN = function (community) {
                  # verteces within the community
                  v_c <- V(g)[membership == community]
                  
                  # diameter
                  g_community <- induced_subgraph(graph = g, vids = as.numeric(v_c))
                  diam <- diameter(graph = g_community, directed = FALSE)
                  
                  # edges within the community
                  e_c <- E(g)[v_c %--% v_c][weight > 0]
                  
                  return(sum(e_c$weight) / diam)
                },
                FUN.VALUE = numeric(1))
  
  return(res)
}