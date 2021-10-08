#' Cluster modularity
#'
#' @param g graph
#'
#' @return cluster modularity vector
#' 
#' @importFrom magrittr "%>%"
#' @importFrom igraph V E "%--%"
#' 
#' @export
#'
#' @examples
cluster_modularity <- function (g) {
  res <- vapply(X = V(g)$membership %>% unique() %>% as.character(),
                FUN = function (community) {
                  # verteces within the community
                  v_c <- V(g)[membership == community]
                  
                  # verteces outside the community
                  v_o <- V(g)[membership != community]
                  
                  # edges within the community
                  e_c <- E(g)[v_c %--% v_c][weight > 0]
                  
                  # edges from the community to outside of the community
                  e_o <- E(g)[v_c %--% v_o][weight > 0]
                  
                  return((sum(e_c$weight) / length(E(g))) - ((sum(e_c$weight) + sum(e_o$weight))^2 / (length(E(g)))^2))
                },
                FUN.VALUE = numeric(1))
  
  return(res)
}