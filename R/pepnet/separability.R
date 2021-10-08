#' Separability
#'
#' @param g graph
#' @param aggregation_fun aggregation function
#'
#' @return separability vector
#' 
#' @importFrom magrittr "%>%"
#' @importFrom igraph V E "%--%"
#' 
#' @export
#'
#' @examples
separability <- function (g,
                          aggregation_fun = mean) {
  
  # for each community
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
                  
                  return(aggregation_fun(e_c$weight) / aggregation_fun(e_o$weight))
                },
                FUN.VALUE = numeric(1))
  
  return(res)
}