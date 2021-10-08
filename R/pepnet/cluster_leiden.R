#' Wrapper for the Leiden community detection algorithm
#'
#' @param graph graph
#'
#' @return communities
#' 
#' @importFrom leiden leiden
#' @importFrom igraph E V vcount
#' 
#' @export
#'
#' @examples
cluster_leiden <- function (graph) {
  
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  
  # run the algorithm
  memberships <- leiden(graph,
                        weights = E(graph)$weight,
                        n_iterations = -1)
  
  # name the memberships
  names(memberships) <- V(graph)$name
  
  # assemble the result object
  res <- list(algorithm = "leiden",
              vcount = vcount(graph),
              membership = memberships,
              merges = NULL)
  
  # set the class
  class(res) <- "communities"
  
  return(res)
}