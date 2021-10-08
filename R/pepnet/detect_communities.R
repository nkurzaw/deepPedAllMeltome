#' Detect graph communities
#'
#' @param graphs \code{list} of graphs
#' @param detect_algorithm algorithm for community detection
#' @param BPPARAM BiocParallel
#' @param verbose flag to enable verbose reporting
#' @param ... additional parameters for the community detecion algorithm
#'
#' @return \code{list} of graphs with additional graph and vertex attributes
#' 
#' @importFrom magrittr "%>%"
#' @importFrom igraph membership set_vertex_attr set_graph_attr
#' 
#' @export
#'
#' @examples
detect_communities <- function (graphs,
                                detect_algorithm = cluster_fast_greedy,
                                BPPARAM = BiocParallel::SerialParam(),
                                verbose = TRUE,
                                ...) {
  
  if (verbose) cat("=== Detect communities...\n")
  
  iois <- graphs %>%
    names() %>%
    set_names(.)
  
  graphs <- bplapply(X = iois,
                     FUN = function (ioi) {
                       # select graph
                       g <- graphs[[ioi]]
                       
                       try({
                         # detect the communities
                         communities <- detect_algorithm(graph = g, ...)
                         
                         memberships <- membership(communities)
                         
                         # store community label in graph
                         g <- set_vertex_attr(graph = g,
                                              name = "membership",
                                              index = names(memberships),
                                              value = memberships)
                         
                         # store the community detection object in the graph
                         g <- set_graph_attr(graph = g,
                                             name = "communities",
                                             value = communities)
                         
                         # store the modularity based on memberships in the graph
                         g <- set_graph_attr(graph = g,
                                             name = "proteoform_modularity",
                                             value = modularity(g, membership = V(g)$membership, weights = E(g)$weight))
                         
                         # determine number of peptides in largest and smalles proteoform
                         # (used as a filter criterion later on)
                         g <- set_graph_attr(graph = g,
                                             name = "num_peptides_in_largest_community",
                                             value = memberships %>% table() %>% max())
                         
                         g <- set_graph_attr(graph = g,
                                             name = "num_peptides_in_smallest_community",
                                             value = memberships %>% table() %>% min())
                         
                         return(g)
                       },
                       silent = FALSE)
                     },
                     BPPARAM = BPPARAM)
  
  is_error <- lapply(X = graphs, FUN = function (g) class(g) == "try-error") %>% unlist()
  
  if (any(is_error)) {
    cat(">>>", "Community detection failed for", names(graphs[is_error]), "\n")
    cat(">>>", "These gene symbols will be eliminated!\n")
    
    graphs <- graphs[!is_error]
  }
  
  return(graphs)
}