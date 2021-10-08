#' Calculate community-level graph quality parameters
#' https://hal.archives-ouvertes.fr/hal-01577343/document
#'
#' @param graphs \code{list} of graphs
#' @param BPPARAM BiocParallel
#'
#' @return \code{data.frame} with parameters
#' 
#' 
#' 
#' @export
#'
#' @examples
calculate_community_parameters <- function (graphs,
                                            BPPARAM = BiocParallel::SerialParam()) {
  
  cat("=== Evaluate community parameters\n")
  
  parameters <- bplapply(X = names(graphs) %>% set_names(.),
                         FUN = function (ioi) {
                           g <- graphs[[ioi]]
                           
                           try({
                             parameters <- list(
                               ioi = rep(ioi, V(g)$membership %>% unique() %>% length()) %>% set_names(V(g)$membership %>% unique()),
                               
                               # separability
                               # ratio between internal connections and external connections of nodes inside a community
                               separability_mean = separability(g = g, aggregation_fun = mean),
                               separability_sum = separability(g = g, aggregation_fun = sum),
                               
                               # embeddedness
                               # reflects how much the direct neighbours of a node belong to its community
                               embeddedness_mean = embeddedness(g = g, aggregation_fun = mean),
                               embeddedness_sum = embeddedness(g = g, aggregation_fun = sum),
                               
                               # density
                               # weights inside community vs. total possible weights inside S
                               density = density(g = g),
                               
                               # compactness
                               # good communities should be dense: internal edges vs. community diameter
                               compactness = compactness(g = g),
                               
                               # clustering coeffcient
                               clustering_coefficient = clustering_coefficient(g = g),
                               
                               # cluster modularity
                               # modularity on community level
                               cluster_modularity = cluster_modularity(g = g)
                             )
                             
                             return(parameters)
                           })
                         },
                         BPPARAM = BPPARAM)
  
  parameters_df <- aggregate_parameters(parameters = parameters)
  
  return(parameters = parameters_df)
}