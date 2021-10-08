#' Calculate graph-level graph quality parameters
#' https://hal.archives-ouvertes.fr/hal-01577343/document
#'
#' @param graphs \code{list} of graphs
#' @param BPPARAM BiocParallel
#'
#' @return \code{data.frame} with parameters
#' 
#' @importFrom BiocParallel bplapply
#' @importFrom magrittr "%>%"
#' @importFrom igraph membership E
#' @import dplyr
#' 
#' @export
#'
#' @examples
calculate_graph_parameters <- function (graphs,
                                        BPPARAM = BiocParallel::SerialParam()) {
  
  cat("=== Evaluate graph parameters\n")
  
  parameters <- bplapply(X = names(graphs) %>% set_names(.),
                         FUN = function (ioi) {
                           g <- graphs[[ioi]]
                           
                           if (!is.null(g$communities$modularity)) {
                             max_modularity <- g$communities$modularity %>% max()
                             num_comms_max_modularity <- length(g$communities$modularity) - which.max(g$communities$modularity) + 1
                           } else {
                             max_modularity <- modularity(g,
                                                          membership = membership(g$communities),
                                                          weights = E(g)$weight)
                             num_comms_max_modularity <- NA
                           }
                           
                           try({
                             return(list(
                               ioi = ioi,
                               
                               # max. modularity
                               max_modularity = max_modularity,
                               
                               # number of communities at max. modularity
                               num_comms_max_modularity = num_comms_max_modularity
                             ))
                           })
                         },
                         BPPARAM = BPPARAM)
  
  parameters_df <- aggregate_parameters(parameters = parameters)
  
  return(parameters = parameters_df)
}