#' Build graphs from similarities
#'
#' @param similarities \code{list} of similarity \code{data.frame}s
#' @param e_set peptides \code{ExpressionSet}
#' @param filter_params \code{list} of filters
#' @param layout_fun \code{igraph} function for graph layouts, used for consistent plotting
#' @param BPPARAM BiocParallel
#' @param verbose flag to enable verbose reporting
#'
#' @return \code{list} of graphs
#' 
#' @importFrom magrittr "%>%"
#' @importFrom Biobase fData
#' @importFrom igraph graph_from_data_frame set_graph_attr
#' 
#' @export
#'
#' @examples
build_graphs <- function (similarities,
                          e_set,
                          filter_params = list(lower_similarity_cutoff = 0,
                                               lower_n_cutoff = -Inf,
                                               upper_q_cutoff = Inf),
                          layout_fun = layout_with_fr,
                          BPPARAM = BiocParallel::SerialParam(),
                          verbose = TRUE) {
  
  if (verbose) cat("=== Build graphs...\n")
  
  iois <- similarities %>%
    names() %>%
    set_names(.)
  
  graphs <- bplapply(X = iois,
                     FUN = function (ioi) {
                       # select the similarities
                       sim_df <- similarities[[ioi]]
                       
                       # transform similarities to weight
                       sim_df$weight <- sim_df$similarity
                       
                       # filter
                       sim_df_filt <- sim_df %>%
                         mutate(use = 1) %>%
                         mutate(use = if_else(similarity > filter_params$lower_similarity_cutoff, use, 0)) %>%
                         mutate(use = if_else(n > filter_params$lower_n_cutoff, use, 0))
                       
                       if ("q" %in% colnames(sim_df_filt)) {
                         sim_df_filt$use_factor = if_else(sim_df_filt$q < filter_params$upper_q_cutoff, sim_df_filt$use, 0)
                       }
                       
                       # assemble the graph edges
                       edge_df <- sim_df_filt %>%
                         mutate(weight = weight * use) %>%
                         filter(weight > 0) %>%
                         filter(source_peptide != target_peptide) %>%
                         dplyr::select(source_peptide, target_peptide, weight, similarity, use)
                       
                       # get vertex properties
                       vertices_df <- e_set %>%
                         fData() %>%
                         dplyr::select(name = peptide, everything()) %>%
                         filter(name %in% c(edge_df$source_peptide, edge_df$target_peptide))
                       
                       # build the graph model
                       g <- graph_from_data_frame(d = edge_df,
                                                  vertices = vertices_df,
                                                  directed = FALSE)
                       
                       # store the goi in the graph
                       g <- set_graph_attr(g, "ioi", ioi)
                       
                       # apply a layout
                       g <- set_graph_attr(g, "layout", layout_fun(g))
                     },
                     BPPARAM = BPPARAM)
  
  is_error <- lapply(X = graphs, FUN = function (g) class(g) == "try-error") %>% unlist()

  if (any(is_error)) {
    cat(">>>", "Graph generation failed for", names(graphs[is_error]), "\n")
    cat(">>>", "These gene symbols will be eliminated!")
    
    graphs <- graphs[!is_error]
  }
  
  if (verbose) cat(">>>", length(graphs), "graphs built\n")
  
  return(graphs)
}