get_similarities_per_proteoform <- function (graphs,
                                             modularity_cutoffs = c(-Inf, 0, 0.05, 0.1, 0.2),
                                             min_peptides_per_proteoform = 10,
                                             do_plot = TRUE,
                                             BPPARAM = BiocParallel::SerialParam(),
                                             verbose = TRUE) {
  
  if (verbose) cat("=== Retrieve similarities per proteoform")
  
  graph_modularities <- lapply(graphs, function (g) get.graph.attribute(g, "proteoform_modularity")) %>% unlist()
  
  all_weights <- bplapply(X = modularity_cutoffs %>% set_names(as.character(.)),
                          FUN = function (mod_cutoff) {
                            graphs_above_cutoff <- graphs[graph_modularities > mod_cutoff]
                            
                            graph_weights <- lapply(X = graphs_above_cutoff,
                                                    FUN = function (g) {
                                                      # get the communities
                                                      communities <- g %>% V() %>% .$membership %>% unique()
                                                      
                                                      community_weights <- lapply(X = communities,
                                                                                  FUN = function (community) {
                                                                                    # verteces within the community
                                                                                    v_c <- V(g)[membership == community]
                                                                                    
                                                                                    if (length(v_c) < min_peptides_per_proteoform) return(NA)
                                                                                    
                                                                                    # edges within the community
                                                                                    e_c <- E(g)[v_c %--% v_c][weight > 0]
                                                                                    
                                                                                    # the weights
                                                                                    w <- e_c$weight
                                                                                    
                                                                                    return(median(w))
                                                                                  })
                                                      
                                                      return(median(unlist(community_weights), na.rm = TRUE))
                                                    })
                            return(graph_weights)
                          },
                          BPPARAM = BPPARAM)
  if (do_plot) {
    cutoff_weights_df <- lapply(X = names(all_weights),
                                FUN = function (cutoff) {
                                  data.frame(cutoff = cutoff,
                                             similarities = all_weights[[cutoff]] %>% unlist(recursive = TRUE),
                                             stringsAsFactors = FALSE)
                                }) %>%
      do.call(rbind, .)
    
    p <- cutoff_weights_df %>%
      ggplot(aes(x = similarities, group = cutoff, color = cutoff)) +
      geom_density() +
      xlab("Median within proteoform similarities") + ylab("Density") +
      theme(legend.position = "bottom")
    
    print(p)
  }
  
  return(all_weights)
}