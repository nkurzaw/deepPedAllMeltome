calculate_confusion_matrices <- function (graphs,
                                          truth_variable = "first_master_protein_id",
                                          BPPARAM = BiocParallel::SerialParam(),
                                          verbose = TRUE) {
  
  if (verbose) cat("=== Calculating confusion matrices...\n")
  
  confusion_matrices <- bplapply(X = graphs,
                                 FUN = function (g) {
                                   # get the proteoform communities
                                   proteoform_memberships <- get.vertex.attribute(graph = g,
                                                                                  name = "membership",
                                                                                  index = V(g)) %>%
                                     set_names(V(g)$name)
                                   
                                   proteoform_communities <- proteoform_memberships %>% unique()
                                   
                                   # get the ground truth communities
                                   truth_memberships <- get.vertex.attribute(graph = g,
                                                                             name = truth_variable,
                                                                             index = V(g)) %>%
                                     set_names(V(g)$name)
                                   
                                   truth_communities <- truth_memberships %>% unique()
                                   
                                   # search for each peptide in the protein id communities
                                   confusion_list <- lapply(X = proteoform_communities %>% set_names(as.character(.)),
                                                            FUN = function (proteoform_community) {
                                                              truth_peptide_memberships <- truth_memberships[proteoform_memberships[proteoform_memberships == proteoform_community] %>% names()]
                                                              table(truth_peptide_memberships)
                                                            })
                                   
                                   confusion_matrix <- lapply(X = names(confusion_list),
                                                              FUN = function (proteoform_membership) {
                                                                data.frame(confusion_list[[proteoform_membership]],
                                                                           stringsAsFactors = FALSE) %>%
                                                                  set_names(gsub("Freq", proteoform_membership, colnames(.)))
                                                              }) %>%
                                     Reduce(f = function (x, y) full_join(x, y, by = "truth_peptide_memberships"), x = .) %>%
                                     column_to_rownames("truth_peptide_memberships")
                                   
                                   return(confusion_matrix)
                                 },
                                 BPPARAM = BPPARAM)
  
  return(confusion_matrices)
}