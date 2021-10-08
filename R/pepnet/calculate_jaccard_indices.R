calculate_jaccard_indices <- function (confusion_matrices,
                                       BPPARAM = BiocParallel::SerialParam(),
                                       verbose = TRUE) {
  
  if (verbose) cat("=== Calculating Jaccard indices...\n")
  
  indices <- bplapply(X = names(confusion_matrices),
                      FUN = function (id) {
                        conf_mat <- as.matrix(confusion_matrices[[id]])
                        
                        # calculate the group unions
                        union_mat <- matrix(data = NA, nrow = nrow(conf_mat), ncol = ncol(conf_mat))
                        for (i in seq_len(nrow(conf_mat))) {
                          for (j in seq_len(ncol(conf_mat))) {
                            union_mat[i, j] <- sum(conf_mat[i, ], na.rm = TRUE) + sum(conf_mat[, j], na.rm = TRUE) - conf_mat[i, j]
                          }
                        }
                        
                        # calculate the Jaccard indices
                        jaccard_mat <- conf_mat / union_mat
                        
                        # aggregate
                        df <- data.frame(id = id,
                                         max_jaccard = max(jaccard_mat, na.rm = TRUE),
                                         mean_max_jaccard = rowMaxs(jaccard_mat, na.rm = TRUE) %>% mean(),
                                         num_communities_row = nrow(jaccard_mat),
                                         num_communities_col = ncol(jaccard_mat),
                                         stringsAsFactors = FALSE)
                        
                        return(df)
                      },
                      BPPARAM = BPPARAM) %>%
    do.call(rbind, .)
  
  return(indices)
}