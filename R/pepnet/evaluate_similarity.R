#' Measure the similarity between peptides per id
#'
#' @param e_set \code{ExpressionSet} with peptide-level quantitative information
#' @param filter_params filter parameters
#' @param method distance metric
#' @param include_ambiguous_ids flag to include ambiguous ids (e.g. ACT1;ACT2)
#' @param transform_fun function to transform distances to similarities
#' @param BPPARAM BiocParallel
#' @param verbose flag to enable verbose reporting
#'
#' @return \code{list} of \code{data.frame}s
#' 
#' @importFrom BiocParallel bplapply
#' @importFrom Biobase exprs
#' 
#' @export
#'
#' @examples
evaluate_similarity <- function (e_set = peptides,
                                 filter_params = list(min_num_peptides_per_ioi = 0,
                                                      max_num_peptides_per_ioi = Inf,
                                                      min_peptides_per_sample = 2,
                                                      min_samples_with_sufficient_peptides = 17),
                                 method = c("pearson", "spearman", "euclidean"),
                                 include_ambiguous_ids = TRUE,
                                 transform_fun = function (x) x,
                                 BPPARAM = BiocParallel::SerialParam(),
                                 verbose = TRUE) {
  
  if (verbose) cat("=== Evaluate similarity...\n")
  if (verbose) cat(">>> Method:", method, "\n")
  
  # determine identifiers of interest
  iois <- get_ids_from_e_set(e_set = e_set)
  
  if (verbose) cat(">>>", length(iois), "IOIs found\n")
  
  # perform similarity analysis on each goi
  sim_dfs <- lapply(X = iois,
                      FUN = function (ioi) {
                        # select the respective goi rows
                        e_set_ioi <- get_ioi_e_set(ioi = ioi,
                                                   e_set = e_set,
                                                   include_ambiguous_ids = include_ambiguous_ids)

                        # filter out gois with insufficient number of peptides
                        if (nrow(e_set_ioi) < filter_params$min_num_peptides_per_ioi) return(NA)
                        if (nrow(e_set_ioi) == 0) return(NA)
                        if (nrow(e_set_ioi) > filter_params$max_num_peptides_per_ioi) return (NA)
                        
                        # calculate and check overlaps
                        values <- e_set_ioi %>% exprs()
                        
                        valid_values_per_column <- apply(X = !is.na(values),
                                                         MARGIN = 2,
                                                         FUN = sum,
                                                         na.rm = TRUE)
                        
                        valid_values_with_min_peptides <- valid_values_per_column >= filter_params$min_peptides_per_sample
                        
                        if (sum(valid_values_with_min_peptides, na.rm = TRUE) < filter_params$min_samples_with_sufficient_peptides) return(NA)
                        
                        # fetch the expression values
                        expr_ioi <- e_set_ioi %>%
                          exprs() %>%
                          as.data.frame() %>%
                          mutate(rownames = fData(e_set_ioi)$peptide) %>%
                          column_to_rownames("rownames") %>%
                          t()
                        
                        # calculate similarities
                        if (method %in% c("pearson", "spearman")) {
                          
                          sim_df <- calculate_correlations(data = expr_ioi,
                                                           method = method,
                                                           feature_data = fData(e_set_ioi),
                                                           verbose = FALSE)
                          
                          if (!is.null(sim_df)) {
                            sim_df$similarity <- transform_fun(sim_df$r)
                          } else {
                            sim_df <- NA
                          }
                          
                        } else if (method == "euclidean") {
                          
                          sim_df <- calculate_euclidean_distances(data = expr_ioi,
                                                                  feature_data = fData(e_set_ioi)) %>%
                            mutate(similarity = transform_fun(d))
                          
                        } else {
                          stop("Similarity evaluation method not known!")
                        }
                        
                        return(sim_df)
                      }) %>%
    .[!is.na(.)]
  
  if (verbose) cat(">>> Similarity evaluated for", length(sim_dfs), "IDs\n")
  
  return(sim_dfs)
}