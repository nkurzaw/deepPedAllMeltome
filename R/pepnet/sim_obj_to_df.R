#' Convert \code{rcorr} object to \code{data.frame}
#'
#' @param sim_obj \code{rcorr} correlation object
#' @param feature_data \code{data.frame} with additional feature information
#' @param which_features vector of feature columns that should be included
#' 
#' @importFrom magrittr "%>%"
#' @import dplyr
#'
#' @return \code{data.frame}
#' @export
#'
#' @examples
sim_obj_to_df <- function(sim_obj,
                          feature_data,
                          which_features = c("protein_ids")) {
  
  # select upper triangles of matrices
  mark_upper_tri <- function (mat) {
    mat[lower.tri(mat, diag = TRUE)] <- NA
    return(mat)
  }
  
  upper_sim_matrix <- mark_upper_tri(mat = sim_obj$r)
  upper_p_matrix <- mark_upper_tri(mat = sim_obj$P)
  upper_n_matrix <- mark_upper_tri(mat = sim_obj$n)
  
  # prepare feature data
  features <- feature_data %>%
    dplyr::select(peptide, !!which_features)
  
  features_source <- features %>%
    set_names(paste0("source_", colnames(.)))
  
  features_target <- features %>%
    set_names(paste0("target_", colnames(.)))
  
  # assemble data frames
  to_df <- function (mat, type) {
    mat %>%
      as.data.frame() %>%
      rownames_to_column("source_peptide") %>%
      gather(key = "target_peptide", value = !!type, -source_peptide)
  }
  
  df_sim <- to_df(upper_sim_matrix, type = "r")
  df_p <- to_df(upper_p_matrix, type = "p")
  df_n <- to_df(upper_n_matrix, type = "n")
  
  # merge
  df <- Reduce(function (x, y) full_join(x, y, by = c("source_peptide", "target_peptide")), list(df_sim, df_p, df_n))
  
  # add features and filter out NA similarities
  df_features <- df %>%
    left_join(features_source, by = "source_peptide") %>%
    left_join(features_target, by = "target_peptide") %>%
    filter(!is.na(r))
  
  return(df_features)
}