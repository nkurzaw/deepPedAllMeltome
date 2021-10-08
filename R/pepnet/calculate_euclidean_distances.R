#' Calculate weighted Euclidean distances
#'
#' @param data \code{data.frame} with quantitative values grouped in columns
#' @param feature_data \code{data.frame} holding feature meta data
#' @param which_features \code{list} of column names in the \code{feature_data} that should be transferred to the output
#'
#' @return long \code{data.frame} with pairwise Euclidean distances
#' 
#' @importFrom magrittr "%>%"
#' @import dplyr
#' 
#' @export
#'
#' @examples
calculate_euclidean_distances <- function (data,
                                           feature_data,
                                           which_features = c("protein_ids")) {
  
  # determine the total number of samples
  num_samples <- nrow(data)
  
  # get the column names
  col_names <- colnames(data)
  
  df <- c()
  already_calc <- c()
  
  for (col_name_a in col_names) {
    for(col_name_b in col_names) {
      if (!col_name_b %in% already_calc) {
        # get the colums
        col_a <- data[, col_name_a]
        col_b <- data[, col_name_b]
        
        # get number of pair-wise valid values
        num_valid_values <- sum(!is.na(col_a) & !is.na(col_b))
        
        # calculate weighted distance
        distance <- sqrt(sum((col_a - col_b)^2, na.rm = TRUE) * num_samples / num_valid_values)
        
        df <- rbind(df, data.frame(source_peptide = col_name_a,
                                   target_peptide = col_name_b,
                                   d = distance,
                                   n = num_valid_values,
                                   stringsAsFactors = FALSE))
      }
    }
    
    already_calc <- c(already_calc, col_name_a)
  }
  
  # prepare feature data
  features <- feature_data %>%
    dplyr::select(peptide, !!which_features)
  
  features_source <- features %>%
    set_names(paste0("source_", colnames(.)))
  
  features_target <- features %>%
    set_names(paste0("target_", colnames(.)))
  
  # add features and filter out NA similarities
  df_features <- df %>%
    left_join(features_source, by = "source_peptide") %>%
    left_join(features_target, by = "target_peptide") %>%
    filter(!is.na(d))
  
  return(df_features)
}