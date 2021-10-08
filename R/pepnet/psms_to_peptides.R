#' Aggregate PSMs to peptides and add sample meta data
#'
#' @param psms \code{ExpressionSet} with PSM information
#' @param summarise_fun summarisation function (default: sum)
#' @param sample_meta_file path to the sample meta data file
#'
#' @return \code{ExpressionSet} with peptide level information
#' 
#' @importFrom Biobase fData exprs
#' @import dplyr
#' 
#' @export
#'
#' @examples
psms_to_peptides <- function(psms,
                             summarise_fun = sum,
                             sample_meta_file,
                             sample_id_col = "sample_id") {
  
  # get the quantification colnames
  quan_cols <- colnames(psms)
  
  # extract the peptide features
  peptide_set_features <- psms %>%
    fData() %>%
    dplyr::select(id, peptide, set, protein_ids, peptide_q_value) %>%
    dplyr::distinct()
  
  # test if peptides and sets are unique
  stopifnot(all(!duplicated(paste0(peptide_set_features$peptide, peptide_set_features$set))))
  
  # take the quantitative data and sum values
  psms_quant_summed <- psms %>%
    exprs() %>%
    as.data.frame() %>%
    cbind(fData(psms) %>% dplyr::select(peptide, set, ms1_area)) %>%
    group_by(set, peptide) %>%
    group_by(num_psms = n(), .add = TRUE) %>% 
    summarise(across(everything(), summarise_fun, na.rm = TRUE), .groups = "drop")
  
  # re-add the feature data
  peptides <- psms_quant_summed %>%
    left_join(peptide_set_features, by = c("peptide", "set"))
  
  # there should be no peptide-set duplicates here
  stopifnot(all(!duplicated(paste0(peptides$peptide, peptides$set))))
  
  # separate data by set
  sets <- peptides %>%
    .$set %>%
    unique() %>%
    sort()
  
  peptides_by_set <- vapply(X = sets,
                            FUN = function (current_set) {
                              peptides_set <- peptides %>%
                                filter(set == current_set) %>%
                                dplyr::select(peptide, num_psms, peptide_q_value, ms1_area, !!quan_cols) %>%
                                set_names(paste0(current_set, "_", colnames(.))) %>%
                                set_names(gsub("Set[0-9]+_peptide$", "peptide", colnames(.)))
                              
                              stopifnot(all(!duplicated(peptides_set$peptide)))
                              
                              return(list(peptides_set))
                            },
                            FUN.VALUE = list(rep(NA, 14)),
                            USE.NAMES = TRUE)
  
  # re-merge them by peptide
  peptides_set <- Reduce(function (x, y) full_join(x, y, by = "peptide"), peptides_by_set)
  
  # check, if the number of final peptides makes sense
  stopifnot(identical(psms %>% fData() %>% .$peptide %>% unique() %>% length(), nrow(peptides_set)))
  
  # re-introduce feature data
  peptides_df <- peptide_set_features %>%
    dplyr::select(-set, -peptide_q_value) %>%
    distinct() %>%
    left_join(peptides_set, by = "peptide")
  
  # check, if total row number has not changed
  stopifnot(identical(nrow(peptides_set), nrow(peptides_df)))
  
  # load sample meta file
  pheno_data <- read_tsv(file = sample_meta_file) %>%
    dplyr::select(sample_id = !!sample_id_col, everything()) %>%
    mutate(rownames = sample_id) %>%
    column_to_rownames("rownames")
  
  # split off feature data
  feature_data <- peptides_df %>%
    dplyr::select(-!!pheno_data$sample_id) %>%
    mutate(rownames = peptide) %>%
    column_to_rownames("rownames")
  
  # split off assay data
  assay_data <- peptides_df %>%
    dplyr::select(peptide, !!row.names(pheno_data)) %>%
    column_to_rownames("peptide") %>%
    as.matrix() %>%
    .[, row.names(pheno_data)]
  
  # build ExpressionSet
  e_set <- ExpressionSet(assayData = assay_data,
                         phenoData = AnnotatedDataFrame(pheno_data),
                         featureData = AnnotatedDataFrame(feature_data))
  
  return(e_set)
}
