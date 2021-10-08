#' Aggregate peptides to proteoforms using a peptide group annotation
#'
#' @param e_set \code{ExpressionSet} with the quantitative data
#' @param graphs \code{list} of graphs with detected communities
#' @param aggregation_fun function to aggregate the peptides, 
#' e.g. \code{sum} or \code{median}
#' @param num_psms_regex regular expression describing the column names 
#' with PSM numbers
#' @param ms1_area_regex regular expression describing the column names 
#' with MS1 areas
#' @param categorical_cols_to_aggregate the categorical columns which 
#' should be aggregated along with the quantitative data
#' @param BPPARAM BiocParallel
#' @param verbose flag to enable verbose reporting
#'
#' @return \code{ExpressionSet}
#' 
#' @importFrom magrittr "%>%"
#' @importFrom Biobase exprs fData pData
#' @importFrom matrixStats rowMeans
#' @import dplyr
#' 
#' @export
#'
#' @examples
aggregate_peptides_to_proteoforms <- function (e_set,
                                               graphs,
                                               aggregation_fun = sum,
                                               num_psms_regex = "num_psms$",
                                               ms1_area_regex = "ms1_area$",
                                               categorical_cols_to_aggregate = c("id", "protein_ids"),
                                               BPPARAM = BiocParallel::SerialParam(),
                                               verbose = TRUE) {
  
  if (verbose) cat("=== Aggregate peptides to proteoforms...\n")
  
  # constants
  quant_identifier <- "__quant_"
  
  # make dataframe from ExpressionSet
  data <- e_set %>%
    exprs() %>%
    as.data.frame() %>%
    set_names(paste0(quant_identifier, colnames(.))) %>%
    cbind(fData(e_set))
  
  # check if peptides are unique
  stopifnot(all(!duplicated(data$peptide)))
  
  # add memberships
  memberships <- get_peptide_memberships(graphs = graphs)
  
  data_memberships <- data %>%
    left_join(memberships, by = "peptide") %>%
    # filter out rows where ioi is NA and the ID contains multiple entries
    filter(!(is.na(ioi) & grepl(";", id))) %>%
    # assign ioi, if there is none
    mutate(ioi = if_else(is.na(ioi), id, ioi)) %>%
    # if the id is not in the graph names, then assign membership to 0
    mutate(membership = if_else(ioi %in% names(graphs), membership, 0)) %>%
    # exclude rows with no id and no membership information
    filter(!is.na(id)) %>%
    filter(!is.na(membership)) %>%
    # create proteoform ID
    mutate(proteoform_id = paste0(ioi, "_", membership))
  
  # aggregate by identifier and membership
  data_aggregated <- data_memberships %>%
    dplyr::select(ioi, membership, proteoform_id, matches(quant_identifier), matches(num_psms_regex), matches(ms1_area_regex)) %>%
    group_by(ioi, membership, proteoform_id) %>%
    group_by(num_peptides = n(), .add = TRUE) %>%
    summarise(across(everything(), list(aggregate = function (x) aggregation_fun(x, na.rm = TRUE),
                                        sd = function (x) sd(x, na.rm = TRUE))),
              .groups = "drop") %>%
    mutate(aggregate_sd = dplyr::select(., ends_with("_sd")) %>% as.matrix() %>% rowMeans(na.rm = TRUE))
  
  # summarise categorical variables for each proteoform
  categorical_data <- data_memberships %>%
    dplyr::select(proteoform_id, !!categorical_cols_to_aggregate)
  
  aggregated_categorical_data <- bplapply(X = categorical_data$proteoform_id %>% unique(),
                                          FUN = function (pf_id) {
                                            df <- categorical_data %>%
                                              filter(proteoform_id == pf_id)
                                            
                                            lapply(X = categorical_cols_to_aggregate,
                                                   FUN = function (col) {
                                                     lapply(df[[col]], str_split, pattern = ";") %>%
                                                       unlist() %>%
                                                       unique() %>%
                                                       paste(collapse = ";")
                                                   }) %>%
                                              do.call(cbind, .) %>%
                                              as.data.frame() %>%
                                              set_names(categorical_cols_to_aggregate) %>%
                                              mutate(proteoform_id = pf_id)
                                          },
                                          BPPARAM = BPPARAM) %>%
    do.call(rbind, .)
  
  # find out how many ambiguous peptides there are per proteoform
  num_ambiguous_peptides_per_proteoform <- data_memberships %>%
    dplyr::select(proteoform_id, id) %>%
    group_by(proteoform_id) %>%
    summarise(num_ambiguous_peptides = sum(grepl(";", id)),
              .groups = "drop")
  
  # add the additional variables to the quantitative data
  data_aggregated_features <- data_aggregated %>%
    left_join(aggregated_categorical_data, by = "proteoform_id") %>%
    left_join(num_ambiguous_peptides_per_proteoform, by = "proteoform_id") %>%
    mutate(ambiguous_peptides_only = num_ambiguous_peptides == num_peptides) %>%
    mutate(ambiguity_ratio = num_ambiguous_peptides / num_peptides)
  
  # rebuild an ExpressionSet
  assay_data <- data_aggregated_features %>%
    dplyr::select(matches(quant_identifier) & ends_with("_aggregate")) %>%
    set_names(gsub(quant_identifier, "", colnames(.))) %>%
    set_names(gsub("_aggregate", "", colnames(.))) %>%
    as.matrix() %>%
    `rownames<-`(data_aggregated_features$proteoform_id)
  
  feature_data <- data_aggregated_features %>%
    dplyr::select(-(matches(quant_identifier) & ends_with("_aggregate"))) %>%
    dplyr::select(proteoform_id, ioi, all_of(categorical_cols_to_aggregate), everything()) %>%
    as.data.frame() %>%
    `rownames<-`(.$proteoform_id)
  
  aggregated_e_set <- ExpressionSet(assayData = assay_data,
                                    phenoData = AnnotatedDataFrame(pData(e_set)),
                                    featureData = AnnotatedDataFrame(feature_data))
  
  if (verbose) cat(">>>", nrow(aggregated_e_set), "proteoform entries generated.\n")
  
  return(aggregated_e_set)
}