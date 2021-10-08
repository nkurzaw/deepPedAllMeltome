#' Import a peptide table from a typical ddamsproteomics Nextflow pipeline 
#'
#' @param file path to the input file; typically target_psmtable.txt
#' @param id_col column name for the protein or gene symbol level identifiers
#' @param peptide_col column name for the peptide identifiers
#' @param protein_id_col column name for protein variant annotation
#' @param quan_regex regular expression that identifies columns holding the quantitative information
#' @param exclude_regex regular expression that identifies quantitative columns which should be excluded
#' @param transform_fun_sample_ids transformation function for sample IDs (the sample column headers)
#' @param extract_sample_meta_from_headers flag to extract sample meta data directly from column headers
#' @param sample_meta_file path to a sample meta data file (tsv), samples are linked by sample_id
#'
#' @return \code{ExpressionSet} containing peptides as rows and samples as columns
#' 
#' @importFrom data.table fread
#' @importFrom magrittr "%>%"
#' @importFrom readr read_tsv
#' @importFrom Biobase AnnotatedDataFrame ExpressionSet
#' @import dplyr
#' 
#' @export
#'
#' @examples
import_peptides_from_nf <- function (file,
                                     id_col = "Associated gene ID(s)",
                                     peptide_col = "Peptide sequence",
                                     protein_id_col = "Protein(s)",
                                     quan_regex = "^X__POOL_",
                                     exclude_regex = "^X__POOL_POOL",
                                     transform_fun_sample_ids = function (x) x,
                                     extract_sample_meta_from_headers = FALSE,
                                     sample_meta_file = NULL) {
  
  # read the peptides file
  data_raw <- fread(file = file,
                    sep = "\t",
                    quote = "",
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    data.table = FALSE)
  
  # filter for relevant columns and some renaming
  data <- data_raw %>%
    dplyr::select(id = !!id_col,
                  peptide = !!peptide_col,
                  protein_ids = !!protein_id_col,
                  matches(quan_regex),
                  -matches(exclude_regex)) %>%
    mutate(id = gsub(";;+", ";", id)) %>%
    mutate(id = gsub("^;", "", id)) %>%
    mutate(id = gsub(";$", "", id))  %>%
    separate(protein_ids,
             into = "first_protein_id",
             sep = ";",
             extra = "drop",
             remove = FALSE)
  
  # extract pheno data from column headers
  # this only works because Jorrit put a lot of info there
  sample_ids <- data %>%
    dplyr::select(matches(quan_regex)) %>%
    colnames() %>%
    transform_fun_sample_ids()
  
  # check, if colnames are still unique
  if (any(duplicated(sample_ids))) {
    stop("Column names (sample IDs) are not unique!")
  }
  
  if (extract_sample_meta_from_headers) {
    pheno_data <- data.frame(sample_id = sample_ids,
                             channel = str_extract(sample_ids, "1[23][0-9][NC]?$"),
                             set = str_extract(sample_ids, "Set[1-9][0-9]?"),
                             sample_name = str_extract(sample_ids, paste0("(?<=", quan_regex, ")\\w+(?=_Set)")),
                             stringsAsFactors = FALSE)
  } else {
    pheno_data <- data.frame(sample_id = sample_ids,
                             sample_name = NA)
  }
  
  # if there is a sample meta file provided, merge this info here
  if (!is.null(sample_meta_file)) {
    
    # if the sample name couldn't be retrieved earlier, delete this column
    if (all(is.na(pheno_data$sample_name))) pheno_data$sample_name <- NULL
    
    sample_meta <- read_tsv(file = sample_meta_file, quote = "")
    
    if ("set" %in% colnames(sample_meta) && extract_sample_meta_from_headers) {
      sample_meta <- sample_meta %>%
        dplyr::select(-set)
    }
    
    if ("channel" %in% colnames(sample_meta) && extract_sample_meta_from_headers) {
      sample_meta <- sample_meta %>%
        dplyr::select(-channel)
    }
    
    pheno_data <- pheno_data %>%
      left_join(sample_meta, by = "sample_id")
  }
  
  # add rownames
  pheno_data <- pheno_data %>%
    mutate(rownames = sample_id) %>%
    column_to_rownames("rownames")
  
  # prepare feature data
  feature_data <- data %>%
    dplyr::select(-matches(quan_regex)) %>%
    separate(id, into = c("first_id"), sep = ";", remove = FALSE, extra = "drop") %>%
    mutate(rownames = peptide) %>%
    column_to_rownames("rownames")
  
  # check that there are no peptide duplicates
  stopifnot(!any(duplicated(feature_data$peptide)))
  
  # extract quantitative values
  quan_data <- data %>%
    dplyr::select(matches(quan_regex)) %>%
    mutate(rownames = row.names(feature_data)) %>%
    column_to_rownames("rownames") %>%
    set_names(transform_fun_sample_ids(colnames(.))) %>%
    as.matrix() %>%
    .[, row.names(pheno_data)]
    
  # build the ExpressionSet
  e_set <- ExpressionSet(assayData = quan_data,
                         phenoData = AnnotatedDataFrame(data = pheno_data),
                         featureData = AnnotatedDataFrame(data = feature_data))
  
  return(e_set)
}