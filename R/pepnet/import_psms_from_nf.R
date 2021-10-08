#' Import a PSM table from a typical ddamsproteomics Nextflow pipeline 
#'
#' @param file path to the input file; typically target_psmtable.txt
#' @param id_col column name with the protein or gene symbol level identifiers
#' @param peptide_col column name with the peptide identifiers
#' @param protein_id_col column name for protein variant annotation
#' @param quan_regex regular expression that identifies columns holding the quantitative information
#'
#' @return \code{ExpressionSet} containing PSMs as rows and samples as columns
#' 
#' @importFrom data.table fread
#' @importFrom magrittr "%>%"
#' @importFrom Biobase AnnotatedDataFrame ExpressionSet
#' @import dplyr
#' 
#' @export
#'
#' @examples
#' psms <- import_psms_from_nf(file = file.path("path", "to", "target_psmtable.txt"),
#'                             id_col = "Gene Symbol)
import_psms_from_nf <- function (file,
                                 id_col = "Gene Symbol",
                                 peptide_col = "Peptide",
                                 protein_id_col = "Protein",
                                 quan_regex = "tmt10plex") {
  
  # read the psm file
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
                  scan_num = ScanNum,
                  charge = Charge,
                  psm_q_value = `PSM q-value`,
                  peptide_q_value = `peptide q-value`,
                  set = `Biological set`,
                  fraction = Fraction,
                  ms1_area = `MS1 area`,
                  matches(quan_regex)) %>%
    mutate(protein_ids = gsub("\\(pre=[A-Z\\-],post=[A-Z\\-]\\)", "", protein_ids)) %>%
    mutate(id = gsub(";;+", ";", id)) %>%
    mutate(id = gsub("^;", "", id)) %>%
    mutate(id = gsub(";$", "", id))  %>%
    separate(id, into = c("first_id"), sep = ";", remove = FALSE, extra = "drop") %>%
    separate(protein_ids,
             into = "first_protein_id",
             sep = ";",
             extra = "drop",
             remove = FALSE)
  
  # split off quantitative data
  quan_data <- data %>%
    dplyr::select(matches(quan_regex)) %>%
    as.matrix()
  
  # split off feature data
  feature_data <- data %>%
    dplyr::select(-matches(quan_regex))
  
  # build ExpressionSet
  e_set <- ExpressionSet(assayData = quan_data,
                         featureData = AnnotatedDataFrame(feature_data))
  
  return(e_set)
}