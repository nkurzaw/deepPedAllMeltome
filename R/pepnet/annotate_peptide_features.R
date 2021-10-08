#' Annotate peptide features with their proteoform membership
#'
#' @param e_set \code{ExpressionSet} with the quantitative data
#' @param graphs \code{list} of graphs with detected communities
#'
#' @return \code{data.frame}
#' 
#' @importFrom magrittr "%>%"
#' @importFrom Biobase fData
#' @import dplyr
#' 
#' @export
#'
#' @examples
annotate_peptide_features <- function (e_set,
                                       graphs) {
  
  # make dataframe from ExpressionSet
  data <- e_set %>%
    fData()
  
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
  
  return(data_memberships)
}