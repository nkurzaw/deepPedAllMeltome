get_goi_e_set <- function (goi = goi,
                           e_set = e_set,
                           id_column = "gene_symbols",
                           include_ambiguous_ids = include_ambiguous_ids,
                           id_sep = id_sep) {
  
  if (include_ambiguous_ids) {
    goi_search_pattern <- paste0(paste0("(^|\\", id_sep, "\\s?)"), goi, paste0("(\\", id_sep, "|$)"))
    e_set_goi <- e_set[grepl(goi_search_pattern, fData(e_set)[[id_column]]), ] 
  } else {
    e_set_goi <- e_set[fData(e_set)[[id_column]] == goi, ] 
  }
  
  return(e_set_goi)
}