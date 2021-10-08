#' Perform group-wise t-tests for proteoforms
#'
#' @param graphs \code{list} of graphs
#' @param e_set \code{ExpressionSet} with quantitative values
#' @param group_by variable to group by
#' @param filter_params \code{list} with filter parameters
#' @param include_ambiguous_proteoforms flag to include ambiguous proteoforms
#' @param include_ambiguous_ids flag to include ambiguous ids (e.g. ACT1;ACT2)
#' @param BPPARAM BiocParallel
#'
#' @return \code{data.frame} with log fold changes and p-values
#' 
#' @import BiocParallel bplapply
#' @importFrom magrittr "%>%"
#' @importFrom igraph get.graph.attribute get.vertex.attribute
#' 
#' @export
#'
#' @examples
group_ttest <- function (graphs,
                         e_set,
                         group_by,
                         filter_params = list(min_num_values_per_group = 10,
                                              min_num_values_per_proteoform = 3),
                         include_ambiguous_proteoforms = TRUE,
                         include_ambiguous_ids = TRUE,
                         BPPARAM = BiocParallel::SerialParam()) {
  
  tests <- bplapply(X = graphs,
                    FUN = function (g) {
                      # get ioi
                      ioi <- get.graph.attribute(graph = g, name = "ioi")
                      
                      # fetch quantitative data
                      e_set_ioi <- get_ioi_e_set(ioi = ioi,
                                                 e_set = e_set,
                                                 include_ambiguous_ids = include_ambiguous_ids)
                      
                      # tidy
                      df <- tidy_e_set(e_set = e_set_ioi)
                      
                      # merge graph data
                      vertex_df <- data.frame(peptide = get.vertex.attribute(g, "name"),
                                              membership = get.vertex.attribute(g, "membership"),
                                              is_unique = !grepl(";", get.vertex.attribute(g, "id")),
                                              stringsAsFactors = FALSE)
                      
                      df_g <- df %>%
                        left_join(vertex_df, by = "peptide") %>%
                        filter(!is.na(value))
                      
                      if (!include_ambiguous_proteoforms) {
                        df_g <- df_g %>% filter(is_unique == TRUE)
                      }
                      
                      # go through the groups
                      groups <- lapply(X = df_g[[group_by]] %>% unique(),
                                       FUN = function (group) {
                                         df_g_group <- df_g %>% filter(.[[group_by]] == group)
                                         
                                         if (nrow(df_g_group) < filter_params$min_num_values_per_group) return(NA)
                                         
                                         # go through each membership
                                         memberships <- lapply(X = unique(df_g_group$membership) %>% .[!is.na(.)],
                                                               FUN = function (m) {
                                                                 test <- df_g_group %>% filter(membership == m)
                                                                 ctrl <- df_g_group %>% filter(membership != m | is.na(membership))
                                                                 
                                                                 if (nrow(test) < filter_params$min_num_values_per_proteoform) return(NA)
                                                                 if (nrow(ctrl) < filter_params$min_num_values_per_proteoform) return(NA)
                                                                 
                                                                 logFC <- abs(median(test$value) - median(ctrl$value))
                                                                 p <- t.test(test$value, ctrl$value)$p.value
                                                                 
                                                                 return(data.frame(ioi = ioi,
                                                                                   proteoform_id = paste0(ioi, "_", m),
                                                                                   membership = m,
                                                                                   group = group,
                                                                                   logFC = logFC,
                                                                                   p = p))
                                                               }) %>%
                                           .[!is.na(.)] %>%
                                           do.call(rbind, .)
                                         
                                         return(memberships)
                                       }) %>%
                        .[!is.na(.)] %>%
                        do.call(rbind, .)
                      
                      return(groups)
                    },
                    BPPARAM = BPPARAM) %>%
    .[!is.na(.)] %>%
    do.call(rbind, .) %>%
    mutate(p_adj = p.adjust(p = p, method = "BH"))
  
  return(tests)
}