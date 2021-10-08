two_way_anova <- function (graphs,
                           e_set,
                           include_ambiguous_ids = TRUE,
                           group1,
                           group2,
                           target = "value",
                           BPPARAM = BiocParallel::SerialParam(),
                           verbose = TRUE) {
  
  dfs <- bplapply(X = graphs,
                  FUN = function (g) {
                    ioi <- get.graph.attribute(g, "ioi")
                    
                    # select the quantitative data
                    e_set_ioi <- get_ioi_e_set(ioi = ioi,
                                               e_set = e_set,
                                               include_ambiguous_ids = include_ambiguous_ids)
                    
                    if (nrow(e_set_ioi) < 3) return(NA)
                    
                    df <- tidy_e_set(e_set = e_set_ioi) %>%
                      dplyr::rename(name = peptide)
                    
                    # merge graph data
                    vertex_df <- data.frame(name = get.vertex.attribute(g, "name"),
                                            membership = get.vertex.attribute(g, "membership"),
                                            is_unique = if_else(grepl(";", get.vertex.attribute(g, "id")), "ambiguous", "unique"),
                                            stringsAsFactors = FALSE)
                    
                    df_g <- df %>%
                      left_join(vertex_df, by = "name") %>%
                      mutate(membership = if_else(is.na(membership), "NA", as.character(membership)))
                    
                    # mutate the parameter columns to factors
                    df_aov <- data.frame(group1 = factor(df_g[[group1]]),
                                         group2 = factor(df_g[[group2]]),
                                         target = df_g[[target]]) %>%
                      filter(!is.na(target))
                    
                    if (length(unique(df_aov$group1)) < 2 || length(unique(df_aov$group2)) < 2) return(NA)
                    
                    anova <- aov(target ~ group1 + group2, data = df_aov)
        
                    unbalanced <- Anova(anova)
                    
                    unbalanced_df <- data.frame(ioi = ioi,
                                                group1_p = unbalanced["group1", "Pr(>F)"],
                                                group2_p = unbalanced["group2", "Pr(>F)"])
                    
                    return(unbalanced_df)
                  },
                  BPPARAM = BPPARAM) %>%
    .[!is.na(.)]
  
  df <- do.call(rbind, dfs)
  
  return(df)
}