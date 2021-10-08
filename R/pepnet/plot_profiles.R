plot_profiles <- function (graphs,
                           e_set,
                           id_column = "gene_symbols",
                           id_sep = ";",
                           peptide_column = "peptide",
                           include_ambiguous_ids = TRUE,
                           x_column = "sample_id",
                           facet_by = NULL,
                           facet_ncol = 1,
                           color_attr = "membership",
                           colors = NULL,
                           BPPARAM = BiocParallel::SerialParam(),
                           output_folder = NULL,
                           verbose = TRUE,
                           ...) {
  
  if (verbose) cat("=== Plotting profiles...\n")
  
  # create folder, if it des not exists
  if (!is.null(output_folder) && !dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  
  gois <- graphs %>%
    names() %>%
    set_names(.)
  
  plots <- bplapply(X = gois,
                    FUN = function (goi) {
                      # select graph
                      g <- graphs[[goi]]
                      
                      # select the quantitative data
                      e_set_goi <- get_goi_e_set(goi = goi,
                                                 e_set = e_set,
                                                 include_ambiguous_ids = include_ambiguous_ids,
                                                 id_sep = id_sep)
                      
                      if (nrow(e_set_goi) == 0) return(NA)
                      
                      # tidied version of quantitative data
                      df <- tidy_e_set(e_set = e_set_goi) %>%
                        dplyr::rename(name = !!peptide_column)
                      
                      # merge graph data
                      vertex_df <- data.frame(name = get.vertex.attribute(g, "name"),
                                              membership = get.vertex.attribute(g, "membership"),
                                              is_unique = if_else(grepl(id_sep, get.vertex.attribute(g, id_column)), "ambiguous", "unique"),
                                              stringsAsFactors = FALSE)
                      
                      df_g <- df %>%
                        left_join(vertex_df, by = "name") %>%
                        mutate(membership = if_else(is.na(membership), "NA", as.character(membership)))
                      
                      # name the x axis column
                      df_g$x_column <- df_g[[x_column]]
                      
                      # name the color column
                      df_g$color <- df_g[[color_attr]]
                      
                      if (is.null(colors)) {
                        colors <- RColorBrewer::brewer.pal(n = max(3, length(unique(df_g$color))), name = "Set3") %>%
                          set_names(as.character(unique(df_g$color)))
                        
                        colors[["NA"]] <- "#494e57"
                      }
                      
                      # peptide profile with colored communities
                      p <- df_g %>%
                        ggplot(aes(x = x_column, y = value, group = name, color = color, alpha = is_unique)) +
                        geom_line() +
                        scale_alpha_manual(values = c("unique" = 0.8, "ambiguous" = 0.2), na.translate = TRUE, na.value = 0.1) +
                        scale_color_manual(values = colors[as.character(df_g$color)], na.value = colors[["NA"]]) +
                        xlab("") + ylab("log2(expression)") +
                        ggtitle(goi) +
                        theme(legend.position = "bottom",
                              axis.text.x = element_text(angle = 90, hjust = 1),
                              legend.title = element_blank())
                      
                      if (!is.null(facet_by)) {
                        p <- p +
                          facet_wrap(facets = facet_by, ncol = facet_ncol)
                      }
                      
                      if (!is.null(output_folder)) pdf(file = file.path(output_folder, paste0(goi, ".pdf")))
                      print(p)
                      if (!is.null(output_folder)) dev.off()
                      
                      return(p)
                    },
                    BPPARAM = BPPARAM)
}