#' Plot peptide profiles for a graph
#'
#' @param g \code{graph}
#' @param e_set \{ExpressionSet} with quantitative peptide information
#' @param include_ambiguous_ids flag to include ambiguous ids (e.g. ACT1;ACT2)
#' @param x_col column to be used for the profiles x axis
#' @param facet_by column to be used for plot facetting
#' @param facet_ncol in case of facetting: number of columns
#' @param color_col column to be used for profile coloring
#' @param colors color pallette as named vector
#' @param legend_position legend position
#' @param main plot title
#'
#' @return plot
#' 
#' @importFrom magrittr "%>%"
#' @importFrom igraph get.vertex.attribute
#' @importFrom RColorBrewer brewer.pal
#' @import dplyr
#' @import ggplot2
#' 
#' @export
#'
#' @examples
plot_profile <- function (g,
                          e_set,
                          include_ambiguous_ids = TRUE,
                          x_col = "sample_id",
                          facet_by = NULL,
                          facet_ncol = 1,
                          color_col = "membership",
                          colors = NULL,
                          legend_position = "bottom",
                          main = "") {
  
  ioi <- get.graph.attribute(g, "ioi")
  
  # select the quantitative data
  e_set_ioi <- get_ioi_e_set(ioi = ioi,
                             e_set = e_set,
                             include_ambiguous_ids = include_ambiguous_ids)
  
  if (nrow(e_set_ioi) > 0) {
    # tidied version of quantitative data
    df <- tidy_e_set(e_set = e_set_ioi) %>%
      dplyr::rename(name = peptide)
    
    # merge graph data
    vertex_df <- data.frame(name = get.vertex.attribute(g, "name"),
                            membership = get.vertex.attribute(g, "membership"),
                            is_unique = if_else(grepl(";", get.vertex.attribute(g, "id")), "ambiguous", "unique"),
                            stringsAsFactors = FALSE)
    
    df_g <- df %>%
      left_join(vertex_df, by = "name") %>%
      mutate(membership = if_else(is.na(membership), "NA", as.character(membership))) %>%
      mutate(is_novel = grepl("nove", protein_ids))
    
    # name the x axis column
    df_g$x_col <- df_g[[x_col]]
    
    # name the color column
    df_g$color <- df_g[[color_col]]
    
    if (is.null(colors)) {
      colors <- brewer.pal(n = max(3, length(unique(df_g$color))), name = "Set3") %>%
        .[seq_len(length(unique(df_g$color)))] %>%
        set_names(as.character(unique(df_g$color)))
      
      colors[["NA"]] <- "#494e57"
    }
    
    # peptide profile with colored communities
    p <- df_g %>%
      ggplot(aes(x = x_col, y = value, group = name, color = color, alpha = is_unique)) +
      geom_line(aes(linetype = is_novel)) +
      scale_alpha_manual(values = c("unique" = 0.7, "ambiguous" = 0.2), na.translate = TRUE, na.value = 0.1) +
      scale_color_manual(values = colors[as.character(df_g$color)], na.value = colors[["NA"]]) +
      xlab("") + ylab("log2(expression)") +
      ggtitle(paste(ioi, main, collapse = ", ")) +
      theme(legend.position = legend_position,
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.title = element_blank())
    
    if (!is.null(facet_by)) {
      p <- p +
        facet_wrap(facets = facet_by, ncol = facet_ncol)
    }
    
    print(p)
  }
}