#' Box plot
#'
#' @param g graph
#' @param e_set \code{ExpressionSet} with quantitative vertex information
#' @param group_by variable for grouping
#' @param colors color pallette
#' @param include_ambiguous_ids flag, if ambiguous IDs should be used
#' @param jitter flag to enable jitter layout
#' @param theme_param additional theme parameter for plotting
#'
#' @return plot
#' 
#' @importFrom igraph get.graph.attribute get.vertex.attribute
#' @importFrom magrittr "%>%"
#' @import dplyr
#' @import ggplot2
#' 
#' @export
#'
#' @examples
plot_box <- function (g,
                      e_set,
                      group_by,
                      colors = NULL,
                      include_ambiguous_ids = TRUE,
                      jitter = FALSE,
                      theme_param = NULL) {
  
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
      mutate(membership = if_else(is.na(membership), "NA", as.character(membership)))
    
    # name the group column
    df_g$group <- df_g[[group_by]] %>% as.factor()
    
    if (jitter) {
      point_draw_fn <- geom_jitter
    } else {
      point_draw_fn <- geom_point
    }
    
    # plot
    p <- df_g %>%
      ggplot(aes(x = group, y = value, group = group, color = membership)) +
      point_draw_fn(stroke = 0, size = 2, alpha = 0.3) +
      geom_boxplot(aes(group = interaction(group, membership), color = membership), outlier.shape = NA, alpha = 0) 
    
    if (!is.null(colors)) {
      p <- p + scale_color_manual(values = colors)
    }
    
    p <- p +
      xlab("Group") +
      ylab("log2(expression)") +
      ggtitle(paste0(ioi, ", peptides by ", group_by, " and membership"))
    
    if (!is.null(theme_param)) {
      p <- p + theme_param
    }
    
    print(p)
  }
}