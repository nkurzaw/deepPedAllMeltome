#' Plot combo plots of graphs
#'
#' @param graphs \code{list} of graphs
#' @param e_set \code{ExpressionSet} with quantitative information for the verteces
#' @param include_ambiguous_ids flag to include ambiguous IDs
#' @param x_col variable specifying the x axis of profile plots
#' @param facet_by variable column
#' @param facet_ncol number of facet columns
#' @param color_attr variable for coloring
#' @param membership_colors color palette for the proteoforms
#' @param box_group_by grouping variable(s) for the boxplots
#' @param width plot width
#' @param output_folder output folder
#' @param verbose flag to enable verbose reporting
#'
#' @return plots
#' 
#' @importFrom magrittr "%>%"
#' @importFrom igraph get.graph.attribute
#' @importFrom ggplot2 theme
#' 
#' @export
#'
#' @examples
plot_combo <- function (graphs,
                        e_set,
                        include_ambiguous_ids = TRUE,
                        x_col = "sample_id",
                        facet_by = NULL,
                        facet_ncol = 1,
                        color_attr = "membership",
                        membership_colors = NULL,
                        box_group_by = "NEJM.group",
                        width = 12,
                        print_modularity = TRUE,
                        output_folder = NULL,
                        verbose = TRUE) {
  
  if (verbose) cat("=== Plotting combos...\n")
  
  # create folder, if it des not exists
  if (!is.null(output_folder) && !dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  
  iois <- graphs %>%
    names() %>%
    set_names(.)
  
  plots <- lapply(X = iois,
                  FUN = function (ioi) {
                    # select graph
                    g <- graphs[[ioi]]
                    
                    # get modularity
                    modularity <- get.graph.attribute(graph = g, name = "proteoform_modularity")
                    
                    if (print_modularity) {
                      file <- paste0(gsub("\\.", "", as.character(sprintf("%.3f", round(modularity, 3)))), "_", ioi, ".pdf")
                    } else {
                      file <- paste0(ioi, ".pdf")
                    }
                    
                    if (!is.null(output_folder)) pdf(file = file.path(output_folder, file),
                                                     width = width)
                    
                    # plot peptide profiles by membership
                    plot_profile(g = g,
                                 e_set = e_set,
                                 include_ambiguous_ids = include_ambiguous_ids,
                                 x_col = x_col,
                                 facet_by = facet_by,
                                 facet_ncol = facet_ncol,
                                 color_col = "membership",
                                 colors = membership_colors,
                                 main = "Peptide profiles, colored by proteoforms")
                    
                    # plot peptide profiles by protein_ids
                    plot_profile(g = g,
                                 e_set = e_set,
                                 include_ambiguous_ids = include_ambiguous_ids,
                                 x_col = x_col,
                                 facet_by = facet_by,
                                 facet_ncol = facet_ncol,
                                 color_col = "protein_ids",
                                 colors = NULL,
                                 legend_position = "none",
                                 main = "Peptide profiles, colored by protein IDs")
                    
                    # box plot by group
                    for (group in box_group_by) {
                      plot_box(g = g,
                               e_set = e_set,
                               group_by = group,
                               colors = membership_colors,
                               include_ambiguous_ids = include_ambiguous_ids,
                               jitter = TRUE,
                               theme_param = theme(legend.position = "bottom",
                                                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)))
                    }
                    
                    # graph with membership coloring
                    plot_graph(g = g,
                               color_attr = "membership",
                               colors = membership_colors)
                    
                    # graph with protein_id coloring
                    plot_graph(g = g,
                               color_attr = "protein_ids",
                               colors = NULL)
                    
                    if (!is.null(output_folder)) dev.off()
                  })
}