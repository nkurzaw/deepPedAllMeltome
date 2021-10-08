#' Plot a graph
#'
#' @param g \code{graph}
#' @param color_attr attribute to be used for coloring
#' @param colors color pallette (named vector)
#' @param scale_edges scale factor for edge weight visualization
#' @param ... additional parameters oto the plot function
#'
#' @return plot
#' 
#' @importFrom igraph get.graph.attribute get.vertex.attribute V E
#' @importFrom magrittr "%>%"
#' @importFrom RColorBrewer brewer.pal
#' 
#' @export
#'
#' @examples
plot_graph <- function (g,
                        color_attr = "membership",
                        colors = NULL,
                        scale_edges = 5,
                        ...) {
  
  ioi <- get.graph.attribute(g, "ioi")
  
  col_attr <- get.vertex.attribute(graph = g,
                                   name = color_attr,
                                   index = V(g)$name)
  
  if (!is.null(colors)) {
    color_vec <- col_attr %>%
      as.character() %>%
      colors[.]
  } else {
    color_vec <- brewer.pal(n = min(12, max(3, length(unique(col_attr)))), name = "Set3") %>%
      .[as.numeric(factor(col_attr))] %>%
      set_names(col_attr)
  }
  
  plot(g,
       vertex.color = color_vec,
       vertex.label = NA,
       edge.width = E(g)$weight * scale_edges,
       main = paste0(ioi, ", colored by ", color_attr),
       ...)
}