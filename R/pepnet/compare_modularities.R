#' Compare different modularities
#'
#' @param graphs graphs
#' @param x_col vertex attribute for the first modularity
#' @param y_col vertex attribute for the second modularity
#' @param x_use_weights flag, if the first modularity calculation should use edge weights
#' @param y_use_weights flag, if the second modularity calculation should use edge weights
#' @param label_fraction fraction of points labeled in the plot
#' @param h_line_intercept horizontal line intercept
#' @param v_line_intercept vertical line intercept
#' @param BPPARAM BiocParallel
#' @param verbose flag to enable verbose reporting
#'
#' @return \code{data.frame}
#' 
#' @importFrom igraph get.vertex.attribute E get.graph.attribute modularity count_components
#' @import ggplot2
#' 
#' @export
#'
#' @examples
compare_modularities <- function (graphs,
                                  x_col = "membership",
                                  y_col = "first_protein_id",
                                  x_use_weights = TRUE,
                                  y_use_weights = FALSE,
                                  label_fraction = 0.01,
                                  h_line_intercept = NULL,
                                  v_line_intercept = NULL,
                                  BPPARAM = BiocParallel::SerialParam(),
                                  verbose = TRUE) {
  
  if (verbose) cat("=== Compare modularities...\n")
  
  extract_numeric_attribute <- function (g, attr) {
    g %>%
      get.vertex.attribute(name = attr) %>%
      factor() %>%
      as.numeric()
  }
  
  # collect the modularities
  modularities <- bplapply(X = graphs,
                           FUN = function (g) {
                             # the x_col memberships
                             x_membership <- extract_numeric_attribute(g = g, attr = x_col)
                             y_membership <- extract_numeric_attribute(g = g, attr = y_col)
                             
                             # set the weights
                             x_weights <- NULL
                             if (x_use_weights) x_weights <- E(g)$weight
                             
                             y_weights <- NULL
                             if (y_use_weights) y_weights <- E(g)$weight
                             
                             data.frame(gene_symbols = get.graph.attribute(g, name = "ioi"),
                                        x_modularity = modularity(g,
                                                                  membership = x_membership,
                                                                  weights = x_weights),
                                        y_modularity = modularity(g,
                                                                  membership = y_membership,
                                                                  weights = y_weights),
                                        num_components = count_components(graph = g),
                                        stringsAsFactors = FALSE)
                           },
                           BPPARAM = BPPARAM) %>%
    do.call(rbind, .)
  
  # visualise
  p <- modularities %>%
    ggplot(aes(x = x_modularity, y = y_modularity, label = gene_symbols, color = as.factor(num_components))) +
    # geom_smooth(method = "lm", formula = "y ~ x", alpha = .3) +
    geom_point(stroke = 0, size = 2, alpha = .6) +
    stat_dens2d_filter(geom = "text_repel", keep.fraction = label_fraction, color = "gray20") +
    xlab(paste0("Modularity by ", x_col)) + 
    ylab(paste0("Modularity by ", y_col)) +
    ggtitle("Modularity comparison") +
    scale_color_discrete(name = "Num. components") +
    theme(legend.position = "bottom")
  
  if (is.numeric(h_line_intercept)) {
    p <- p + geom_hline(yintercept = h_line_intercept)
  }
  
  if (is.numeric(v_line_intercept)) {
    p <- p + geom_vline(xintercept = v_line_intercept)
  }
  
  print(p)
  
  return(modularities)
}