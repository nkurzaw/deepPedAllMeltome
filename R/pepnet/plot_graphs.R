plot_graphs <- function (graphs,
                         color_attr = "membership",
                         colors = NULL,
                         BPPARAM = BiocParallel::SerialParam(),
                         output_folder = NULL,
                         verbose = TRUE,
                         ...) {
  
  if (verbose) cat("=== Plotting graphs...\n")
  
  if (!is.null(output_folder) && !dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  
  gois <- graphs %>%
    names() %>%
    set_names(.)
  
  plots <- bplapply(X = gois,
                    FUN = function (goi) {
                      g <- graphs[[goi]]
                      
                      col_attr <- get.vertex.attribute(graph = g,
                                                       name = color_attr,
                                                       index = V(g)$name)
                      
                      if (!is.null(colors)) {
                        color_vec <- col_attr %>%
                          as.character() %>%
                          colors[.]
                      } else {
                        color_vec <- RColorBrewer::brewer.pal(n = min(12, max(3, length(unique(col_attr)))), name = "Set3") %>%
                          .[as.numeric(factor(col_attr))] %>%
                          set_names(col_attr)
                      }
                      
                      if (!is.null(output_folder)) pdf(file = file.path(output_folder, paste0(goi, ".pdf")))
                      plot(g,
                           vertex.color = color_vec,
                           vertex.label = NA,
                           main = goi)
                      if (!is.null(output_folder)) dev.off()
                    },
                    BPPARAM = BPPARAM)
}