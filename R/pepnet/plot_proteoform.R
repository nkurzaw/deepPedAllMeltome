plot_proteoform <- function (ioi,
                             e_set,
                             x_col = "sample_name",
                             facet_col = NULL,
                             color_col = "membership",
                             colors = c("1" = "#FF5376",
                                        "2" = "#72AFD9",
                                        "3" = "#E3D26F",
                                        "4" = "#A288E3",
                                        "5" = "#1B5299",
                                        "6" = "#68D8D6",
                                        "NA" = "#494e57"),
                             add_splines = FALSE,
                             hide_lines = FALSE,
                             x_label = "Sample name",
                             y_label = "log2 proteoform level",
                             y_limits = c(-0.2, 1.2),
                             custom_theme = NULL) {
  
  ioi_proteoforms <- e_set[fData(e_set)$ioi == ioi, ]
  
  if (nrow(ioi_proteoforms) == 0) stop("No data found for this ID of interest.")
  
  data <- ioi_proteoforms %>%
    tidy.ExpressionSet(addPheno = TRUE) %>%
    set_names(gsub("^gene$", "proteoform_id", colnames(.))) %>%
    left_join(fData(ioi_proteoforms), by = "proteoform_id") %>%
    mutate(membership = as.factor(membership))
  
  data$x_col <- data[[x_col]]
  data$color_col <- data[[color_col]]
  
  p <- data %>%
    ggplot(aes(x = x_col,
               y = value,
               group = interaction(proteoform_id, sample_name),
               color = color_col)) +
    geom_point() +
    scale_color_manual(values = colors, name = color_col) +
    ylim(y_limits)
  
  if (!hide_lines) {
    p <- p + geom_line(alpha = .7)
  }
  
  if (add_splines) {
    p <- p + geom_smooth(method = "lm", formula = 'y ~ splines::ns(x, df = 4)', se = FALSE, alpha = .4, size = .7)
  }
  
  if (!is.null(facet_col)) {
    p <- p + facet_wrap(facets = facet_col)
  }
  
  if (!is.null(custom_theme)) {
    p <- p + custom_theme()
  }
  
  p <- p +
    xlab(x_label) + ylab(y_label) +
    ggtitle(ioi) +
    theme(legend.position = "bottom")
  
  print(p)
  
  return(p)
}