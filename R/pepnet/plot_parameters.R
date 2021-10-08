plot_parameters <- function (parameters_df,
                             output_folder = NULL) {
  
  if (!is.null(output_folder) && !dir.exists(output_folder)) dir.create(path = output_folder, recursive = TRUE)
  
  plots <- lapply(X = colnames(parameters_df),
                  FUN = function (parameter_name) {
                    ggplot(parameters_df) +
                      geom_density(aes_string(x = parameter_name))
                  })
  if (!is.null(output_folder)) pdf(file = file.path(output_folder, "community_parameters.pdf"))
  print(cowplot::plot_grid(plotlist = plots))
  if (!is.null(output_folder)) dev.off()
}