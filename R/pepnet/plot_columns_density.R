#' Plot columns density, e.g. quality parameters
#'
#' @param df \code{data.frame}
#'
#' @return plot
#' 
#' @importFrom cowplot plot_grid
#' @import ggplot2
#' 
#' @export
#'
#' @examples
plot_columns_density <- function (df) {
  
  plots <- lapply(X = colnames(df),
                  FUN = function (parameter_name) {
                    ggplot(df) +
                      geom_density(aes_string(x = parameter_name))
                  })
  
  print(plot_grid(plotlist = plots))
}