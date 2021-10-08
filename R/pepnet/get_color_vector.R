#' Create a color vector
#'
#' @param colors color pallette
#' @param vec vector of labels
#' @param gray standard gray color for NA values
#'
#' @return named vector
#' 
#' @importFrom magrittr "%>%"
#' 
#' @export
#'
#' @examples
get_color_vector <- function (colors,
                              vec,
                              gray = "#494e57") {
  
  n <- length(unique(vec))
  
  suppl_needed <- max(n - length(colors), 0)
  col_suppl <- c(colors[1:min(n, length(colors))], rep(gray, suppl_needed))
  
  # count appearances of each group
  vec_accum <- table(vec) %>% sort(decreasing = TRUE)
  
  names(col_suppl) <- names(vec_accum) %>% as.character()
  
  # add a gray missing value
  col_suppl <- c(col_suppl, c("NA" = gray))
  
  return(col_suppl)
}
