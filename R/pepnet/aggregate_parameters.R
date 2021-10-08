#' Aggregate parameters
#'
#' @param parameters \code{list} of parameters
#'
#' @return \code{data.frame} of parameters
#' 
#' @importFrom magrittr "%>%"
#' @importFrom purrr set_names
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr full_join
#' 
#' @export
#'
#' @examples
aggregate_parameters <- function (parameters) {
  
  # aggregate parameters in a dataframe
  parameter_names <- parameters[[1]] %>% names()
  
  parameters_df <- lapply(X = parameter_names,
                          FUN = function (parameter_name) {
                            lapply(X = parameters,
                                   FUN = function (param) param[[parameter_name]]) %>%
                              unlist() %>%
                              data.frame(parameter_name = ., stringsAsFactors = FALSE) %>%
                              set_names(parameter_name) %>%
                              rownames_to_column("item_id")
                          }) %>%
    Reduce(f = function (x, y) full_join(x, y, by = "item_id"), x = .)
  
  return(parameters_df)
}