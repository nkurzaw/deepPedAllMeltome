#' Plot the number of proteoforms per id
#'
#' @param graphs graphs
#' @param BPPARAM BiocParallel
#'
#' @return plot
#' 
#' @importFrom BiocParallel bplapply
#' @importFrom magrittr "%>%"
#' @importFrom igraph V
#' @import ggplot2
#' 
#' @export
#'
#' @examples
plot_num_proteoforms_per_id <- function (graphs,
                                         BPPARAM = BiocParallel::SerialParam()) {
  
  num_proteoforms <- bplapply(X = graphs,
                              FUN = function (g) {
                                g %>%
                                  V() %>%
                                  .$membership %>%
                                  unique() %>%
                                  length()
                              },
                              BPPARAM = BPPARAM)
  
  p <- num_proteoforms %>%
    unlist() %>%
    data.frame(num_proteoforms = ., stringsAsFactors = FALSE) %>%
    arrange(desc(num_proteoforms)) %>%
    mutate(rank = seq_len(nrow(.)) %>% as.numeric()) %>%
    ggplot(aes(x = rank, y = num_proteoforms)) +
    geom_point(size = 3, stroke = 0, alpha = .7, color = "steelblue") +
    xlab("ID rank") + ylab("Proteoforms per ID")
  
  print(p)
}