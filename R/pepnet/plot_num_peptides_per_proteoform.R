#' Plot the number of peptides per proteoform
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
plot_num_peptides_per_proteoform <- function (graphs,
                                              BPPARAM = BiocParallel::SerialParam()) {
  
  num_peptides <- bplapply(X = graphs,
                           FUN = function (g) {
                             g %>%
                               V() %>%
                               .$membership %>%
                               table()
                           },
                           BPPARAM = BPPARAM)
  
  p <- num_peptides %>%
    unlist() %>%
    data.frame(num_peptides = ., stringsAsFactors = FALSE) %>%
    arrange(desc(num_peptides)) %>%
    mutate(rank = seq_len(nrow(.)) %>% as.numeric()) %>%
    ggplot(aes(x = rank, y = num_peptides)) +
    geom_point(size = 3, stroke = 0, alpha = .7, color = "steelblue") +
    scale_y_log10() +
    xlab("ID rank") + ylab("Peptides per proteoform")
  
  print(p)
}