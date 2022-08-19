#' Compute F statistic for PPI coaggregation differences across cell lines
#' 
#' @param fc_mat_list list of fold changes matrices from different cell lines
#' with temperature attribute
#' @param ppi_anno data frame of annotating known PPIs from e.g. StringDb
#' @param min_n minimal number of cell lines to compare
#' @param BPPARAM BiocParallel variable defining whether
#' to parallelize NPARC fits
#' 
#' @return data frame with F-statistic for PPI coaggregation group differences
#' 
#' @import Rtpca
#' 
#' @export
compute_Fstat_PPI_coaggregation_differences <- function(
    fc_mat_list, ppi_anno, min_n = 3,
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)){

    # loop over all fc mats / cell lines and create distance matrices
    message("Computing distance matrices...")
    dist_mat_list <- lapply(fc_mat_list, function(fc_mat){
        dist_mat <- createDistMat(list(fc_mat))
        Rtpca:::.filterDistMat(
            dist_mat,
            ppi_anno)
    })
    
    # loop over all dist matrices and create annotated ppi df
    message("Computing F-statistics...")
    dist_combo_df <- bind_rows(
        bplapply(dist_mat_list, 
                 BPPARAM = BPPARAM,
                 function(dist_mat){
                     Rtpca:::.distMat2AnnotatedPPiDf(
                         dist_mat,
                         ppi_anno
                         )
                     })) %>% 
        group_by(pair) %>% 
        filter(n() >= min_n) %>% 
        dplyr::summarize(min_value = min(value, na.rm = TRUE),
                  max_value = max(value, na.rm = TRUE),
                  .groups = "keep") %>% 
        ungroup %>% 
        mutate(min_rss = min_value^2,
               max_rss = max_value^2,
               f_stat = (max_rss - min_rss)/min_rss)
    
    return(dist_combo_df)
}

#' Compute robust F statistic for PPI coaggregation differences across cell lines
#' 
#' @param fc_mat_list list of fold changes matrices from different cell lines
#' with temperature attribute
#' @param ppi_anno data frame of annotating known PPIs from e.g. StringDb
#' @param min_n minimal number of cell lines to compare
#' @param BPPARAM BiocParallel variable defining whether
#' to parallelize NPARC fits
#' 
#' @return data frame with F-statistic for PPI coaggregation group differences
#' 
#' @import Rtpca
#' 
#' @export
compute_robust_Fstat_PPI_coaggregation_differences <- function(
    fc_mat_list, ppi_anno, min_n = 3,
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)){
    
    # loop over all fc mats / cell lines and create distance matrices
    message("Computing distance matrices...")
    dist_mat_list <- lapply(fc_mat_list, function(fc_mat){
        dist_mat <- createDistMat(list(fc_mat))
        Rtpca:::.filterDistMat(
            dist_mat,
            ppi_anno)
    })
    
    # loop over all dist matrices and create annotated ppi df
    message("Computing F-statistics...")
    dist_combo_df <- bind_rows(
        bplapply(dist_mat_list, 
                 BPPARAM = BPPARAM,
                 function(dist_mat){
                     Rtpca:::.distMat2AnnotatedPPiDf(
                         dist_mat,
                         ppi_anno
                     )
                 })) %>% 
        group_by(pair) %>% 
        filter(n() >= min_n) %>% 
        dplyr::summarize(min_value = sort(value)[2],
                  max_value = sort(value)[length(value) - 1],
                  .groups = "keep") %>% 
        ungroup %>% 
        mutate(min_rss = min_value^2,
               max_rss = max_value^2,
               f_stat = (max_rss - min_rss)/min_rss)
    
    return(dist_combo_df)
}
