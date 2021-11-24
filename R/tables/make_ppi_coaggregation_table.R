library(here)
library(tidyverse)

# define ppi coaggregation analysis output folder
ppi_coaggregation_folder <- here("ppi_coaggregation/output")

# read in multi-cell line tpca data frame
multi_cell_line_rtpca_robust_df <- readRDS(
    file.path(ppi_coaggregation_folder, "multi_cell_line_rtpca_robust_df.RDS"))

# make final table
multi_tpca_out_df <- multi_cell_line_rtpca_robust_df %>% 
    dplyr::select(proteoform_pair = pair, min_rss, max_rss, F_statistic = f_stat) %>% 
    mutate(significant_differential_coaggregation = case_when(
        F_statistic < quantile(F_statistic, 0.9) ~ TRUE,
        TRUE ~ FALSE
    ))

# write table
write_csv(multi_tpca_out_df, here("R/tables/suppl_table_ppi_coaggregation.csv"))
