library(tidyverse)
library(here)

# define nparc output folder
nparc_folder <- here("nparc/output/standard/")

# read in NPARC result tables
nparc_res_hq_df <- readRDS(file.path(nparc_folder, "nparc_res_hq_df.RDS"))
nparc_fstat_df <- readRDS(file.path(nparc_folder, "nparc_fstat_df.RDS"))

# select relevant columns
nparc_fstat_sub_df <- nparc_fstat_df %>% 
    dplyr::select(proteoform_id = id,
                  n_fitted = nFitted,
                  rss_null = rss,
                  rss_alternative = rssAlternative,
                  d1, d2,
                  F_statistic)

# get aucs in different cell lines
auc_df <- nparc_res_hq_df %>% 
    dplyr::select(sample = sample_name, id, aumc) %>% 
    mutate(sample = paste(sample, "AUMC", sep = "_")) %>% 
    spread(sample, aumc)

# join tables
nparc_full_suppl_df <- nparc_fstat_sub_df %>% 
    left_join(auc_df, by = c("proteoform_id" = "id")) %>% 
    mutate(significant = case_when(
        F_statistic >= quantile(nparc_fstat_df$F_statistic, 0.9) ~ TRUE,
        TRUE ~ FALSE))

write_tsv(nparc_full_suppl_df, path = here("R/tables/suppl_table_nparc.txt"))
