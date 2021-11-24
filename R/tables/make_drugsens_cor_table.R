library(here)
library(tidyverse)

# define drugsens cor analysis output folder
drugsens_cor_folder <- here("drugsens_cor/output/")

# read in multi-cell line tpca data frame
limma_out_df <- readRDS(
    file.path(drugsens_cor_folder, "limma_out_df.RDS"))

# make final table
drugsens_cor_out_df <- limma_out_df %>% 
    dplyr::select(drug_name = rowname, 
                  tested_proteoform = testedProt, 
                  average_drugsens = AveExpr,
                  lm_slope = logFC,
                  t_statistic = t,
                  p_value = P.Value,
                  p_adj) %>% 
    mutate(significant = case_when(
        p_adj < 0.1 ~ TRUE,
        TRUE ~ FALSE
    ))

# write table
write_csv(drugsens_cor_out_df, here("R/tables/suppl_table_drugsens_cor.csv"))
