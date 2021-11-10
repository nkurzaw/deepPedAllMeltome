library(tidyverse)
library(here)

# get sample meta data
sample_meta_raw <- read_tsv(here("meta/sample_meta.txt"))
sample_meta_stages <- read_tsv(here("meta/sample_meta_stages.txt"))

# join sample meta files
full_sample_meta <- left_join(
    sample_meta_raw %>% 
        dplyr::select(tmt_set = set,
                      labels,
                      sample_name,
                      subtype,
                      sample_name_machine),
    sample_meta_stages %>% 
        dplyr::select(sample_name,
                      sex = gender,
                      developmental_stage = Assigned_stage), 
    by = "sample_name") %>% 
    mutate(subtype = sub("\\.", "-", subtype)) %>% 
    within(subtype[subtype == "PTMA-TMSB4X"] <- "other") %>%
    within(developmental_stage[developmental_stage == "pre-B other"] <- "pre-B") %>%
    mutate(developmental_stage = factor(
        developmental_stage, levels = c("pre-pro-B", "pro-B", "early pre-B", 
                                        "pre-B", "late pre-B"))) %>% 
    mutate(sample_name = gsub("_", "-", sample_name_machine)) %>% 
    dplyr::select(tmt_set, labels, sample_name, subtype, sex, developmental_stage)
    
write_tsv(full_sample_meta, path = here("R/tables/suppl_table_full_sample_meta.txt"))
