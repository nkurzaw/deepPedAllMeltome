library(tidyverse)
library(readxl)
library(limma)
library(here)

# define data path
data_path <- here("data/eps8l2_ip/")

# read in metadata
metadata_df <- read_xlsx(list.files(data_path, full.names = TRUE, pattern = "metadata")) %>% 
    dplyr::select(label, lysate, sample, data_file = `text file`)

# read in MS data
all_df <- bind_rows(lapply(list.files(data_path, full.names = TRUE, pattern = "MSdata"), 
                            function(file_path){
    readin_df <- read_tsv(file_path) %>% 
        dplyr::select(protein_id, gene_name, matches("signal_sum")) %>% 
        filter(!grepl("##", protein_id)) %>% 
        gather(key, value, -protein_id, -gene_name) %>% 
        mutate(label = sub("signal_sum_", "", key),
               data_file = basename(file_path))
}))

# join meta data
full_df <- left_join(all_df, metadata_df, by = c("data_file", "label"))

# create data matrix
full_mat_df <- full_df %>% 
    filter(value > 0, !is.na(value), !is.na(lysate), !is.na(sample)) %>% 
    dplyr::select(-key, -label, -data_file) %>% 
    unite("sample_name", lysate, sample, sep = "_") %>% 
    spread(sample_name, value)

# create expression set
