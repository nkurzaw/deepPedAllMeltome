library(tidyverse)
library(Biobase)
library(igraph)
library(here)
library(biobroom)

# define proteoform detection output folder
proteoform_detection_folder <- here("proteoform_detection/output/standard/")

# load peptides
peptides <- readRDS(
    file.path(proteoform_detection_folder, "peptides.RDS"))

# convert peptide level data to data frame    
peptides_df <- biobroom::tidy.ExpressionSet(peptides) %>% 
    left_join(fData(peptides) %>% 
                  as.data.frame() %>%
                  rownames_to_column() %>% 
                  as_tibble() %>%                  
                  dplyr::select(gene = rowname, id, first_protein_id),
              by = "gene") %>% 
    left_join(pData(peptides) %>% 
                  as_tibble() %>% 
                  dplyr::select(sample = sample_id, temperature, sample_name_machine),
              by = "sample")

# load graphs
graphs <- readRDS(file.path(proteoform_detection_folder, "graphs_comms.RDS"))

# filter graphs for 0 modularity
graphs_gt0 <- graphs[(lapply(graphs, get.graph.attribute, 
                            name = "proteoform_modularity") > 0) %>% unlist()]

graphs_lte0 <- graphs[(lapply(graphs, get.graph.attribute, 
                            name = "proteoform_modularity") <= 0) %>% unlist()]

# load proteoform ratios
proteoforms <- readRDS(
    file.path(proteoform_detection_folder, "proteoforms_narrow_range_focused.RDS"))

# create proteoform data frame
proteoform_df <- biobroom::tidy.ExpressionSet(
    proteoforms, 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

# convert graphs_gt0 to data frame
graph_gt0_df <- bind_rows(lapply(names(graphs_gt0), function(gr_nm){
    membership_vec <- graphs_gt0[[gr_nm]]$communities$membership
    tibble(gene = gr_nm,
           psms = names(membership_vec),
           proteoform = membership_vec,
           proteoform_id = paste(gr_nm, membership_vec, sep = "_")) %>% 
        arrange(proteoform) %>% 
        group_by(gene, proteoform, proteoform_id) %>% 
        summarize(peptides = toString(psms), .groups = "keep") %>% 
        ungroup()
})) %>% 
    filter(proteoform_id %in% proteoform_df$gene) %>% 
    mutate(group = "modularity > 1e-13")

# convert graphs_lte0 to data frame
graph_lte0_df <- bind_rows(lapply(names(graphs_lte0), function(gr_nm){
    membership_vec <- graphs_lte0[[gr_nm]]$communities$membership
    tibble(gene = gr_nm,
           psms = names(membership_vec),
           proteoform = 0,
           proteoform_id = paste(gr_nm, 0, sep = "_")) %>% 
        group_by(gene, proteoform, proteoform_id) %>% 
        summarize(peptides = toString(psms), .groups = "keep") %>% 
        ungroup()
})) %>% 
    filter(proteoform_id %in% proteoform_df$gene) %>% 
    mutate(group = "modularity <= 1e-13")

# get peptides of proteoforms with too few peptides for proteoform detection
not_yet_covered_peptides <- unique(
    proteoform_df %>% 
        filter(!gene %in% c(graph_gt0_df$proteoform_id, graph_lte0_df$proteoform_id)) %>% 
        pull(gene))

# make data frame with remaining psm assignments
non_graph_df <- peptides_df %>% 
    dplyr::select(gene = id, psms = gene) %>% 
    filter(gene %in% sub("_.+", "", not_yet_covered_peptides)) %>% 
    mutate(proteoform = 0, proteoform_id = paste(gene, "0", sep = "_")) %>% 
    group_by(gene, proteoform, proteoform_id) %>% 
    summarize(peptides = toString(psms), .groups = "keep") %>% 
    ungroup() %>% 
    mutate(group = "proteoform detection criteria not met")

# assemble full table
full_suppl_df <- bind_rows(
    graph_gt0_df,
    graph_lte0_df,
    non_graph_df
)

# get cell line fold changes
proteofrom_spread_df <- proteoform_df %>% 
    dplyr::select(proteoform_id = gene,
                  sample_name_machine,
                  temperature, 
                  rel_value) %>% 
    unite(new_col, sample_name_machine, temperature) %>% 
    spread(new_col, rel_value)

# join annotation and fold changes
full_suppl_fc_df <- left_join(
    full_suppl_df, proteofrom_spread_df,
    by = "proteoform_id")

write_delim(full_suppl_fc_df, file = here("R/tables/suppl_table_2_proteoform_detection.txt"),
            delim = "\t")
