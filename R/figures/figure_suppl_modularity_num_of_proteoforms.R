library(tidyverse)
library(biomaRt)
library(Biostrings)
library(Biobase)
library(readxl)
library(here)

theme_paper <- theme_bw(base_size = 6) +
    theme(legend.background = element_blank(), 
          legend.key = element_blank(), 
          panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 8),
          plot.background = element_blank(), 
          complete = TRUE,
          axis.line = element_line(color = "black", size = 0.25),
          text = element_text(size = 8),
          axis.ticks = element_line(color = "black", size = 0.25),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8),
          legend.text = element_text(size = 8))

# define proteoform output folder
proteoform_detection_folder <- here("proteoform_detection/output/standard/")

# read in peptide data
peptides <- readRDS(file.path(proteoform_detection_folder, "peptides.RDS"))
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

# read in proteoform data
proteoforms <- readRDS(file.path(proteoform_detection_folder, 
                                 "proteoforms_narrow_range_focused.RDS"))
proteoform_df <- biobroom::tidy.ExpressionSet(
    proteoforms, 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

number_of_proteoforms_df <- proteoform_df %>% 
    summarize(proteoform = unique(gene)) %>% 
    mutate(gene = sub("_.+", "", proteoform)) %>% 
    group_by(gene) %>% 
    summarize(n = n()) %>% 
    ungroup 

num_of_proteoforms_df <- read_csv(
    file.path(proteoform_detection_folder, "num_of_components.csv"))

number_of_proteoforms_df <- number_of_proteoforms_df %>% 
    left_join(num_of_proteoforms_df, by = c("gene" = "key")) %>% 
    mutate(enough_data = ifelse(is.na(value), FALSE, TRUE)) %>% 
    arrange(desc(n), desc(enough_data)) %>% 
    filter(!duplicated(gene)) %>% 
    mutate(rank = 1:n())

ggplot(number_of_proteoforms_df, aes(rank, n)) + 
    geom_point(aes(color = enough_data)) +
    labs(y = "Number of proteoforms per gene",
         x = "Gene rank") +
    scale_color_manual("Enough peptides for proteoform detection",
                       values = c("gray", "black")) +
    theme_paper +
    theme(legend.position = "bottom")

ggsave(filename = "number_of_proteoforms_freq.pdf", 
       width = 21, height = 12, units = "cm")

modularity_comparison_df <- read_delim(
    file.path(proteoform_detection_folder, "modularity_comparison.txt"), delim = "\t")

modularities_num_proteoforms_df <- modularity_comparison_df %>% 
    left_join(num_of_proteoforms_df, by = c("gene_symbols" = "key"))

ggplot(modularities_num_proteoforms_df, aes(x_modularity, y_modularity)) +
    geom_point(aes(color = as.factor(value)), alpha = 0.5) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    #ggrepel::geom_text_repel() +
    ggpmisc::stat_dens2d_filter(aes(label = gene_symbols), geom = "text_repel", 
                                keep.fraction = 0.01, color = "gray20") +
    scale_color_discrete("Number of detected proteoforms") +
    labs(x = "Modularity by melting curve cluster membership",
         y = "Modularity by protein ID") +
    theme_paper +
    theme(legend.position = "bottom")

ggsave(filename = "modularity_comparison.pdf", 
       width = 21, height = 12, units = "cm")