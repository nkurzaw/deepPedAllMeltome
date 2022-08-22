library(tidyverse)
library(Biobase)
library(readxl)
library(cowplot)
library(here)

# define plotting style for manuscript
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

proteoforms <- readRDS(
    here("proteoform_detection/output/standard/proteoforms_narrow_range_focused.RDS"))

proteoform_fdata <- fData(proteoforms)[,1:6]

proteoform_fdata_gene_level <- proteoform_fdata %>% 
    group_by(ioi) %>% 
    dplyr::summarize(num_peptides = sum(num_peptides),
              num_proteoforms = max(membership),
              .groups = "drop") %>% 
    within(num_proteoforms[num_proteoforms == 0] <- 1)

uniprot_features <- read_tsv(here("data/uniprot_human_proteome_features.tsv")) %>% 
    mutate(gene_name = sub(" .+", "", `Gene Names`),
           known_isoforms = sub(";.+", "", sub(".+Named isoforms=", "", `Alternative products (isoforms)`)),
           cc_loc = case_when(grepl("embrane", `Subcellular location [CC]`) ~ "Membrane",
                              grepl("Nucleus", `Subcellular location [CC]`) ~ "Nucleus",
                              grepl("Mitochondrion", `Subcellular location [CC]`) ~ "Mitochondrion",
                              grepl("Golgi|Endoplasmatic|ER", `Subcellular location [CC]`) ~ "Golgi/ER",
                              grepl("Cytoplasm", `Subcellular location [CC]`) ~ "Cytoplasm",
                              TRUE ~ "Other"),
           small_molecule_binding = grepl("small molecule binding|vitamin binding|ATP binding|GTP binding|alcohol binding|organic acid binding",
                                          `Gene Ontology (molecular function)`)) %>% 
    within(known_isoforms[known_isoforms == "ALTERNATIVE PRODUCTS:"] <- 1) %>% 
    dplyr::select(uniprot_id = Entry, gene_name, length = Length, known_isoforms, 
                  cc_loc, small_molecule_binding) %>% 
    group_by(gene_name) %>% 
    summarize(length = max(length),
              known_isoforms = max(as.numeric(known_isoforms)),
              cc_loc = case_when(any(grepl("Membrane", cc_loc)) ~ "Membrane",
                                 any(grepl("Nucleus", cc_loc)) ~ "Nucleus",
                                 any(grepl("Mitochondrion", cc_loc)) ~ "Mitochondrion",
                                 any(grepl("Golgi|Endoplasmatic|ER", cc_loc)) ~ "Golgi/ER",
                                 any(grepl("Cytoplasm", cc_loc)) ~ "Cytoplasm",
                                 TRUE ~ "Other"),
              small_molecule_binding = any(small_molecule_binding),
              .groups = "drop")

halflife_df <- read_xlsx(here("data/41467_2018_3106_MOESM5_ESM.xlsx")) %>% 
    dplyr::select(gene_name, matches("Bcells.+half_life")) %>% 
    na.omit() %>% 
    gather(key, value, -gene_name) %>% 
    group_by(gene_name) %>% 
    summarise(half_life_log10 = log10(mean(value))) %>% 
    ungroup()

proteoform_fdata_gene_level <- left_join(proteoform_fdata_gene_level, uniprot_features, 
                                         by = c("ioi" = "gene_name")) %>% 
    left_join(halflife_df, by = c("ioi" = "gene_name"))

# plot number of found proteoforms stratified by known proteoforms
proteoform_fdata_gene_level %>% 
    mutate(known_isoforms_binned = cut(known_isoforms, breaks = c(0, 1, 2, 4, 6, 32))) %>% 
    filter(!is.na(known_isoforms_binned)) %>% 
    ggplot(aes(num_proteoforms)) + 
    geom_bar() + 
    facet_wrap(~known_isoforms_binned)

# plot number of found proteoforms stratified by cellular location
proteoform_fdata_gene_level %>% 
    filter(!is.na(cc_loc)) %>% 
    ggplot(aes(num_proteoforms)) + 
    geom_bar() + 
    facet_wrap(~cc_loc)

# plot number of found proteoforms stratified by cellular location
proteoform_fdata_gene_level %>% 
    filter(!is.na(cc_loc)) %>% 
    filter(!is.na(known_isoforms)) %>% 
    ggplot(aes(known_isoforms)) + 
    geom_bar() + 
    facet_wrap(~cc_loc)

# plot number of found proteoforms vs. turnover in B-cells
proteoform_fdata_gene_level %>% 
    filter(!is.na(half_life_log10)) %>% 
    ggplot(aes(as.factor(num_proteoforms), half_life_log10)) +
    geom_violin() +
    geom_boxplot(width = 0.15, outlier.colour = NA)

proteoform_fdata_gene_level %>% 
    filter(!is.na(half_life_log10)) %>% 
    mutate(known_isoforms_binned = cut(known_isoforms, breaks = c(0, 1, 2, 4, 6, 32))) %>% 
    filter(!is.na(known_isoforms_binned)) %>% 
    ggplot(aes(as.factor(known_isoforms_binned), half_life_log10)) +
    geom_violin() +
    geom_boxplot(width = 0.15, outlier.colour = NA)

# plot number of found proteoforms vs. protein length
proteoform_fdata_gene_level %>% 
    filter(!is.na(length)) %>% 
    ggplot(aes(as.factor(num_proteoforms), log2(length))) +
    geom_violin() +
    geom_boxplot(width = 0.15, outlier.colour = NA)

proteoform_fdata_gene_level %>% 
    filter(!is.na(num_peptides)) %>% 
    ggplot(aes(as.factor(num_proteoforms), log2(num_peptides))) +
    geom_violin() +
    geom_boxplot(width = 0.15, outlier.colour = NA)


proteoform_fdata_gene_level %>% 
    filter(!is.na(length)) %>% 
    ggplot(aes(as.factor(num_proteoforms), log2(num_peptides/length))) +
    geom_violin() +
    geom_boxplot(width = 0.15, outlier.colour = NA)

proteoform_fdata_gene_level %>% 
    filter(!is.na(cc_loc)) %>% 
    ggplot(aes(cc_loc, num_proteoforms)) +
    geom_violin() +
    geom_boxplot(width = 0.15, outlier.colour = NA)

