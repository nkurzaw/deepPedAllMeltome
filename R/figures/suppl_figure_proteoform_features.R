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

proteoform_fdata <- fData(proteoforms)[,c(1:7, grep("126_sd|134N_sd", colnames(fData(proteoforms))))]

proteoform_fdata_gene_level <- proteoform_fdata %>% 
    group_by(ioi) %>% 
    dplyr::summarize(num_peptides = sum(num_peptides),
                     num_proteoforms = max(membership),
                     RCH_ACV_abundance_log2 = log2(sum(`__quant_Set1_tmt16plex_126_sd`, na.rm = TRUE)),
                     LC4_1_abundance_log2 = log2(sum(`__quant_Set1_tmt16plex_134N_sd`, na.rm = TRUE)),
                     REH_abundance_log2 = log2(sum(`__quant_Set2_tmt16plex_126_sd`, na.rm = TRUE)),
                     P30_OHKUBO_abundance_log2 = log2(sum(`__quant_Set2_tmt16plex_134N_sd`, na.rm = TRUE)),
                     KASUMI_9_abundance_log2 = log2(sum(`__quant_Set3_tmt16plex_126_sd`, na.rm = TRUE)),
                     SEM_abundance_log2 = log2(sum(`__quant_Set3_tmt16plex_134N_sd`, na.rm = TRUE)),
                     COG_355_abundance_log2 = log2(sum(`__quant_Set4_tmt16plex_126_sd`, na.rm = TRUE)),
                     X697_abundance_log2 = log2(sum(`__quant_Set4_tmt16plex_134N_sd`, na.rm = TRUE)),
                     KASUMI_2_abundance_log2 = log2(sum(`__quant_Set5_tmt16plex_126_sd`, na.rm = TRUE)),
                     ALL_PO_abundance_log2 = log2(sum(`__quant_Set5_tmt16plex_134N_sd`, na.rm = TRUE)),
                     SUP_15_abundance_log2 = log2(sum(`__quant_Set6_tmt16plex_126_sd`, na.rm = TRUE)),
                     MHH_CALL_3_abundance_log2 = log2(sum(`__quant_Set6_tmt16plex_134N_sd`, na.rm = TRUE)),
                     TMD5_abundance_log2 = log2(sum(`__quant_Set7_tmt16plex_126_sd`, na.rm = TRUE)),
                     KOPN_8_abundance_log2 = log2(sum(`__quant_Set7_tmt16plex_134N_sd`, na.rm = TRUE)),
                     NALL_1_abundance_log2 = log2(sum(`__quant_Set8_tmt16plex_126_sd`, na.rm = TRUE)),
                     COG_319_abundance_log2 = log2(sum(`__quant_Set8_tmt16plex_134N_sd`, na.rm = TRUE)),
                     MHH_CALL_2_abundance_log2 = log2(sum(`__quant_Set9_tmt16plex_126_sd`, na.rm = TRUE)),
                     COG_394_abundance_log2 = log2(sum(`__quant_Set9_tmt16plex_134N_sd`, na.rm = TRUE)),
                     HAL_01_abundance_log2 = log2(sum(`__quant_Set10_tmt16plex_126_sd`, na.rm = TRUE)),
                     MHH_CALL_4_abundance_log2 = log2(sum(`__quant_Set10_tmt16plex_134N_sd`, na.rm = TRUE)),
                     .groups = "drop") %>% 
    within(num_proteoforms[num_proteoforms == 0] <- 1)

functional_phosphosites <- read_xlsx(here("data/41587_2019_344_MOESM5_ESM.xlsx")) %>% 
    group_by(uniprot) %>% 
    summarize(phosphosite_count = n(),
              max_functional_score = max(functional_score),
              .groups = "drop")

disorder_predictions <- read_tsv(here("data/all_uniprot_instrinsically_disordered_annotation_swissprot_march21_min_5_aas_full.txt")) %>% 
    mutate(rank = as.numeric(rank)) %>% 
    dplyr::select(uniprot_id = protein, fraction_disordered = rank)

uniprot_features <- read_tsv(here("data/uniprot_human_proteome_features.tsv")) %>% 
    mutate(gene_name = sub(" .+", "", `Gene Names`),
           known_isoforms = sub(";.+", "", sub(".+Named isoforms=", "", `Alternative products (isoforms)`)),
           cc_loc = case_when(grepl("embrane", `Subcellular location [CC]`) ~ "Membrane",
                              grepl("Nucleus", `Subcellular location [CC]`) ~ "Nucleus",
                              grepl("Mitochondrion", `Subcellular location [CC]`) ~ "Mitochondria",
                              grepl("Golgi|Endoplasmatic|ER", `Subcellular location [CC]`) ~ "Golgi/ER",
                              grepl("Cytoplasm", `Subcellular location [CC]`) ~ "Cytoplasm",
                              TRUE ~ "Other"),
           small_molecule_binding = grepl("small molecule binding|vitamin binding|ATP binding|GTP binding|alcohol binding|organic acid binding",
                                          `Gene Ontology (molecular function)`)) %>% 
    within(known_isoforms[known_isoforms == "ALTERNATIVE PRODUCTS:"] <- 1) %>% 
    dplyr::select(uniprot_id = Entry, gene_name, length = Length, known_isoforms, 
                  cc_loc, small_molecule_binding) %>% 
    left_join(functional_phosphosites, by = c("uniprot_id" = "uniprot")) %>% 
    left_join(disorder_predictions, by = "uniprot_id") %>% 
    group_by(gene_name) %>% 
    summarize(length = max(length),
              known_isoforms = max(as.numeric(known_isoforms)),
              cc_loc = case_when(any(grepl("Membrane", cc_loc)) ~ "Membrane",
                                 any(grepl("Nucleus", cc_loc)) ~ "Nucleus",
                                 any(grepl("Mitochondria", cc_loc)) ~ "Mitochondria",
                                 any(grepl("Golgi|Endoplasmatic|ER", cc_loc)) ~ "Golgi/ER",
                                 any(grepl("Cytoplasm", cc_loc)) ~ "Cytoplasm",
                                 TRUE ~ "Other"),
              small_molecule_binding = any(small_molecule_binding),
              phosphosite_count = sum(phosphosite_count),
              max_functional_score = max(max_functional_score),
              fraction_disorder = max(fraction_disordered, na.rm = TRUE),
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
detected_vs_known_proteoforms <- proteoform_fdata_gene_level %>% 
    mutate(known_isoforms_binned = cut(known_isoforms, breaks = c(0, 1, 2, 4, 6, 32))) %>% 
    filter(!is.na(known_isoforms_binned)) %>% 
    group_by(known_isoforms_binned) %>% 
    mutate(median_num_proteoforms = median(num_proteoforms),
              .groups = "drop") %>% 
    ggplot(aes(num_proteoforms)) + 
    geom_bar() + 
    #geom_vline(aes(xintercept = median_num_proteoforms), color = "orange") +
    facet_wrap(~known_isoforms_binned) +
    labs(x = "Number of detected proteoforms",
         y = "Count by bin of number known isoforms") +
    #ggtitle("Number of detected proteoforms per gene stratified by number of known isoforms") +
    theme_paper

# plot number of found proteoforms stratified by cellular location
cc_loc_detected_proteoforms <- proteoform_fdata_gene_level %>% 
    filter(!is.na(cc_loc)) %>% 
    group_by(cc_loc) %>% 
    mutate(median_num_proteoforms = median(num_proteoforms),
           .groups = "drop") %>% 
    ggplot(aes(num_proteoforms)) + 
    geom_bar() + 
    #geom_vline(aes(xintercept = median_num_proteoforms), color = "orange") +
    facet_wrap(~cc_loc) +
    labs(x = "Number of detected proteoforms",
         y = "Count by bin of cellular location") +
    #ggtitle("Number of detected proteoforms per gene stratified by cellular location") +
    theme_paper

# plot number of found proteoforms stratified by cellular location
cc_loc_known_isoforms <- proteoform_fdata_gene_level %>% 
    filter(!is.na(cc_loc)) %>% 
    filter(!is.na(known_isoforms)) %>% 
    group_by(cc_loc) %>% 
    mutate(median_known_isoforms = median(known_isoforms),
           .groups = "drop") %>% 
    ggplot(aes(known_isoforms)) + 
    geom_bar(color = "lightgray", size = 0.05) + 
    #geom_vline(aes(xintercept = median_known_isoforms), color = "orange") +
    facet_wrap(~cc_loc) +
    labs(x = "Number of known isoforms",
         y = "Count by bin of cellular location") +
    #ggtitle("Number of known isoforms per gene stratified by cellular location") +
    theme_paper

# plot number of found proteoforms vs. turnover in B-cells
half_life_detected_proteoforms <- proteoform_fdata_gene_level %>% 
    filter(!is.na(half_life_log10)) %>% 
    ggplot(aes(as.factor(num_proteoforms), half_life_log10)) +
    geom_violin(fill = "gray", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.15, outlier.colour = NA) +
    geom_hline(aes(yintercept = median(half_life_log10)), 
               color = "black", alpha = 0.25, size = 2) +
    labs(x = "Number of detected proteoforms",
         y = "log10 half life") +
    #ggtitle("Half lives of proteins with different numbers of detected proteoforms") +
    theme_paper

half_life_known_isoforms <- proteoform_fdata_gene_level %>% 
    filter(!is.na(half_life_log10)) %>% 
    mutate(known_isoforms_binned = cut(known_isoforms, breaks = c(0, 1, 2, 4, 6, 32))) %>% 
    filter(!is.na(known_isoforms_binned)) %>% 
    ggplot(aes(as.factor(known_isoforms_binned), half_life_log10)) +
    geom_violin(fill = "gray", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.15, outlier.colour = NA) +
    geom_hline(aes(yintercept = median(half_life_log10)), 
               color = "black", alpha = 0.25, size = 2) +
    labs(x = "Binned number of known isoforms",
         y = "log10 half life") +
    #ggtitle("Half lives of proteins with different numbers of known isoforms") +
    theme_paper

# plot number of found proteoforms vs. protein length
length_detected_proteoforms <- proteoform_fdata_gene_level %>% 
    filter(!is.na(length)) %>% 
    ggplot(aes(as.factor(num_proteoforms), log2(length))) +
    geom_violin(fill = "gray", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.15, outlier.colour = NA) +
    geom_hline(aes(yintercept = median(log2(length))), 
               color = "black", alpha = 0.25, size = 2) +
    labs(x = "Number of detected proteoforms",
         y = "log2 max. isoform\nlength (AAs) per protein") +
    #ggtitle("Max. lengths of proteins with different numbers of detected proteoforms") +
    theme_paper

length_known_isoforms <- proteoform_fdata_gene_level %>% 
    filter(!is.na(length)) %>% 
    mutate(known_isoforms_binned = cut(known_isoforms, breaks = c(0, 1, 2, 4, 6, 32))) %>% 
    filter(!is.na(known_isoforms_binned)) %>% 
    ggplot(aes(as.factor(known_isoforms_binned), log2(length))) +
    geom_violin(fill = "gray", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.15, outlier.colour = NA) +
    geom_hline(aes(yintercept = median(log2(length))), 
               color = "black", alpha = 0.25, size = 2) +
    labs(x = "Binned number of known isoforms",
         y = "log2 max. isoform\nlength (AAs) per protein") +
    #ggtitle("Max. lengths of proteins with different known numbers of isoforms") +
    theme_paper

identified_peptides_detected_proteoforms <- proteoform_fdata_gene_level %>% 
    filter(!is.na(num_peptides)) %>% 
    ggplot(aes(as.factor(num_proteoforms), log2(num_peptides))) +
    geom_violin(fill = "gray", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.15, outlier.colour = NA) +
    geom_hline(aes(yintercept = median(log2(num_peptides))), 
               color = "black", alpha = 0.25, size = 2) +
    labs(x = "Number of detected proteoforms",
         y = "log2 number of unique quantified\npeptides identified per protein") +
    #ggtitle("Number of uniquely quantified peptides per gene with different numbers of detected proteoforms") +
    theme_paper

identified_peptides_known_isoforms <- proteoform_fdata_gene_level %>% 
    filter(!is.na(length)) %>% 
    mutate(known_isoforms_binned = cut(known_isoforms, breaks = c(0, 1, 2, 4, 6, 32))) %>% 
    filter(!is.na(known_isoforms_binned)) %>% 
    ggplot(aes(as.factor(known_isoforms_binned), log2(num_peptides))) +
    geom_violin(fill = "gray", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.15, outlier.colour = NA) +
    geom_hline(aes(yintercept = median(log2(num_peptides))), 
               color = "black", alpha = 0.25, size = 2) +
    labs(x = "Binned number of known isoforms",
         y = "log2 number of unique quantified\npeptides identified per protein") +
    #ggtitle("Number of uniquely quantified peptides per gene with different known numbers of isoforms") +
    theme_paper

norm_peptides_detected_proteoforms <- proteoform_fdata_gene_level %>% 
    filter(!is.na(length)) %>% 
    ggplot(aes(as.factor(num_proteoforms), log2(num_peptides/length))) +
    geom_violin(fill = "gray", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.15, outlier.colour = NA) +
    geom_hline(aes(yintercept = median(log2(num_peptides/length))), 
               color = "black", alpha = 0.25, size = 2) +
    labs(x = "Number of detected proteoforms",
         y = "log2 length normalized number of\nunique quantified peptides per protein") +
    #ggtitle("Protein length normalized number of uniquely quantified peptides per gene with different numbers of detected proteoforms") +
    theme_paper

norm_peptides_known_isoforms <- proteoform_fdata_gene_level %>% 
    filter(!is.na(length)) %>% 
    mutate(known_isoforms_binned = cut(known_isoforms, breaks = c(0, 1, 2, 4, 6, 32))) %>% 
    filter(!is.na(known_isoforms_binned)) %>% 
    ggplot(aes(as.factor(known_isoforms_binned), log2(num_peptides/length))) +
    geom_violin(fill = "gray", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.15, outlier.colour = NA) +
    geom_hline(aes(yintercept = median(log2(num_peptides/length))), 
               color = "black", alpha = 0.25, size = 2) +
    labs(x = "Binned number of known isoforms",
         y = "log2 length normalized number of\nunique quantified peptides per protein") +
    #ggtitle("Protein length normalized number of uniquely quantified peptides per gene with different numbers of detected proteoforms") +
    theme_paper

# proteoform_fdata_gene_level %>% 
#     filter(!is.na(cc_loc)) %>% 
#     ggplot(aes(cc_loc, num_proteoforms)) +
#     geom_violin() +
#     geom_boxplot(width = 0.15, outlier.colour = NA)

# phosphosite count
phosphosites_detected_proteoforms <- proteoform_fdata_gene_level  %>% 
    filter(!is.na(phosphosite_count)) %>% 
    ggplot(aes(as.factor(num_proteoforms), log2(phosphosite_count))) +
    geom_violin(fill = "gray", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.15, outlier.colour = NA) +
    geom_hline(aes(yintercept = median(log2(phosphosite_count))), 
               color = "black", alpha = 0.25, size = 2) +
    labs(x = "Number of detected proteoforms",
         y = "log2 number of phosphosites") +
    #ggtitle("Number of phosphosites") +
    theme_paper

# correct for length
norm_phosphosites_detected_proteoform <- proteoform_fdata_gene_level  %>% 
    filter(!is.na(phosphosite_count)) %>% 
    ggplot(aes(as.factor(num_proteoforms), log2(phosphosite_count/length))) +
    geom_violin(fill = "gray", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.15, outlier.colour = NA) +
    geom_hline(aes(yintercept = median(log2(phosphosite_count/length))), 
               color = "black", alpha = 0.25, size = 2) +
    labs(x = "Number of detected proteoforms",
         y = "log2 number of phosphosites\nnormalized by protein length") +
    #ggtitle("Number of phosphosites normalized by protein length") +
    theme_paper

phosphosites_known_isoforms <- proteoform_fdata_gene_level  %>% 
    filter(!is.na(phosphosite_count)) %>% 
    mutate(known_isoforms_binned = cut(known_isoforms, breaks = c(0, 1, 2, 4, 6, 32))) %>% 
    filter(!is.na(known_isoforms_binned)) %>% 
    ggplot(aes(as.factor(known_isoforms_binned), log2(phosphosite_count))) +
    geom_violin(fill = "gray", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.15, outlier.colour = NA) +
    geom_hline(aes(yintercept = median(log2(phosphosite_count))), 
               color = "black", alpha = 0.25, size = 2) +
    labs(x = "Binned number of known isoforms",
         y = "log2 number of phosphosites") +
    #ggtitle("Number of phosphosites") +
    theme_paper

# correct for length
norm_phosphosites_known_isoforms <- proteoform_fdata_gene_level  %>% 
    filter(!is.na(phosphosite_count)) %>% 
    mutate(known_isoforms_binned = cut(known_isoforms, breaks = c(0, 1, 2, 4, 6, 32))) %>% 
    filter(!is.na(known_isoforms_binned)) %>% 
    ggplot(aes(as.factor(known_isoforms_binned), log2(phosphosite_count/length))) +
    geom_violin(fill = "gray", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.15, outlier.colour = NA) +
    geom_hline(aes(yintercept = median(log2(phosphosite_count/length))), 
               color = "black", alpha = 0.25, size = 2) +
    labs(x = "Binned number of known isoforms",
         y = "log2 number of phosphosites\nnormalized by protein length") +
    #ggtitle("Number of phosphosites normalized by protein length") +
    theme_paper

# fraction disorder
disorder_detected_proteoforms <- proteoform_fdata_gene_level  %>% 
    filter(!is.na(fraction_disorder)) %>% 
    ggplot(aes(as.factor(num_proteoforms), fraction_disorder)) +
    geom_violin(fill = "gray", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.15, outlier.colour = NA) +
    geom_hline(aes(yintercept = median(fraction_disorder)), 
               color = "black", alpha = 0.25, size = 2) +
    labs(x = "Number of detected proteoforms",
         y = "Fraction of disorder\nin protein sequence") +
    theme_paper

# fraction disorder known isoforms
disorder_known_isoforms <- proteoform_fdata_gene_level  %>% 
    filter(!is.na(fraction_disorder)) %>% 
    mutate(known_isoforms_binned = cut(known_isoforms, breaks = c(0, 1, 2, 4, 6, 32))) %>% 
    filter(!is.na(known_isoforms_binned)) %>% 
    ggplot(aes(as.factor(known_isoforms_binned), fraction_disorder)) +
    geom_violin(fill = "gray", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.15, outlier.colour = NA) +
    geom_hline(aes(yintercept = median(fraction_disorder)), 
               color = "black", alpha = 0.25, size = 2) +
    labs(x = "Binned number of known isoforms",
         y = "Fraction of disorder\nin protein sequence") +
    #ggtitle("Number of phosphosites normalized by protein length") +
    theme_paper

# abundance
abundance_detected_proteoforms <- proteoform_fdata_gene_level %>% 
    dplyr::select(ioi, num_proteoforms, matches("abundance")) %>% 
    gather(key, value, matches("abundance")) %>% 
    group_by(ioi, num_proteoforms) %>% 
    summarize(mean_value = mean(value, na.rm = TRUE),
              .groups = "drop") %>% 
    filter(is.finite(mean_value)) %>% 
    ggplot(aes(as.factor(num_proteoforms), mean_value)) +
    geom_violin(fill = "gray", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.15, outlier.colour = NA) +
    geom_hline(aes(yintercept = median(mean_value)), 
               color = "black", alpha = 0.25, size = 2) +
    labs(x = "Number of detected proteoforms",
         y = "log2 average protein intensity\nat 41 degree C") +
    #ggtitle("Average protein intensity of proteins with different numbers of detected proteoforms") +
    theme_paper

abundance_known_isoforms <- proteoform_fdata_gene_level %>% 
    mutate(known_isoforms_binned = cut(known_isoforms, breaks = c(0, 1, 2, 4, 6, 32))) %>% 
    dplyr::select(ioi, known_isoforms_binned, matches("abundance")) %>% 
    filter(!is.na(known_isoforms_binned)) %>%
    gather(key, value, matches("abundance")) %>% 
    group_by(ioi, known_isoforms_binned) %>% 
    summarize(mean_value = mean(value, na.rm = TRUE),
              .groups = "drop") %>% 
    filter(is.finite(mean_value)) %>% 
    ggplot(aes(as.factor(known_isoforms_binned), mean_value)) +
    geom_violin(fill = "gray", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.15, outlier.colour = NA) +
    geom_hline(aes(yintercept = median(mean_value)), 
               color = "black", alpha = 0.25, size = 2) +
    labs(x = "Binned number of known isoforms",
         y = "log2 average protein intensity\nat 41 degree C") +
    #ggtitle("Average protein intensity of proteins with different numbers of detected proteoforms") +
    theme_paper

# assemble figure
plot_grid(
    plot_grid(detected_vs_known_proteoforms, cc_loc_detected_proteoforms, cc_loc_known_isoforms,
              nrow = 1, labels = letters[1:3]),
    plot_grid(length_detected_proteoforms, length_known_isoforms, 
              abundance_detected_proteoforms, abundance_known_isoforms,
              nrow = 1, labels = letters[4:7]),
    plot_grid(identified_peptides_detected_proteoforms, identified_peptides_known_isoforms, 
              norm_peptides_detected_proteoforms, norm_peptides_known_isoforms,
              nrow = 1, labels = letters[8:11]),
    plot_grid(phosphosites_detected_proteoforms, phosphosites_known_isoforms, 
              norm_phosphosites_detected_proteoform, norm_phosphosites_known_isoforms,
              nrow = 1, labels = letters[12:15]),
    plot_grid(half_life_detected_proteoforms, half_life_known_isoforms, 
              disorder_detected_proteoforms, disorder_known_isoforms,
              nrow = 1, labels = letters[16:19]),
    nrow = 5, rel_heights = c(1, 0.8, 0.8, 0.8, 0.8)
)
ggsave(filename = here("R/figures/suppl_fig_proteoform_features.pdf"), 
       width = 21, height = 31, units = "cm")
