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

uniprot_length_df <- read_delim(
    here("data/uniprot-filtered-reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B96--.tab"),
    delim = "\t") %>% 
    dplyr::select(`Gene names`, Length) %>% 
    mutate(gene = sub(" .+", "", `Gene names`)) %>% 
    dplyr::select(gene, Length) %>% 
    group_by(gene) %>% 
    summarize(max_length = max(Length)) %>% 
    ungroup

num_prot_length_df <- left_join(number_of_proteoforms_df, uniprot_length_df, by = "gene")

ggplot(num_prot_length_df, aes(n, log2(max_length))) +
    geom_point(alpha = 0.2) +
    #ggpubr::stat_cor(cor.coef.name = "rho") +
    ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE) +
    theme_paper
    

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

# load graphs
graphs <- readRDS(file.path(proteoform_detection_folder, "graphs_comms.RDS"))
# filter graphs for 0 modularity
graphs_gt0 <- graphs[(lapply(graphs, get.graph.attribute, 
                             name = "proteoform_modularity") > 1e-13) %>% unlist()]

# convert graphs_gt0 to data frame
graph_gt0_df <- bind_rows(lapply(names(graphs_gt0), function(gr_nm){
    membership_vec <- graphs_gt0[[gr_nm]]$communities$membership
    tibble(gene = gr_nm,
           psms = names(membership_vec),
           proteoform = membership_vec,
           proteoform_id = paste(gr_nm, membership_vec, sep = "_")) %>% 
        arrange(proteoform)
})) %>% 
    filter(proteoform_id %in% proteoform_df$gene) %>% 
    mutate(group = "modularity > 0")

# define axis labels
x_label <- expression("Temperature"* " " * "("*degree*C*")")
y_label <- "Fraction non-denatured"

# TMPO peptides assigned to proteoforms plot
tmpo_proteofom_peptide_df <- graph_gt0_df %>% 
    filter(gene == "TMPO") %>% 
    left_join(peptides_df %>% 
                  filter(id == "TMPO") %>% 
                  na.omit(),
              by = c("psms" = "gene")) %>% 
    filter(!grepl("_BR", sample_name_machine))

ggplot(tmpo_proteofom_peptide_df, aes(temperature, value)) +
    geom_line(aes(group = psms, color = proteoform_id), alpha = 0.25) +
    scale_color_manual(values = c("#FF5376", "#72AFD9")) +
    facet_wrap(~sample_name_machine) +
    coord_cartesian(ylim = c(0, 1.75)) +
    labs(x = x_label, y = y_label) +
    theme_paper +
    theme(legend.position = "bottom")

ggsave(filename = here("R/figures/figure_tmpo_proteofom_peptide_df.pdf"), 
       width = 18, height = 9, units = "cm")

# akap1 peptides assigned to proteoforms plot
akap1_proteofom_peptide_df <- graph_gt0_df %>% 
    filter(gene == "AKAP1") %>% 
    left_join(peptides_df %>% 
                  filter(id == "AKAP1") %>% 
                  na.omit(),
              by = c("psms" = "gene")) %>% 
    filter(!grepl("_BR", sample_name_machine))

ggplot(akap1_proteofom_peptide_df, aes(temperature, value)) +
    geom_line(aes(group = psms, color = proteoform_id), alpha = 0.25) +
    scale_color_manual(values = c("#FF5376", "#72AFD9", "#E3D26F")) +
    facet_wrap(~sample_name_machine) +
    coord_cartesian(ylim = c(0, 1.75)) +
    labs(x = x_label, y = y_label) +
    ggtitle("AKAP1 peptides colored by detected proteoforms") +
    theme_paper +
    theme(legend.position = "bottom")

ggsave(filename = here("R/figures/figure_aka1_proteofom_peptide_df.pdf"), 
       width = 18, height = 9, units = "cm")

## run COPF on same examples to compare results
tmpo_akap1_copf_input_df <- bind_rows(tmpo_proteofom_peptide_df, akap1_proteofom_peptide_df) %>%
    mutate(intensity = value * 1000) %>% 
    group_by(psms) %>% 
    filter(n() == 160) %>% 
    ungroup() %>% 
    dplyr::select(protein_id = gene,
                  peptide_id = psms,
                  temperature,
                  sample_name_machine,
                  intensity) %>% 
    unite("filename", c("temperature", "sample_name_machine"))

write.csv(tmpo_akap1_copf_input_df, here("R/benchmark/tmpo_akap1_copf_input.csv"))

fraction_annotation <- tibble(filename = unique(tmpo_akap1_copf_input_df$filename))  %>%
    mutate(fraction_number = 1:n())

pepTraces <- importPCPdata(input_data = as.data.table(df),
                           fraction_annotation = as.data.table(fraction_annotation),
                           rm_decoys = FALSE)

traces_corr <- CCprofiler::calculateGeneCorrMatrices(pepTraces)

traces_clustered <- clusterPeptides(
    traces_corr,
    method = "average", plot = F, PDF=F,
    name="ProteoformClusters_deepmeltome_TMPO_AKAP1_")

traces_clusteredInN <- cutClustersInNreal(traces_clustered, clusterN = 2,
                                          min_peptides_per_cluster = 5)

write_csv(traces_clusteredInN$trace_annotation, here("R/benchmark/copf_peptide_clustering_tmpo_akap1.csv"))

copf_tmpo_akap_res <- read_csv(here("R/benchmark/copf_peptide_clustering_tmpo_akap1.csv"))


tmpo_proteofom_peptide_copf_df <- tmpo_proteofom_peptide_df %>% 
    left_join(copf_tmpo_akap_res, by = c("psms" = "id"))

ggplot(tmpo_proteofom_peptide_copf_df, aes(temperature, value)) +
    geom_line(aes(group = psms, color = as.factor(cluster)), alpha = 0.25) +
    scale_color_manual("COPF cluster", values = c("#FF5376", "#72AFD9", "#E3D26F")) +
    facet_wrap(~sample_name_machine) +
    coord_cartesian(ylim = c(0, 1.75)) +
    labs(x = x_label, y = y_label) +
    ggtitle("AKAP1 peptides colored by detected COPF cluster") +
    theme_paper +
    theme(legend.position = "bottom")

ggsave(filename = here("R/figures/figure_tmpo_copf_proteofom_peptide_df.pdf"), 
       width = 18, height = 9, units = "cm")

akap1_proteofom_peptide_copf_df <- akap1_proteofom_peptide_df %>% 
    left_join(copf_tmpo_akap_res, by = c("psms" = "id"))

ggplot(akap1_proteofom_peptide_copf_df, aes(temperature, value)) +
    geom_line(aes(group = psms, color = as.factor(cluster)), alpha = 0.25) +
    scale_color_manual(values = c("#FF5376", "#72AFD9", "#E3D26F")) +
    facet_wrap(~sample_name_machine) +
    coord_cartesian(ylim = c(0, 1.75)) +
    labs(x = x_label, y = y_label) +
    ggtitle("AKAP1 peptides colored by detected COPF cluster") +
    theme_paper +
    theme(legend.position = "bottom")


ggsave(filename = here("R/figures/figure_akap1_copf_proteofom_peptide_df.pdf"), 
       width = 18, height = 9, units = "cm")

# global peptides assigned to proteoforms plot
global_proteofom_peptide_df <- graph_gt0_df %>% 
    left_join(peptides_df %>% 
                  filter(temperature == 41, sample_name_machine == "697") %>% 
                  na.omit() %>% 
                  dplyr::select(gene, first_protein_id),
              by = c("psms" = "gene"))

global_proteoform_summary_df <- global_proteofom_peptide_df %>% 
    filter(!is.na(first_protein_id)) %>% 
    group_by(gene, first_protein_id) %>% 
    dplyr::summarize(proteoform_1_count = sum(proteoform == 1), 
                     proteoform_2_count = sum(proteoform == 2),
                     proteoform_3_count = sum(proteoform == 3),
                     proteoform_4_count = sum(proteoform == 4)) %>% 
    # ungroup() %>% 
    # group_by(gene) %>% 
    # filter(n() > 1) %>% 
    ungroup

global_proteoform_summary_df <- global_proteoform_summary_df %>% 
    group_by(gene) %>% 
    dplyr::summarise(known_proteoforms_reflected = 
                         (any(proteoform_1_count > proteoform_2_count) & 
                         any(proteoform_2_count > proteoform_1_count)) | 
                         (any(proteoform_1_count > proteoform_3_count) & 
                         any(proteoform_3_count > proteoform_1_count)) |
                         (any(proteoform_1_count > proteoform_4_count) & 
                         any(proteoform_4_count > proteoform_1_count)) | 
                         (any(proteoform_2_count > proteoform_4_count) & 
                        any(proteoform_4_count > proteoform_2_count)) | 
                         (any(proteoform_2_count > proteoform_3_count) & 
                        any(proteoform_3_count > proteoform_2_count)) |
                         (any(proteoform_3_count > proteoform_4_count) & 
                        any(proteoform_4_count > proteoform_3_count))) %>% 
    ungroup()

global_proteoform_summary_df
# # A tibble: 5,620 × 2
# gene  known_proteoforms_reflected
# <chr> <lgl>                      
#     1 A2M   FALSE                      
# 2 AAAS  FALSE                      
# 3 AACS  FALSE                      
# 4 AAGAB FALSE                      
# 5 AAK1  TRUE                       
# 6 AAR2  FALSE                      
# 7 AARS1 FALSE                      
# 8 AARS2 FALSE                      
# 9 AASDH FALSE                      
# 10 AASS  FALSE                      
# # … with 5,610 more rows

global_proteoform_summary_df %>% filter(known_proteoforms_reflected)
# A tibble: 1,307 × 2
# gene       known_proteoforms_reflected
# <chr>      <lgl>                      
#     1 AAK1       TRUE                       
# 2 ABCE1      TRUE                       
# 3 ABCF1      TRUE                       
# 4 ABL2       TRUE                       
# 5 ABLIM1     TRUE                       
# 6 ABR        TRUE                       
# 7 AC004151.1 TRUE                       
# 8 AC005747.1 TRUE                       
# 9 AC022966.1 TRUE                       
# 10 AC093012.1 TRUE                       
# # … with 1,297 more rows
