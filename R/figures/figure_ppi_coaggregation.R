library(here)
library(tidyverse)
library(DESeq2)
library(fgsea)
library(reactome.db)

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

# read in sample meta info
sample_meta_stages <- read_tsv(here("meta/sample_meta_stages.txt"))

# define ppi coaggregation analysis output folder
ppi_coaggregation_folder <- here("ppi_coaggregation/output")

# define output folder for figures
figure_output_folder <- here("R/figures")

# load proteoform ratios
proteoforms <- readRDS(
    here("proteoform_detection/output/standard/proteoforms_narrow_range_focused.RDS"))

# create proteoform data frame
proteoform_df <- biobroom::tidy.ExpressionSet(
    proteoforms, 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

cl_colors <-  c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
                "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", 
                "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", 
                "#771122", "#AA4455")
# Peter Carl - R-Bloggers 02/2013 "The Paul Tol 21-color salute"
names(cl_colors) <- sort(
    unique(gsub("_BR2", "", proteoform_df$sample_name_machine)))

# load individual cell line tpca results
tpca_result_list <- readRDS(
    file.path(ppi_coaggregation_folder,
              "tpca_result_list_narrow_range_focused_hq_filtered.RDS"))

# load multi-cell line tpca results
multi_cell_line_rtpca_df <- readRDS(
    file.path(ppi_coaggregation_folder,
              "multi_cell_line_rtpca_robust_df.RDS"))

# read in and preprocess drugsens data
dss_df <- read_delim(
    here("data/2020-09-09_sDSS_2D_matrix_20_with_annotation_Meltome.txt"), 
    delim = "\t") %>% 
    dplyr::select(-Mechanism.Targets, -Putative.Target.Protein,
                  -Class.explained, -High.phase.Approval.status,
                  -Res_code, -Alias, -Activity.modifier,
                  -Active_inactive.in.clinic, -Solvent,
                  -High.conc..nM., -InChI, -ChEMBL.ID) %>% 
    gather(sample, dss, -FIMM.ID, -DRUG.NAME) %>% 
    dplyr::select(id = FIMM.ID, drug_name = DRUG.NAME, sample, dss) %>% 
    mutate(sample = gsub("\\.", "_", sample)) %>% 
    mutate(sample = gsub("X_", "", sample)) %>% 
    mutate(sample = gsub("Kasumi", "KASUMI", sample)) %>% 
    # filter out COG_319 due to bad quality
    filter(sample != "COG_319")

# filter for minimal dss effect
dss_df <- filter_drugsens_for_minimal_effect(dss_df, min_effect = 6)

# volcano plot of multi-cell line tpca results
ppi_volcano <- ggplot(multi_cell_line_rtpca_df, aes(max_rss - min_rss, f_stat)) + 
    geom_point(alpha = 0.25, color = "gray") + 
    geom_point(color = "black", alpha = 0.5,
               data = filter(multi_cell_line_rtpca_df,
                             f_stat > quantile(multi_cell_line_rtpca_df$f_stat, 0.9))) + 
    geom_text_repel(
        aes(label = pair),  
        nudge_y = 200,
        direction = "y",
        segment.size = 0.25,
        color = "black", 
        data = filter(multi_cell_line_rtpca_df, pair %in% c("CXXC1_2:SETD1A_3",
                                                            "POLR3A_1:POLR3B_3",
                                                            "EXOC2_3:EXOC4_3",
                                                            "PSMA3_2:PSMA7_1",
                                                            "EIF3F_2:EIF3I_1",
                                                            "MDM2_2:TP53_1")))+ #, pair == "CXXC1_2:SETD1A_3")) +
    scale_x_log10() +
    coord_cartesian(xlim = c(0.005, 3)) +
    labs(x =  bquote('RSS'['(n-1)']~' - RSS'['(2)']~''),
         y = expression(''*italic(F)*'-statistic')) +
    theme_paper

ggsave(ppi_volcano, filename = file.path(figure_output_folder, "fig4_ppi_volcano.pdf"),
       width = 10, height = 8, units = "cm")


# plot profiles of CXXC1_2 and SETD1A_3 across cell lines
x_label <- expression("Temperature"* " " * "("*degree*C*")")
y_label <- "Fraction non-denatured"

cxxc1_setd1a_profiles <- proteoform_df %>% 
    filter(grepl("CXXC1_2|SETD1A_3", gene),
           !grepl("BR2", sample_name_machine)) %>%
    mutate(sample_name_machine = gsub("_", "-", sample_name_machine)) %>% 
    ggplot(aes(temperature, rel_value)) + 
    geom_point(aes(color = gene), alpha = 0.5) + 
    facet_wrap(~sample_name_machine) +
    scale_color_manual("Proteoform", values = c("darkorange", "steelblue")) +
    labs(x = x_label, y = y_label) +
    theme_paper +
    theme(legend.position = "bottom")

ggsave(cxxc1_setd1a_profiles, filename = 
           file.path(figure_output_folder, "fig4_cxxc1_setd1a_profiles.pdf"),
       width = 12, height = 12, units = "cm")


# RNA-seq analysis
# read in RNA-seq counts
rna_seq_df <- read_tsv(here("data/ALL_raw_counts_with_samples_info.txt")) %>% 
    column_to_rownames(var = "Gene.Name") %>% 
    dplyr::select(-matches("BR")) %>% 
    as.matrix()
colnames(rna_seq_df) <- sub("h$", "", sub("X_", "", gsub("\\.", "_", colnames(rna_seq_df))))

# define coldata
coldata <- data.frame(
    sample_name = colnames(rna_seq_df),
    condition = ifelse(colnames(rna_seq_df) %in% 
                           c("HAL_01", "MHH_CALL_4", "KASUMI_9", "KASUMI_2", 
                             "ALL_PO", "TMD5", "NALL_1", "COG_394"), 
                       "differential", "coaggregation")
) %>% 
    left_join(sample_meta_stages %>% 
                  mutate(sample = gsub("-", "_", sample_name)) %>% 
                  dplyr::select(sample_name=sample, sex=gender, stage=Assigned_stage),
              by = "sample_name")

rownames(coldata) <- colnames(rna_seq_df)

# perform DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = rna_seq_df,
                              colData = coldata,
                              design= ~sex + condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "differential", "coaggregation"))

# get DESeq2 results
plot_df <- res %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    as_tibble() %>% 
    arrange(padj)

# GSEA
# map gene symbols to entrez ids
entrez_map_df <- bitr(geneID = plot_df$rowname, fromType = "SYMBOL", 
                      toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db)

plot_df <- left_join(plot_df, entrez_map_df, by = c("rowname" = "SYMBOL"))
deseq_ranks <- filter(plot_df, !is.na(ENTREZID), !is.na(stat))$stat
names(deseq_ranks) <- filter(plot_df, !is.na(ENTREZID), !is.na(stat))$ENTREZID

# get reactome pathways
pathways <- reactomePathways(names(deseq_ranks))

# perform gene set enrichment
fgseaRes <- fgsea(pathways = pathways, 
                  stats    = deseq_ranks,
                  minSize  = 15,
                  maxSize  = 500)

# plot p53-indep. dna damage enrichment plot
plotEnrichment(pathways[["p53-Independent DNA Damage Response"]],
               deseq_ranks) + labs(title="p53-Independent DNA Damage Response") +
    theme_paper

ggsave(file.path(figure_output_folder, "fig4_dna_damage_enrichment.pdf"), 
       width = 8, height = 5, units = "cm")

# get Eucl. distances for CXXC1_2:SETD1A_3 in different cell lines
cxxc1_2_setd1a_3_dist_df <- bind_rows(lapply(names(tpca_result_list), function(nm){
    res_df <- tpca_result_list[[nm]]
    filter(res_df, complex_name == "CXXC1_2:SETD1A_3") %>% 
        mutate(sample = nm)
}))

# test for differences in drug sensitivity of nucleosid analagon drugs
cxxc1_2_setd1a_3_dist_dss_df <- cxxc1_2_setd1a_3_dist_df %>% 
    left_join(dss_df %>% filter(drug_name %in% c("Decitabine", "Azacitidine")), 
              by = "sample") %>% 
    na.omit()

ggplot(cxxc1_2_setd1a_3_dist_dss_df, aes(mean_dist < 0.1, dss)) + 
    geom_boxplot() + 
    ggsignif::geom_signif(comparisons = list(c("TRUE", "FALSE")),
                          test = t.test) + 
    geom_jitter(aes(color = sample), width = 0.05) +
    scale_color_manual("Cell line", values = cl_colors) +
    labs(x = "Eucledian distance CXXC1_2:SETD1A_3 < 0.1",
         y = "Drug sensitivity score") +
    facet_wrap(~drug_name) +
    theme_paper +
    theme(legend.position = "bottom")

ggsave(file.path(figure_output_folder, "fig4_cxxc1_2_setd1a_3_drug_sens.pdf"), 
       width = 8, height = 10, units = "cm")
