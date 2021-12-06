library(here)
library(tidyverse)
library(cowplot)

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

sourceDir <- function(path, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        source(file.path(path, nm), ...)
    }
}
sourceDir(path = file.path(here("R/drugsens_cor")))

# define ppi coaggregation analysis output folder
drugsens_cor_folder <- here("drugsens_cor/output")

# define output folder for figures
figure_output_folder <- here("R/figures")

# load nparc hq results
nparc_res_hq_df <- readRDS(here("nparc/output/standard/nparc_res_hq_df.RDS"))

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

# get ABL1 fusion subtypes
abl_fusion_anno <- 
    sample_meta_raw %>% filter(grepl("ABL", subtype))

cl_colors <-  c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
                "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", 
                "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", 
                "#771122", "#AA4455")
# Peter Carl - R-Bloggers 02/2013 "The Paul Tol 21-color salute"
names(cl_colors) <- sort(
    unique(gsub("_", "-", sub("_BR2", "", proteoform_df$sample_name_machine))))

# read in preprocessed drugsens data
dss_mat <- readRDS(here("drugsens_cor/output/dss_mat.RDS"))

# read in preprocessed AUC matrix
auc_mat_norm <- readRDS(here("drugsens_cor/output/auc_mat_norm.RDS"))

# create full hq auc df
auc_full_hq_df <- nparc_res_hq_df %>% 
    dplyr::select(sample = sample_name, id, aumc) 

# convert to matrix
auc_full_hq_mat <- auc_full_hq_df %>% 
    dplyr::select(id, aumc, sample) %>% 
    spread(sample, aumc) %>% 
    column_to_rownames("id") %>% 
    as.matrix()

# normalize AUC values
auc_full_hq_mat_norm <- preprocessCore::normalize.quantiles(auc_full_hq_mat)
dimnames(auc_full_hq_mat_norm) <- dimnames(auc_full_hq_mat)

# define axis lavels
x_label <- expression("Temperature"* " " * "("*degree*C*")")
y_label <- "Fraction non-denatured"

# get fignl specific proteoform dataset
fignl_proteoform_df <- proteoform_df %>% 
    filter(gene == "FIGNL1_1", !grepl("BR", sample_name_machine)) %>% 
    left_join(nparc_res_hq_df %>% dplyr::select(id, sample_name, conv),
              by = c("gene" = "id", "sample_name_machine" = "sample_name")) %>% 
    filter(conv == TRUE)  %>% 
    mutate(sample_name_machine = gsub("_", "-", sample_name_machine))

# plot FIGNL1_1 proteoforms across cell lines
fignl_proteoform_profile_plot <- 
    ggplot(fignl_proteoform_df, 
           aes(temperature, rel_value, color = sample_name_machine)) +
    geom_smooth(aes(group = sample_name_machine), method = "lm",
                formula = 'y ~ splines::ns(x, df = 4)',
                se = FALSE, alpha = 0.5, size = 0.75) +
    # stat_smooth(aes(group = gene, linetype = gene),
    #             method = "nls", se = FALSE,
    #             formula = y ~ (1-a)/(1 + exp(-(b/x - c))) + a,
    #             method.args = list(start = c(a = 0.1, b = 1550, c = 40),
    #                                algorithm = 'port'),
    #             geom = "line",
    #             size = 1) +
    scale_color_manual("Subtype", 
                       values = cl_colors) +
    theme_paper +
    theme(legend.position = "none") +
    labs(x = x_label,
         y = y_label) +
    ggtitle("FIGNL1 proteoform 1")

ggsave(filename = here("R/figures/fig_fignl1_spline_fit.pdf"), 
       width = 7, height = 7, units = "cm")

# get fignl specific proteoform dataset
fignl2_proteoform_df <- proteoform_df %>% 
    filter(gene == "FIGNL1_2", !grepl("BR", sample_name_machine)) %>% 
    left_join(nparc_res_hq_df %>% dplyr::select(id, sample_name, conv),
              by = c("gene" = "id", "sample_name_machine" = "sample_name")) %>% 
    filter(conv == TRUE)  %>% 
    mutate(sample_name_machine = gsub("_", "-", sample_name_machine))

# plot FIGNL1_2 proteoforms across cell lines
fignl2_proteoform_profile_plot <- 
    ggplot(fignl2_proteoform_df, 
           aes(temperature, rel_value, color = sample_name_machine)) +
    geom_smooth(aes(group = sample_name_machine), method = "lm",
                formula = 'y ~ splines::ns(x, df = 4)',
                se = FALSE, alpha = 0.5, size = 0.75) +
    # stat_smooth(aes(group = gene, linetype = gene),
    #             method = "nls", se = FALSE,
    #             formula = y ~ (1-a)/(1 + exp(-(b/x - c))) + a,
    #             method.args = list(start = c(a = 0.1, b = 1550, c = 40),
    #                                algorithm = 'port'),
    #             geom = "line",
    #             size = 1) +
    scale_color_manual("Subtype", 
                       values = cl_colors) +
    theme_paper +
    theme(legend.position = "none") +
    labs(x = x_label,
         y = y_label) +
    ggtitle("FIGNL1 proteoform 2")

# read in results from limma analysis
limma_out_df <- readRDS(file.path(drugsens_cor_folder, "limma_out_df.RDS"))

limma_out_anno_df <- limma_out_df %>% 
    mutate(id = 1:n()) 

# loop over significant hits and make scatter plots of thermal stability and drug sens.
pList <- lapply(seq(nrow(filter(limma_out_anno_df, p_adj < 0.1))), function(i) {
    protn <- limma_out_anno_df$testedProt[i] 
    drugn <- limma_out_anno_df$rowname[i]
    p_adj <- limma_out_anno_df$p_adj[i]
    
    title_string <- paste(drugn, protn, sep = " vs. ")
    
    auc_fil_df <- auc_mat_norm %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "id") %>% 
        filter(id == protn) %>% 
        gather(sample, aumc, -id) %>% 
        filter(!grepl("_BR", sample)) 
    
    
    dss_fil_df <- dss_mat %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "sample") %>% 
        gather(drug, dss, -sample) %>% 
        filter(drug == drugn) 
    
    plot_df <- left_join(
        auc_fil_df,
        dss_fil_df,
        by = "sample") 
    
    ggplot(plot_df, aes(aumc, dss)) +
        geom_point(aes(color = sample)) +
        geom_smooth(method = "lm", color = "black") +
        scale_color_manual(values = cl_colors) +
        labs(x = "Area under the melting curve",
             y = "Drug sensitivity") +
        ggpubr::stat_cor(method = "pearson",
                         label.x.npc = "left",
                         cor.coef.name = "rho") +
        # geom_text(aes(x, y), label = round(p_adj, 3), 
        #           data = tibble(x = 0.9 * max(plot_df$aumc, na.rm = TRUE),
        #                         y = 0.9 * max(plot_df$dss, na.rm = TRUE))) +
        ggtitle(title_string) +
        theme_paper +
        theme(legend.position = "none")
})

# make combo plot of FIGNL1_1 profile and drug associations
plot_grid(fignl_proteoform_profile_plot, pList[[1]], pList[[3]],
          pList[[4]], pList[[6]], pList[[7]],
          ncol = 3, nrow = 2)

ggsave(filename = here("R/figures/suppl_fig_drugsens_cor_fignl.pdf"), 
       width = 21, height = 14, units = "cm")

# read in qMS data and make abundance boxplot for FBP1
qms_eset <- readRDS(here("data/proteins.B.RDS"))

# make proteoform data frame
qms_df <- biobroom::tidy.ExpressionSet(
    qms_eset) 

# get qms sample annotation
qms_pdata_df <- pData(qms_eset) %>% 
    as_tibble() %>% 
    dplyr::select(sample = proteomics_id,
                  sample_name = Cell_Line_Name_Paper) %>% 
    mutate(sample_name_machine = sub("h", "", sub("_LL", "", gsub("-", "_", sample_name)))) %>% 
    mutate(cell_line = sub("_BR2", "", sample_name_machine))

# join anno with qms data
qms_anno_df <- left_join(qms_df, qms_pdata_df, by = "sample")


# FIGNL1_1 thermal stability abundance comparison
fignl1_1_auc_df <- auc_df %>% 
    filter(id == "FIGNL1_1") %>% 
    gather(key, value, -id) %>% 
    separate(key, c("sample", "remove"), sep = "_AUM") %>% 
    dplyr::select(-remove)

fignl1_qms_fignl1_1_auc_df <- qms_anno_df %>% 
    filter(gene == "FIGNL1") %>% 
    dplyr::select(gene, value, sample = sample_name_machine) %>% 
    left_join(fignl1_1_auc_df, by = "sample") %>% 
    na.omit() %>% 
    mutate(sample = gsub("_", "-", sample))

fignl1_qms_fignl1_1_auc_scatter <- 
    ggplot(fignl1_qms_fignl1_1_auc_df, aes(value.x, value.y)) +
    geom_smooth(method = "lm", color = "black") +
    geom_point(aes(color = sample)) +
    ggpubr::stat_cor(method = "pearson") +
    scale_color_manual(values = cl_colors) +
    labs(x = bquote('log'[2]*' relative protein abundance'),
         y = "Drug sensitivity") +
    ggtitle("FIGNL1 abundance vs. FIGNL1_1 thermal stability") +
    theme_paper +
    theme(legend.position = "none")

# FIGNL1 qMS Eribulin scatter
fignl1_qms_dss_df <- qms_anno_df %>% 
    filter(gene == "FIGNL1") %>% 
    dplyr::select(gene, value, sample = sample_name_machine) %>% 
    left_join(dss_df %>% filter(drug_name == "Eribulin"),
              by = "sample") %>% 
    na.omit() %>% 
    mutate(sample = gsub("_", "-", sample))

fignl1_qms_dss_scatter <- 
    ggplot(fignl1_qms_dss_df, aes(value, dss)) +
    geom_smooth(method = "lm", color = "black") +
    geom_point(aes(color = sample)) +
    ggpubr::stat_cor(method = "pearson") +
    scale_color_manual(values = cl_colors) +
    labs(x = bquote('log'[2]*' relative protein abundance'),
         y = "Drug sensitivity") +
    ggtitle("FIGNL1 abundance vs. Erubilin") +
    theme_paper +
    theme(legend.position = "none")

# FIGNL1 qMS Vinorelbine scatter
fignl1_qms_vino_dss_df <- qms_anno_df %>% 
    filter(gene == "FIGNL1") %>% 
    dplyr::select(gene, value, sample = sample_name_machine) %>% 
    left_join(dss_df %>% filter(drug_name == "Vinorelbine"),
              by = "sample") %>% 
    na.omit() %>% 
    mutate(sample = gsub("_", "-", sample))

fignl1_qms_vino_dss_scatter <- 
    ggplot(fignl1_qms_vino_dss_df, aes(value, dss)) +
    geom_smooth(method = "lm", color = "black") +
    geom_point(aes(color = sample)) +
    ggpubr::stat_cor(method = "pearson") +
    scale_color_manual(values = cl_colors) +
    labs(x = bquote('log'[2]*' relative protein abundance'),
         y = "Drug sensitivity") +
    ggtitle("FIGNL1 abundance vs. Vinorelbine") +
    theme_paper +
    theme(legend.position = "none")

fignl2_eribulin_scatter <- plot_auc_dss_per_protein_drug_scatter(
    auc_mat = auc_full_hq_mat_norm, dss_mat = dss_mat, protein_id = "FIGNL1_2", 
    drug_name = "Eribulin", sample_anno = sample_meta_raw) + 
    scale_color_manual(values = cl_colors) +
    ggtitle("Eribulin vs. FIGNL1_2") +
    labs(x = "Area under the melting curve",
         y = "Drug sensitivity") +
    theme_paper +
    theme(legend.position = "none") 

fignl2_vino_scatter <- plot_auc_dss_per_protein_drug_scatter(
    auc_mat = auc_full_hq_mat_norm, dss_mat = dss_mat, protein_id = "FIGNL1_2", 
    drug_name = "Vinorelbine", sample_anno = sample_meta_raw) + 
    scale_color_manual(values = cl_colors) +
    ggtitle("Vinorelbine vs. FIGNL1_2") +
    labs(x = "Area under the melting curve",
         y = "Drug sensitivity") +
    theme_paper +
    theme(legend.position = "none") 

# make combo plot of FIGNL1_1 profile and drug associations
plot_grid(fignl_proteoform_profile_plot, pList[[1]], pList[[3]],
          pList[[4]], pList[[6]], pList[[7]],
          fignl1_qms_fignl1_1_auc_scatter, 
          fignl1_qms_dss_scatter,
          fignl1_qms_vino_dss_scatter,
          fignl2_proteoform_profile_plot, 
          fignl2_eribulin_scatter, 
          fignl2_vino_scatter,
          labels = letters[1:9],
          ncol = 3, nrow = 4)

ggsave(filename = here("R/figures/suppl_fig_drugsens_cor_fignl_extended.pdf"), 
       width = 21, height = 28, units = "cm")

# CRKL figures
# get crkl_1 specific proteoform dataset
crkl_1_proteoform_df <- proteoform_df %>% 
    filter(gene == "CRKL_1", !grepl("BR", sample_name_machine)) %>% 
    left_join(nparc_res_hq_df %>% dplyr::select(id, sample_name, conv),
              by = c("gene" = "id", "sample_name_machine" = "sample_name")) %>% 
    filter(conv == TRUE)  %>% 
    mutate(sample_name_machine = gsub("_", "-", sample_name_machine)) %>% 
    mutate(selected_subtype = 
               case_when(subtype %in% abl_fusion_anno$subtype ~ subtype,
                         TRUE ~ "other")) %>% 
    mutate(selected_subtype = gsub("\\.", "-", selected_subtype)) %>% 
    mutate(selected_subtype = 
               factor(selected_subtype, 
                      levels = c("BCR-ABL1", "ZMIZ1-ABL1", "other")))

# plot CRKL_1 proteoforms across cell lines
crkl_1_proteoform_profile_plot <- 
    ggplot(crkl_1_proteoform_df, 
           aes(temperature, rel_value, color = selected_subtype)) +
    # geom_smooth(aes(group = sample_name_machine), method = "lm",
    #             formula = 'y ~ splines::ns(x, df = 4)',
    #             se = FALSE, alpha = 0.5, size = 0.75) +
    stat_smooth(aes(group = sample_name_machine),
                method = "nls", se = FALSE,
                formula = y ~ (1-a)/(1 + exp(-(b/x - c))) + a,
                method.args = list(start = c(a = 0.1, b = 1550, c = 40),
                                   algorithm = 'port'),
                geom = "line",
                size = 1) +
    scale_color_manual("Subtype", 
                       values = c("BCR-ABL1" = "red",
                                  "ZMIZ1-ABL1" = "orange",
                                  "other" = "gray")) +
    theme_paper +
    theme(legend.position = "bottom") +
    labs(x = x_label,
         y = y_label) +
    ggtitle("CRKL proteoform 1")

# get crkl_2 specific proteoform dataset
crkl_2_proteoform_df <- proteoform_df %>% 
    filter(gene == "CRKL_2", !grepl("BR", sample_name_machine)) %>% 
    left_join(nparc_res_hq_df %>% dplyr::select(id, sample_name, conv),
              by = c("gene" = "id", "sample_name_machine" = "sample_name")) %>% 
    filter(conv == TRUE)  %>% 
    mutate(sample_name_machine = gsub("_", "-", sample_name_machine)) %>% 
    mutate(selected_subtype = 
               case_when(subtype %in% abl_fusion_anno$subtype ~ subtype,
                         TRUE ~ "other")) %>% 
    mutate(selected_subtype = gsub("\\.", "-", selected_subtype)) %>% 
    mutate(selected_subtype = 
               factor(selected_subtype, 
                      levels = c("BCR-ABL1", "ZMIZ1-ABL1", "other")))

# plot CRKL_2 proteoforms across cell lines
crkl_2_proteoform_profile_plot <- 
    ggplot(crkl_2_proteoform_df, 
           aes(temperature, rel_value, color = selected_subtype)) +
    geom_smooth(aes(group = sample_name_machine), method = "lm",
                formula = 'y ~ splines::ns(x, df = 4)',
                se = FALSE, alpha = 0.5, size = 0.75) +
    stat_smooth(aes(group = sample_name_machine),
                method = "nls", se = FALSE,
                formula = y ~ (1-a)/(1 + exp(-(b/x - c))) + a,
                method.args = list(start = c(a = 0.1, b = 1550, c = 40),
                                   algorithm = 'port'),
                geom = "line",
                size = 1) +
    scale_color_manual("Subtype", 
                       values = c("BCR-ABL1" = "red",
                                  "ZMIZ1-ABL1" = "orange",
                                  "other" = "gray")) +
    theme_paper +
    theme(legend.position = "bottom") +
    labs(x = x_label,
         y = y_label) +
    ggtitle("CRKL proteoform 2")

# combined CRKL figure
crkl_1_2_profile <- proteoform_df %>% 
    filter(gene %in% c("CRKL_1", "CRKL_2"), !grepl("BR", sample_name_machine)) %>% 
    left_join(nparc_res_hq_df %>% dplyr::select(id, sample_name, conv),
              by = c("gene" = "id", "sample_name_machine" = "sample_name")) %>% 
    filter(conv == TRUE)  %>% 
    bind_rows(bind_rows(lapply(unique(filter(proteoform_df, !grepl("BR", sample_name_machine))$sample_name_machine), function(samp){
        proteoform_df %>% 
            filter(gene %in% c("CRKL_1", "CRKL_2"), !grepl("BR", sample_name_machine)) %>% 
            group_by(temperature) %>% 
            summarise(rel_value = mean(rel_value), .groups = "keep") %>% 
            ungroup() %>% 
            mutate(gene = "mean", sample_name_machine = samp, subtype = "mean")
    }))) %>% 
    mutate(sample_name_machine = gsub("_", "-", sample_name_machine)) %>% 
    mutate(selected_subtype = 
               case_when(subtype %in% abl_fusion_anno$subtype ~ subtype,
                         subtype == "mean" ~ "mean",
                         TRUE ~ "other")) %>% 
    mutate(selected_subtype = gsub("\\.", "-", selected_subtype)) %>% 
    mutate(selected_subtype = 
               factor(selected_subtype, 
                      levels = c("BCR-ABL1", "ZMIZ1-ABL1", "other", "mean"))) %>% 
    ggplot(aes(temperature, rel_value, color = selected_subtype)) + 
    stat_smooth(aes(group = gene, linetype = gene),
                method = "nls", se = FALSE,
                formula = y ~ (1-a)/(1 + exp(-(b/x - c))) + a,
                method.args = list(start = c(a = 0.1, b = 1550, c = 40),
                                   algorithm = 'port'),
                geom = "line",
                size = 1) +
    scale_color_manual("Subtype", 
                       values = c("BCR-ABL1" = "red",
                                  "ZMIZ1-ABL1" = "orange",
                                  "other" = "gray",
                                  "mean" = "black")) +
    labs(x = x_label,
         y = y_label) +
    ggtitle("CRKL proteoforms 1 & 2") +
    facet_wrap(~sample_name_machine) +
    theme_paper

# CRKL qMS Imatinib scatter
crkl_qms_dss_df <- qms_anno_df %>% 
    filter(gene == "CRKL") %>% 
    dplyr::select(gene, value, sample = sample_name_machine) %>% 
    left_join(dss_df %>% filter(drug_name == "Imatinib"),
              by = "sample") %>% 
    na.omit() %>% 
    mutate(sample = gsub("_", "-", sample))

crkl_qms_dss_scatter <- 
    ggplot(eps8l2_qms_dss_df, aes(value, dss)) +
    geom_smooth(method = "lm", color = "black") +
    geom_point(aes(color = sample)) +
    ggpubr::stat_cor(method = "pearson") +
    scale_color_manual(values = cl_colors) +
    labs(x = bquote('log'[2]*' relative protein abundance'),
         y = "Drug sensitivity") +
    ggtitle("CRKL abundance vs. Imatinib") +
    theme_paper +
    theme(legend.position = "none") 


crkl_2_imatinib_scatter <- plot_auc_dss_per_protein_drug_scatter(
    auc_mat = auc_full_hq_mat_norm, dss_mat = dss_mat, protein_id = "CRKL_2", 
    drug_name = "Imatinib", sample_anno = sample_meta_raw) + 
    scale_color_manual(values = cl_colors) +
    ggtitle("Imatinib vs. CRKL_2") +
    labs(x = "Area under the melting curve",
         y = "Drug sensitivity") +
    theme_paper +
    theme(legend.position = "none") 

# assemble crkl suppl. figure
plot_grid(
    crkl_1_2_profile,
    plot_grid(crkl_qms_dss_scatter,
              crkl_2_imatinib_scatter,
              ncol = 1, nrow = 2),
    rel_widths = c(2, 1),
    ncol = 2, nrow = 1)

ggsave(filename = here("R/figures/suppl_fig_drugsens_cor_crkl_extended.pdf"), 
       width = 21, height = 14, units = "cm")

# EPS8L2 qMS Eltanexor scatter
eps8l2_qms_dss_df <- qms_anno_df %>% 
    filter(gene == "EPS8L2") %>% 
    dplyr::select(gene, value, sample = sample_name_machine) %>% 
    left_join(dss_df %>% filter(drug_name == "Eltanexor"),
              by = "sample") %>% 
    na.omit() %>% 
    mutate(sample = gsub("_", "-", sample))

eps8l2_qms_dss_scatter <- 
    ggplot(eps8l2_qms_dss_df, aes(value, dss)) +
    geom_smooth(method = "lm", color = "black") +
    geom_point(aes(color = sample)) +
    ggpubr::stat_cor(method = "pearson") +
    scale_color_manual(values = cl_colors) +
    labs(x = bquote('log'[2]*' relative protein abundance'),
         y = "Drug sensitivity") +
    ggtitle("EPS8L2 abundance vs. Eltanexor") +
    theme_paper +
    theme(legend.position = "none") 

# get EPS8L2 specific proteoform dataset
eps8l2_proteoform_df <- proteoform_df %>% 
    filter(grepl("EPS8L2_", gene), !grepl("BR", sample_name_machine)) %>% 
    left_join(nparc_res_hq_df %>% dplyr::select(id, sample_name, conv),
              by = c("gene" = "id", "sample_name_machine" = "sample_name")) %>% 
    filter(conv == TRUE) %>% 
    mutate(sample_name_machine = gsub("_", "-", sample_name_machine))

# plot EPS8L2 proteoforms across cell lines
eps8l2_proteoform_profile_plot <- 
    ggplot(eps8l2_proteoform_df, 
           aes(temperature, rel_value, color = sample_name_machine)) +
    #geom_point() +
    geom_smooth(aes(group = gene, linetype = gene), method = "lm",
                formula = 'y ~ splines::ns(x, df = 4)',
                se = FALSE, alpha = 0.5, size = 0.75) +
    # stat_smooth(aes(group = sample_name_machine),
    #             method = "nls", se = FALSE,
    #             formula = y ~ (1-a)/(1 + exp(-(b/x - c))) + a,
    #             method.args = list(start = c(a = 0, b = 550, c = 10),
    #                                algorithm = 'port'),
    #             geom = "line",
    #             alpha = 0.5) +
    scale_color_manual("Cell line", 
                       values = cl_colors) +
    facet_wrap(~sample_name_machine) +
    theme_paper +
    theme(legend.position = "bottom") +
    labs(x = x_label,
         y = y_label) +
    ggtitle("EPS8L2 proteoforms 1 and 2") 


eps8l2_1_eltanexor_scatter <- plot_auc_dss_per_protein_drug_scatter(
    auc_mat = auc_mat_norm, dss_mat = dss_mat, protein_id = "EPS8L2_1", 
    drug_name = "Eltanexor", sample_anno = sample_meta_raw) + 
    scale_color_manual(values = cl_colors) +
    theme_paper +
    theme(legend.position = "none") 


plot_grid(eps8l2_proteoform_profile_plot,
          plot_grid(
              eps8l2_qms_dss_scatter,
              eps8l2_1_eltanexor_scatter,
              labels = letters[2:3],
              ncol = 1
          ), rel_widths = c(2, 1),
          labels = c("a", NA))

ggsave(filename = here("R/figures/suppl_fig_drugsens_cor_eps8l2_combo.pdf"), 
       width = 21, height = 14, units = "cm")
