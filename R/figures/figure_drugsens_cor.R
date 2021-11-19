library(tidyverse)
library(Biobase)
library(NPARC)
library(limma)
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

sourceDir <- function(path, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        source(file.path(path, nm), ...)
    }
}
sourceDir(path = file.path(here("R/drugsens_cor")))

# get sample meta data
sample_meta_raw <- read_tsv(here("meta/sample_meta.txt"))

# define drugsens analysis output folder
drugsens_cor_folder <- here("drugsens_cor/output")

# define output folder for figures
figure_output_folder <- here("R/figures")

# load proteoform data
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

# get ABL1 fusion subtypes
abl_fusion_anno <- 
    sample_meta_raw %>% filter(grepl("ABL", subtype))

abl1_proteoform_df <- proteoform_df %>% 
    filter(gene == "ABL1_1") %>% 
    mutate(selected_subtype = 
               case_when(subtype %in% abl_fusion_anno$subtype ~ subtype,
                         TRUE ~ "other")) %>% 
    mutate(selected_subtype = gsub("\\.", "-", selected_subtype)) %>% 
    mutate(selected_subtype = 
               factor(selected_subtype, 
                      levels = c("BCR-ABL1", "ZMIZ1-ABL1", "other")))

# define axis lavels
x_label <- expression("Temperature"* " " * "("*degree*C*")")
y_label <- "Fraction non-denatured"

# plot ABL1_1 proteoforms highlighting cell lines with fusion variants
abl1_proteoform_profile_plot <- 
    ggplot(abl1_proteoform_df, 
           aes(temperature, rel_value, color = selected_subtype)) +
    # geom_smooth(aes(group = sample_name_machine), method = "lm",
    #             formula = 'y ~ splines::ns(x, df = 4)',
    #             se = FALSE, alpha = 0.5) +
    stat_smooth(aes(group = sample_name_machine),
                method = "nls", se = FALSE,
                formula = y ~ (1-a)/(1 + exp(-(b/x - c))) + a,
                method.args = list(start = c(a = 0, b = 550, c = 10),
                                   algorithm = 'port'),
                geom = "line",
                alpha = 0.5) +
    scale_color_manual("Subtype", 
                       values = c("BCR-ABL1" = "red",
                                  "ZMIZ1-ABL1" = "orange",
                                  "other" = "gray")) +
    theme_paper +
    theme(legend.position = "bottom") +
    labs(x = x_label,
         y = y_label) +
    ggtitle("ABL1 proteoform 1") 

# create data frame to plot abl1 thermal stability against imatinib sensitivity
abl1_imatinib_scatter_df <- 
    plot_auc_dss_per_protein_drug_scatter(
        auc_mat = nparc_res_df %>% 
            dplyr::select(sample = sample_name, id, aumc) %>% 
            dplyr::select(id, aumc, sample) %>% 
            spread(sample, aumc) %>% 
            column_to_rownames("id") %>% 
            as.matrix(),
        dss_mat = dss_mat, 
        sample_anno = sample_meta_raw,
        protein_id = "ABL1_1", 
        drug_name = "Imatinib")$data %>% 
    mutate(selected_subtype = 
               case_when(subtype %in% abl_fusion_anno$subtype ~ subtype,
                         TRUE ~ "other")) %>% 
    mutate(selected_subtype = gsub("\\.", "-", selected_subtype)) %>% 
    mutate(selected_subtype = 
               factor(selected_subtype, 
                      levels = c("BCR-ABL1", "ZMIZ1-ABL1", "other")))

# plot abl1 thermal stability against imatinib sensitivity
abl1_imatinib_scatter <- 
    ggplot(abl1_imatinib_scatter_df, aes(aumc, dss)) +
    geom_point(aes(color = selected_subtype), alpha = 0.5) +
    geom_smooth(method = "lm") +
    scale_color_manual("", values = c("BCR-ABL1" = "red",
                                      "ZMIZ1-ABL1" = "orange",
                                      "other" = "gray")) +
    ggpubr::stat_cor(method = "pearson",
                     label.x.npc = "left",
                     cor.coef.name = "rho") +
    theme_paper + 
    theme(legend.position = "bottom") +
    labs(x = "Area under the melting curve of ABL1_1",
         y = "Imatinib drug sensitivity")

# create data frame to plot abl1 thermal stability against nilotinib sensitivity
abl1_nilotinib_scatter_df <- 
    plot_auc_dss_per_protein_drug_scatter(
        auc_mat = nparc_res_df %>% 
            dplyr::select(sample = sample_name, id, aumc) %>% 
            dplyr::select(id, aumc, sample) %>% 
            spread(sample, aumc) %>% 
            column_to_rownames("id") %>% 
            as.matrix(),
        dss_mat = dss_mat, 
        sample_anno = sample_meta_raw,
        protein_id = "ABL1_1", 
        drug_name = "Nilotinib")$data %>% 
    mutate(selected_subtype = 
               case_when(subtype %in% abl_fusion_anno$subtype ~ subtype,
                         TRUE ~ "other")) %>% 
    mutate(selected_subtype = gsub("\\.", "-", selected_subtype)) %>% 
    mutate(selected_subtype = 
               factor(selected_subtype, 
                      levels = c("BCR-ABL1", "ZMIZ1-ABL1", "other")))

# plot abl1 thermal stability against nilotinib sensitivity
abl1_nilotinib_scatter <- 
    ggplot(abl1_nilotinib_scatter_df, aes(aumc, dss)) +
    geom_point(aes(color = selected_subtype), alpha = 0.5) +
    geom_smooth(method = "lm") +
    scale_color_manual("", values = c("BCR-ABL1" = "red",
                                      "ZMIZ1-ABL1" = "orange",
                                      "other" = "gray")) +
    ggpubr::stat_cor(method = "pearson",
                     label.x.npc = "left",
                     cor.coef.name = "rho") +
    theme_paper + 
    theme(legend.position = "bottom") +
    labs(x = "Area under the melting curve of ABL1_1",
         y = "Nilotinib drug sensitivity")

abl1_legend <- get_legend(
    proteoform_profile_plot +
        geom_point(aes(color = selected_subtype))
)

abl1_plots <-  
    plot_grid(
        plot_grid(
            abl1_proteoform_profile_plot + theme(legend.position = "none"),
            abl1_imatinib_scatter + theme(legend.position = "none"),
            abl1_nilotinib_scatter + theme(legend.position = "none"),
            labels = letters[1:3], ncol = 3),
        abl1_legend, rel_heights = c(9, 1), ncol = 1)

abl1_plots

ggsave(filename = here("R/figures/fig_abl1_imatinib_nilotinib_cor.pdf"), 
       width = 21, height = 7, units = "cm")

# read in results from limma analysis
limma_out_df <- readRDS(file.path(drugsens_cor_folder, "limma_out_df.RDS"))

limma_out_anno_df <- limma_out_df %>% 
    mutate(id = 1:n()) 

# make volcano plot
limma_volcano <- ggplot(filter(limma_out_anno_df, p_adj >= 0.1), 
                        aes(logFC, -log10(P.Value))) +
    stat_binhex(aes(alpha = log10(..count..)), bins = 100, fill = "black") +
    geom_point(alpha = 0.5, color = "black", 
               data = filter(limma_out_anno_df, p_adj < 0.1)) +
    labs(x = bquote('log'[2]*'(fold change)'),
         y = bquote('-log'[10]*'('*italic(p)*'-value)')) +
    geom_point(alpha = 0.5, color = "green", 
               data = filter(limma_out_anno_df, p_adj < 0.1, grepl("FIGNL", testedProt))) +
    theme_paper +
    theme(legend.position = "none")

# loop over significant hits and make scatter plots of thermal stability and drug sens.
pList <- lapply(seq(nrow(filter(all_limma_out_df, p_adj < 0.1))), function(i) {
    protn <- all_limma_out_df$testedProt[i] 
    drugn <- all_limma_out_df$rowname[i]
    p_adj <- all_limma_out_df$p_adj[i]
    
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
        geom_text(aes(x, y), label = round(p_adj, 3), 
                  data = tibble(x = 0.9 * max(plot_df$aumc, na.rm = TRUE),
                                y = 0.9 * max(plot_df$dss, na.rm = TRUE))) +
        ggtitle(title_string) +
        theme_paper +
        theme(legend.position = "none")
})

# plot plots in a row
plot_grid(
    limma_volcano, pList[[8]], pList[[26]], ncol = 3)

ggsave(filename = here("R/figures/fig_drugsens_cor_volcano_plus_examples.pdf"), 
       width = 21, height = 7, units = "cm")