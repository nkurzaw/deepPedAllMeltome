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

POLR3A_1_POLR3B_3_profiles <- proteoform_df %>% 
    filter(grepl("POLR3A_1|POLR3B_3", gene),
           !grepl("BR2", sample_name_machine)) %>%
    mutate(sample_name_machine = gsub("_", "-", sample_name_machine)) %>% 
    ggplot(aes(temperature, rel_value)) + 
    geom_point(aes(color = gene), alpha = 0.5) + 
    facet_wrap(~sample_name_machine) +
    scale_color_manual("Proteoform", values = c("darkorange", "steelblue")) +
    scale_x_continuous(breaks = c(40, 50, 60)) +
    labs(x = x_label, y = y_label) +
    theme_paper +
    theme(legend.position = "bottom")

ggsave(POLR3A_1_POLR3B_3_profiles, filename = 
           file.path(figure_output_folder, "suppl_fig_POLR3A_1_POLR3B_3_profiles.pdf"),
       width = 10.5, height = 10, units = "cm")

PSMA3_2_PSMA7_1_profiles <- proteoform_df %>% 
    filter(grepl("PSMA3_2|PSMA7_1", gene),
           !grepl("BR2", sample_name_machine)) %>%
    mutate(sample_name_machine = gsub("_", "-", sample_name_machine)) %>% 
    ggplot(aes(temperature, rel_value)) + 
    geom_point(aes(color = gene), alpha = 0.5) + 
    facet_wrap(~sample_name_machine) +
    scale_color_manual("Proteoform", values = c("darkorange", "steelblue")) +
    scale_x_continuous(breaks = c(40, 50, 60)) +
    labs(x = x_label, y = y_label) +
    theme_paper +
    theme(legend.position = "bottom")

ggsave(PSMA3_2_PSMA7_1_profiles, filename = 
           file.path(figure_output_folder, "suppl_fig_PSMA3_2_PSMA7_1_profiles.pdf"),
       width = 10.5, height = 10, units = "cm")

EIF3F_2_EIF3I_1_profiles <- proteoform_df %>% 
    filter(grepl("EIF3F_2|EIF3I_1", gene),
           !grepl("BR2", sample_name_machine)) %>%
    mutate(sample_name_machine = gsub("_", "-", sample_name_machine)) %>% 
    ggplot(aes(temperature, rel_value)) + 
    geom_point(aes(color = gene), alpha = 0.5) + 
    facet_wrap(~sample_name_machine) +
    scale_color_manual("Proteoform", values = c("darkorange", "steelblue")) +
    scale_x_continuous(breaks = c(40, 50, 60)) +
    labs(x = x_label, y = y_label) +
    theme_paper +
    theme(legend.position = "bottom")

ggsave(EIF3F_2_EIF3I_1_profiles, filename = 
           file.path(figure_output_folder, "suppl_fig_EIF3F_2_EIF3I_1_profiles.pdf"),
       width = 10.5, height = 10, units = "cm")

EXOC2_3_EXOC4_3_profiles <- proteoform_df %>% 
    filter(grepl("EXOC2_3|EXOC4_3", gene),
           !grepl("BR2", sample_name_machine)) %>%
    mutate(sample_name_machine = gsub("_", "-", sample_name_machine)) %>% 
    ggplot(aes(temperature, rel_value)) + 
    geom_point(aes(color = gene), alpha = 0.5) + 
    facet_wrap(~sample_name_machine) +
    scale_color_manual("Proteoform", values = c("darkorange", "steelblue")) +
    scale_x_continuous(breaks = c(40, 50, 60)) +
    labs(x = x_label, y = y_label) +
    theme_paper +
    theme(legend.position = "bottom")

ggsave(EXOC2_3_EXOC4_3_profiles, filename = 
           file.path(figure_output_folder, "suppl_fig_EXOC2_3_EXOC4_3_profiles.pdf"),
       width = 10.5, height = 10, units = "cm")

MDM2_2_TP53_1_proteoform_df <- proteoform_df %>% 
    filter(grepl("MDM2_2|TP53_1", gene),
           !grepl("BR2", sample_name_machine)) %>%
    na.omit() %>% 
    group_by(sample_name_machine) %>% 
    filter(n() == 16) %>% 
    ungroup() 

MDM2_2_TP53_1_profiles <- MDM2_2_TP53_1_proteoform_df %>% 
    mutate(sample_name_machine = gsub("_", "-", sample_name_machine)) %>% 
    ggplot(aes(temperature, rel_value)) + 
    geom_point(aes(color = gene), alpha = 0.5) + 
    facet_wrap(~sample_name_machine, nrow = 3) +
    scale_color_manual("Proteoform", values = c("darkorange", "steelblue")) +
    scale_x_continuous(breaks = c(40, 50, 60)) +
    labs(x = x_label, y = y_label) +
    theme_paper +
    theme(legend.position = "bottom")

ggsave(MDM2_2_TP53_1_profiles, filename = 
           file.path(figure_output_folder, "suppl_fig_MDM2_2_TP53_1_profiles.pdf"),
       width = 8, height = 8, units = "cm")


# get Eucl. distances for MDM2_2:TP53_1 in different cell lines
MDM2_2_TP53_1_dist_df <- bind_rows(lapply(names(tpca_result_list), function(nm){
    res_df <- tpca_result_list[[nm]]
    filter(res_df, complex_name == "MDM2_2:TP53_1") %>% 
        mutate(sample = nm)
}))

# test for differences in drug sensitivity of nucleosid analagon drugs
MDM2_2_TP53_1_dist_dss_df <- MDM2_2_TP53_1_dist_df %>% 
    left_join(dss_df %>% filter(drug_name %in% c("Idasanutlin")), 
              by = "sample") %>% 
    na.omit() %>% 
    filter(sample %in% MDM2_2_TP53_1_proteoform_df$sample_name_machine)

MDM2_2_TP53_1_dist_Idasa_boxpl <- 
    ggplot(MDM2_2_TP53_1_dist_dss_df, aes(mean_dist < 0.1, dss)) + 
    geom_boxplot() + 
    ggsignif::geom_signif(comparisons = list(c("TRUE", "FALSE"))) + 
    geom_jitter(aes(color = sample), width = 0.05) +
    scale_color_manual("Cell line", values = cl_colors) +
    labs(x = "Eucledian distance MDM2_2:TP53_1 < 0.1",
         y = "Drug sensitivity score") +
    facet_wrap(~drug_name) +
    theme_paper +
    theme(legend.position = "none")

ggsave(MDM2_2_TP53_1_dist_Idasa_boxpl, filename = 
           file.path(figure_output_folder, "suppl_fig_MDM2_2_TP53_1_dist_Idasa_boxpl.pdf"),
       width = 8, height = 7, units = "cm")
