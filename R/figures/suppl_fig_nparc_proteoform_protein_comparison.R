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

# define nparc analysis output folder
proteoform_nparc_folder <- here("nparc/output/standard/")
protein_nparc_folder <- here("nparc/output/protein_level/")

proteoform_fstat_df <- readRDS(here(proteoform_nparc_folder, "nparc_fstat_df.RDS")) %>% 
    mutate(proteoform_id = id) %>% 
    mutate(id = sub("_.+", "", id)) %>% 
    mutate(proteoform_90_quantile = F_statistic > quantile(F_statistic, 0.9))

protein_fstat_df <- readRDS(here(protein_nparc_folder, "nparc_fstat_protein_level_df.RDS")) %>% 
    mutate(protein_90_quantile = F_statistic > quantile(F_statistic, 0.9))

combo_fstat_df <- left_join(
    proteoform_fstat_df %>% dplyr::select(id, proteoform_id, proteoform_tm = tm, 
                                          proteoform_F_statistic = F_statistic,
                                          proteoform_90_quantile),
    protein_fstat_df %>% dplyr::select(id, tm, F_statistic, protein_90_quantile),
    by = "id"
)


ggplot(combo_fstat_df, aes(F_statistic, proteoform_F_statistic)) +
    geom_point(alpha = 0.2) +
    geom_abline(slope = 1, linetype = "dashed", color = "gray") +
    labs(x = "F statistic, protein-level", y = "F statistic, proteoform-level") +
    coord_fixed() +
    theme_paper

ggsave(filename = "R/figures/suppl_fig_nparc_proteoform_vs_protein_scatter.pdf", 
       width = 10, height = 10, units = "cm")
