library(tidyverse)
library(Biobase)
library(ggthemes)
library(NPARC)
library(BiocParallel)
library(cowplot)
library(GGally)
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

naprc_result_path <- here("nparc/output/standard/")

nparc_res_files <- list.files(naprc_result_path, pattern = "nparc_res_[A-Z,6]", full.names = TRUE)
nparc_res_list <- lapply(nparc_res_files, function(file){
    readRDS(file)
})
names(nparc_res_list) <- sub("nparc_res_", "", sub(".RDS", "", basename(nparc_res_files)))

nparc_res_df <- bind_rows(nparc_res_list, .id = "sample_name") %>% 
    filter(conv, resid_sd < 0.1)

nparc_res_auc_spread_df <- nparc_res_df %>% 
    dplyr::select(sample_name, id, aumc) %>% 
    spread(sample_name, aumc)

cor_mat <- cor(nparc_res_auc_spread_df[,-1], use = "complete.obs")

lowerFn <- function(data, mapping, dot.size) {
    p <- ggplot(data = data, mapping = mapping) +
        geom_point(alpha = 0.15, size = 0.5) +
        geom_line(aes(x,y), data = tibble(x = 0:35, y = 0:35),
                  color = "darkgray") +
        ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE, size = 2.5) +
        theme_bw() +
        theme(axis.text = element_text(size = 6))
    p
}

GGally::ggpairs(nparc_res_auc_spread_df[,c(2, 5, 18, 19)],
                upper = "blank", 
                lower = list(continuous = GGally::wrap(lowerFn)))

ggsave("R/figures/suppl_figure_AUC_scatterplots_replicates_other.pdf", 
       width = 12, height = 12, units = "cm")