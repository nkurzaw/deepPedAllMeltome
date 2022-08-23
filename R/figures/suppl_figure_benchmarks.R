library(tidyverse)
library(Biobase)
library(readxl)
library(cowplot)
library(here)
library(igraph)
library(pROC)

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


# load pre-run PPI datasets
reh_ori_ppi_tpca_roc_df <- readRDS(here("R/benchmark/reh_ori_ppi_tpca_roc_df.RDS"))
heusel_sub_inter_ori_ppi_tpca_roc_df <-  readRDS(here("R/benchmark/heusel_sub_inter_ori_ppi_tpca_roc_df.RDS"))

# compute AUCs
#deepmeltome_auc <- auc(roc(reh_ori_ppi_tpca_roc_df$annotated ~ reh_ori_ppi_tpca_roc_df$eucl_dist))
#sec_auc <- auc(roc(heusel_sub_inter_ori_ppi_tpca_roc_df$annotated ~ heusel_sub_inter_ori_ppi_tpca_roc_df$))

# plot Deepmeltome PPI ROC curve
set.seed(123)
deepmeltome_roc <- 
    ggplot(reh_ori_ppi_tpca_roc_df[sort(sample(seq(nrow(reh_ori_ppi_tpca_roc_df)), 10000)),], aes(FPR, TPR)) + 
    geom_path() +
    geom_abline(slope = 1, linetype = "dashed", color = "gray") +
    geom_text(label = "AUC = 0.698", aes(x, y), data = tibble(x = 0.75, y = 0.25)) +
    ggtitle("Deepmeltome cell line REH") +
    theme_paper +
    coord_fixed()

# plot SEC PPI ROC curve
sec_roc <- 
    ggplot(heusel_sub_inter_ori_ppi_tpca_roc_df[sort(sample(seq(nrow(heusel_sub_inter_ori_ppi_tpca_roc_df)), 10000)),], aes(FPR, TPR)) + 
    geom_path() +
    geom_abline(slope = 1, linetype = "dashed", color = "gray") +
    geom_text(label = "AUC = 0.744", aes(x, y), data = tibble(x = 0.75, y = 0.25)) +
    ggtitle("Heusel et al. interphase downsampled SEC-MS") +
    theme_paper +
    coord_fixed()

plot_grid(deepmeltome_roc, sec_roc, labels = letters[1:2], nrow = 1)
ggsave(filename = here("R/figures/suppl_fig_ppi_benchmark.pdf"), 
       width = 21, height = 8, units = "cm")

# load pre-run proteoform detection datasets

## 15 peptides per protein
### pepnet
graphs <- readRDS(here("R/benchmark/graphs_comms_pep_cov_50_20_intra_noise.RDS"))

eval_df <- tibble(
    protein_name = names(graphs),
    modularity = sapply(graphs, get.graph.attribute, name = "proteoform_modularity"),
    detected_proteoforms = sapply(graphs, function(x) max(get.graph.attribute(x, "communities")$membership))
) %>% 
    arrange(desc(modularity)) %>% 
    mutate(tp = as.numeric(grepl("tp", protein_name)))

roc_obj <- roc(eval_df$tp ~ eval_df$modularity)
auc(roc_obj)
plot(roc_obj)

pepnet_15_roc_df <- ggroc(roc_obj)$data %>% 
    mutate(method = "pepnet")

### COPF
copf_15_scores_df <- read_csv(here("R/benchmark/cc_profiler_proteoform_scores_20_intra_noise.csv")) %>%
    within(proteoform_score[is.na(proteoform_score)] <- 0) %>%
    mutate(tp = as.numeric(grepl("tp", protein_id)))

copf_15_roc_obj <- roc(copf_15_scores_df$tp ~ copf_15_scores_df$proteoform_score)
auc(copf_15_roc_obj)
plot(copf_15_roc_obj)

copf_15_roc_df <- ggroc(copf_roc_obj)$data %>% 
    mutate(method = "COPF")

combo_15_roc_df <- bind_rows(pepnet_15_roc_df, copf_15_roc_df)

benchmark_15_pep_roc <- ggplot(combo_15_roc_df, aes(1-specificity, sensitivity, color = method)) +
    geom_line() +
    geom_abline(slope = 1, linetype = "dashed", color="gray") +
    theme_bw()

## 50 peptides per protein

