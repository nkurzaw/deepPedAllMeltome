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
       width = 21, height = 10, units = "cm")

# load pre-run proteoform detection datasets

## 15 peptides per protein
### pepnet
graphs <- readRDS(here("R/benchmark/graphs_comms_20_intra_noise.RDS"))

eval_df <- tibble(
    protein_name = names(graphs),
    modularity = sapply(graphs, get.graph.attribute, name = "proteoform_modularity"),
    detected_proteoforms = sapply(graphs, function(x) max(get.graph.attribute(x, "communities")$membership))
) %>% 
    arrange(desc(modularity)) %>% 
    mutate(tp = as.numeric(grepl("tp", protein_name)))

roc_15_obj <- roc(eval_df$tp ~ eval_df$modularity)
auc(roc_15_obj)
plot(roc_15_obj)

pepnet_15_roc_df <- ggroc(roc_15_obj)$data %>% 
    mutate(method = "pepnet")

### COPF
copf_15_scores_df <- read_csv(here("R/benchmark/cc_profiler_proteoform_scores_20_intra_noise.csv")) %>%
    within(proteoform_score[is.na(proteoform_score)] <- 0) %>%
    mutate(tp = as.numeric(grepl("tp", protein_id)))

copf_15_roc_obj <- roc(copf_15_scores_df$tp ~ copf_15_scores_df$proteoform_score)
auc(copf_15_roc_obj)
plot(copf_15_roc_obj)

copf_15_roc_df <- ggroc(copf_15_roc_obj)$data %>% 
    mutate(method = "COPF")

combo_15_roc_df <- bind_rows(pepnet_15_roc_df, copf_15_roc_df)

benchmark_15_pep_roc <- ggplot(combo_15_roc_df, aes(1-specificity, sensitivity, color = method)) +
    geom_line() +
    geom_abline(slope = 1, linetype = "dashed", color="gray") +
    geom_text(aes(x, y, label = label, color = method), 
              data = tibble(x = c(0.75, 0.75),
                            y = c(0.175, 0.25), 
                            method = c("pepnet", "COPF"), 
                            label = c(paste("AUC =", round(auc(roc_15_obj)[1], 3)), 
                                      paste("AUC =", round(auc(copf_15_roc_obj)[1], 3))))) +
    theme_paper +
    coord_fixed() +
    ggtitle("15 peptides per protein") +
    labs(x = "FPR", y = "TPR") +
    theme(legend.position = "bottom")

## 50 peptides per protein
### pepnet
graphs50 <- readRDS(here("R/benchmark/graphs_comms_pep_cov_50_20_intra_noise.RDS"))

eval_50_df <- tibble(
    protein_name = names(graphs50),
    modularity = sapply(graphs50, get.graph.attribute, name = "proteoform_modularity"),
    detected_proteoforms = sapply(graphs50, function(x) max(get.graph.attribute(x, "communities")$membership))
) %>% 
    arrange(desc(modularity)) %>% 
    mutate(tp = as.numeric(grepl("tp", protein_name)))

roc_50_obj <- roc(eval_50_df$tp ~ eval_50_df$modularity)
auc(roc_50_obj)
plot(roc_50_obj)

pepnet_50_roc_df <- ggroc(roc_50_obj)$data %>% 
    mutate(method = "pepnet")

### COPF
copf_50_scores_df <- read_csv(here("R/benchmark/cc_profiler_proteoform_scores_pep_cov_50.csv")) %>%
    within(proteoform_score[is.na(proteoform_score)] <- 0) %>%
    mutate(tp = as.numeric(grepl("tp", protein_id)))

copf_50_roc_obj <- roc(copf_50_scores_df$tp ~ copf_50_scores_df$proteoform_score)
auc(copf_50_roc_obj)
plot(copf_50_roc_obj)

copf_50_roc_df <- ggroc(copf_50_roc_obj)$data %>% 
    mutate(method = "COPF")

combo_50_roc_df <- bind_rows(pepnet_50_roc_df, copf_50_roc_df)

benchmark_50_pep_roc <- ggplot(combo_50_roc_df, aes(1-specificity, sensitivity, color = method)) +
    geom_line() +
    geom_abline(slope = 1, linetype = "dashed", color="gray") +
    geom_text(aes(x, y, label = label, color = method), 
              data = tibble(x = c(0.75, 0.75),
                            y = c(0.175, 0.25), 
                            method = c("pepnet", "COPF"), 
                            label = c(paste0(paste("AUC =", round(auc(roc_50_obj)[1], 3)), "0"), 
                                      paste("AUC =", round(auc(copf_50_roc_obj)[1], 3))))) +
    theme_paper +
    coord_fixed() +
    ggtitle("50 peptides per protein") +
    labs(x = "FPR", y = "TPR") +
    theme(legend.position = "bottom")

plot_grid(
    plot_grid(benchmark_15_pep_roc + theme(legend.position = "none"), 
              benchmark_50_pep_roc + theme(legend.position = "none"),
              labels = letters[1:2]),
    get_legend(benchmark_15_pep_roc), nrow = 2, rel_heights = c(9, 1)
)


ggsave(filename = here("R/figures/suppl_fig_pepnet_benchmark.pdf"), 
       width = 21, height = 10, units = "cm")

# compute FDR, sensitivity for different proteoforms
## 15 peptides
nrow(filter(eval_df, detected_proteoforms == 2))
# 87
nrow(filter(eval_df, detected_proteoforms == 2, grepl("tn_", protein_name)))
# 0
nrow(filter(eval_df, detected_proteoforms == 2, grepl("tp_", protein_name)))
# 87

nrow(filter(eval_df, detected_proteoforms == 2, grepl("tp_protein_4", protein_name)))
# 50
nrow(filter(eval_df, detected_proteoforms == 2, grepl("tp_protein_3", protein_name)))
# 37


## 50 peptides
nrow(filter(eval_50_df, detected_proteoforms == 2))
# 164
nrow(filter(eval_50_df, detected_proteoforms == 2, grepl("tn_", protein_name)))
# 12
nrow(filter(eval_50_df, detected_proteoforms == 2, grepl("tp_", protein_name)))
# 152
#fdr = 12/164 = 0.073

nrow(filter(eval_50_df, detected_proteoforms == 2, grepl("tp_protein_4", protein_name)))
# 50
nrow(filter(eval_50_df, detected_proteoforms == 2, grepl("tp_protein_3", protein_name)))
# 50
nrow(filter(eval_50_df, detected_proteoforms == 2, grepl("tp_protein_2", protein_name)))
# 45
nrow(filter(eval_50_df, detected_proteoforms == 2, grepl("tp_protein_1", protein_name)))
# 7
