library(data.table)
library(tidyverse)
library(Biobase)
library(matrixStats)
library(Hmisc)
library(igraph)
library(BiocParallel)
library(biobroom)
library(cowplot)
library(ggrepel)
library(ggpmisc)
library(RColorBrewer)
library(here)
library(pROC)

library(leiden)

sourceDir <- function(path, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        source(file.path(path, nm), ...)
    }
}
sourceDir(path = file.path(here("R/pepnet")))

#simulated_peptides_pep_cov_15 <- readRDS(here("R/benchmark/simulated_peptides_pep_cov_15_eset.RDS"))
#simulated_peptides_pep_cov_15 <- readRDS(here("R/benchmark/simulated_peptides_pep_cov_15_02_intra_noise.RDS"))
#simulated_peptides_pep_cov_15 <- readRDS(here("R/benchmark/simulated_peptides_pep_cov_15_08_intra_noise.RDS"))
#simulated_peptides_pep_cov_15 <- readRDS(here("R/benchmark/simulated_peptides_pep_cov_15_15_intra_noise.RDS"))
#simulated_peptides_pep_cov_15 <- readRDS(here("R/benchmark/simulated_peptides_pep_cov_15_20_intra_noise.RDS"))
#simulated_peptides_pep_cov_15 <- readRDS(here("R/benchmark/simulated_peptides_pep_cov_15_15_intra_noise_more_hard_cases.RDS"))
simulated_peptides_pep_cov_15 <- readRDS(here("R/benchmark/simulated_peptides_pep_cov_15_15_intra_noise_more_intermediate_hard_cases.RDS"))



BPPARAM <- BiocParallel::MulticoreParam(workers = 4)

# similarity analysis (euclidean distance)
sim_similarities <- evaluate_similarity(e_set = simulated_peptides_pep_cov_15,
                                    filter_params = list(min_num_peptides_per_ioi = 10,
                                                         max_num_peptides_per_ioi = Inf,
                                                         min_peptides_per_sample = 2,
                                                         min_samples_with_sufficient_peptides = 20),
                                    include_ambiguous_ids = TRUE,
                                    method = "euclidean",
                                    transform_fun = function (x) 1 / (1 + x),
                                    BPPARAM = BPPARAM)

#saveRDS(object = sim_similarities, file = here("R/benchmark/sim_similarities.RDS"))
#saveRDS(object = sim_similarities, file = here("R/benchmark/sim_similarities_02_intra_noise.RDS"))
#saveRDS(object = sim_similarities, file = here("R/benchmark/sim_similarities_08_intra_noise.RDS"))
#saveRDS(object = sim_similarities, file = here("R/benchmark/sim_similarities_15_intra_noise.RDS"))
#saveRDS(object = sim_similarities, file = here("R/benchmark/sim_similarities_20_intra_noise.RDS"))
#saveRDS(object = sim_similarities, file = here("R/benchmark/sim_similarities_15_intra_noise_more_hard_cases.RDS"))
saveRDS(object = sim_similarities, file = here("R/benchmark/sim_similarities_15_intra_noise_more_intermediate_hard_cases.RDS"))

# build graph
graphs <- build_graphs(similarities = sim_similarities,
                       e_set = simulated_peptides_pep_cov_15,
                       filter_params = list(lower_similarity_cutoff = 0,
                                            lower_n_cutoff = 20,
                                            upper_q_cutoff = Inf),
                       BPPARAM = BPPARAM)

# store
#saveRDS(object = graphs, file = here("R/benchmark/graphs.RDS"))
#saveRDS(object = graphs, file = here("R/benchmark/graphs_02_intra_noise.RDS"))
#saveRDS(object = graphs, file = here("R/benchmark/graphs_08_intra_noise.RDS"))
#saveRDS(object = graphs, file = here("R/benchmark/graphs_15_intra_noise.RDS"))
#saveRDS(object = graphs, file = here("R/benchmark/graphs_20_intra_noise.RDS"))
#saveRDS(object = graphs, file = here("R/benchmark/graphs_15_intra_noise_more_hard_cases.RDS"))
saveRDS(object = graphs, file = here("R/benchmark/graphs_15_intra_noise_more_intermediate_hard_cases.RDS"))

# detect communities
graphs <- detect_communities(graphs = graphs,
                             detect_algorithm = cluster_leiden,
                             BPPARAM = BPPARAM)

# store
#saveRDS(object = graphs, file = here("R/benchmark/graphs_comms.RDS"))
#saveRDS(object = graphs, file = here("R/benchmark/graphs_comms_02_intra_noise.RDS"))
#saveRDS(object = graphs, file = here("R/benchmark/graphs_comms_08_intra_noise.RDS"))
#saveRDS(object = graphs, file = here("R/benchmark/graphs_comms_15_intra_noise.RDS"))
#saveRDS(object = graphs, file = here("R/benchmark/graphs_comms_20_intra_noise.RDS"))
#saveRDS(object = graphs, file = here("R/benchmark/graphs_comms_15_intra_noise_more_hard_cases.RDS"))
saveRDS(object = graphs, file = here("R/benchmark/graphs_comms_15_intra_noise_more_intermediate_hard_cases.RDS"))

# filter graphs for 0 modularity
graphs_01 <- graphs[(lapply(graphs, get.graph.attribute, name = "proteoform_modularity") > 0) %>% unlist()]

graphs_001 <- graphs[(lapply(graphs, get.graph.attribute, name = "proteoform_modularity") > 1e-14) %>% unlist()]

eval_df <- tibble(
    protein_name = names(graphs),
    modularity = sapply(graphs, get.graph.attribute, name = "proteoform_modularity")
) %>% 
    arrange(desc(modularity)) %>% 
    mutate(tp = as.numeric(grepl("tp", protein_name)))

roc_obj <- roc(eval_df$tp ~ eval_df$modularity)
auc(roc_obj)
plot(roc_obj)

pepnet_roc_df <- ggroc(roc_obj)$data %>% 
    mutate(method = "pepnet")

# COPF benchmark
copf_scores_df <- read_csv(here("R/benchmark/cc_profiler_proteoform_scores_20_intra_noise.csv")) %>% 
    within(proteoform_score[is.na(proteoform_score)] <- 0) %>% 
    mutate(tp = as.numeric(grepl("tp", protein_id)))

# copf_scores_df <- read_csv(here("R/benchmark/cc_profiler_proteoform_scores.csv")) %>% 
#     within(proteoform_score[is.na(proteoform_score)] <- 0) %>% 
#     mutate(tp = as.numeric(grepl("tp", protein_id)))

roc_obj <- roc(copf_scores_df$tp ~ copf_scores_df$proteoform_score)
auc(roc_obj)
plot(roc_obj)

copf_roc_df <- ggroc(roc_obj)$data %>% 
    mutate(method = "COPF")

combo_roc_df <- bind_rows(pepnet_roc_df, copf_roc_df)

ggplot(combo_roc_df, aes(1-specificity, sensitivity, color = method)) +
    geom_line() +
    geom_abline(slope = 1, linetype = "dashed", color="gray") +
    theme_bw()

