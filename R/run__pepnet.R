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

# workaround to install leiden on tamarindo
# library(reticulate)
# py_install("python-igraph")
# py_install("leidenalg", forge = TRUE)

# seems that the leiden package only works from the github repo
# devtools::install_github("TomKellyGenetics/leiden", ref = "master")

library(leiden)

# workaround for vsn incompatibility with other packages
# BiocManager::install("vsn")

sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = file.path("R/pepnet"))

# set up parallelisation
BPPARAM <- BiocParallel::MulticoreParam(workers = 2, progressbar = TRUE)
# BPPARAM <- BiocParallel::MulticoreParam(workers = 20, progressbar = TRUE)


# set the membership colors
palette <- c("#FF5376", "#72AFD9", "#E3D26F", "#A288E3", "#1B5299", "#68D8D6", "#B78DA3")
membership_colors <- get_color_vector(colors = palette, vec = seq_len(8))

# the output folder
output_folder <- here("proteoform_detection", "output", "standard")
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)


############### IMPORT ###############

# import PSMs
psms <- import_psms_from_nf(file = here("data", "target_psmtable.txt"),
                            id_col = "Gene Name",
                            protein_id_col = "Protein",
                            peptide_col = "Peptide",
                            quan_regex = "tmt16plex")

# sum to peptides
peptides_raw <- psms_to_peptides(psms = psms,
                                 summarise_fun = sum,
                                 sample_meta_file = here("meta", "meltome_sample_meta.txt"),
                                 sample_id_col = "sample_id")

saveRDS(object = peptides_raw, file = file.path(output_folder, "peptides_raw.RDS"))

# peptides_raw <- readRDS(file = file.path(output_folder, "peptides_raw.RDS"))

# # VSN normalisation
peptides_vsn_norm <- vsn_normalize_by_temperature(e_set = peptides_raw)

# ratios to lowest temperature
peptides <- build_ratios_to_lowest_temperature(e_set = peptides_vsn_norm, 
                                               sample_col = "sample_name")

# extract first protein IDs
fData(peptides) <- peptides %>%
  fData() %>%
  separate(protein_ids,
           sep = ";",
           into = "first_protein_id",
           extra = "drop",
           remove = FALSE)

# combine sample name and temperature
pData(peptides)$sample_name_temp <- paste0(pData(peptides)$sample_name, "_", pData(peptides)$temperature)

# store
saveRDS(object = peptides, file = file.path(output_folder, "peptides.RDS"))

# peptides <- readRDS(file = file.path(output_folder, "peptides.RDS"))


############### ANALYSIS ###############

# similarity analysis (euclidean distance)
similarities <- evaluate_similarity(e_set = peptides,
                                    filter_params = list(min_num_peptides_per_ioi = 10,
                                                         max_num_peptides_per_ioi = Inf,
                                                         min_peptides_per_sample = 2,
                                                         min_samples_with_sufficient_peptides = 20),
                                    include_ambiguous_ids = TRUE,
                                    method = "euclidean",
                                    transform_fun = function (x) 1 / (1 + x),
                                    BPPARAM = BPPARAM)

saveRDS(object = similarities, file = file.path(output_folder, "similarities.RDS"))

# similarities <- readRDS(file = file.path(output_folder, "similarities.RDS"))

# build graph
graphs <- build_graphs(similarities = similarities,
                       e_set = peptides,
                       filter_params = list(lower_similarity_cutoff = 0,
                                            lower_n_cutoff = 20,
                                            upper_q_cutoff = Inf),
                       BPPARAM = BPPARAM)

# store
saveRDS(object = graphs, file = file.path(output_folder, "graphs.RDS"))

# detect communities
graphs <- detect_communities(graphs = graphs,
                             detect_algorithm = cluster_leiden,
                             BPPARAM = BPPARAM)

# store
saveRDS(object = graphs, file = file.path(output_folder, "graphs_comms.RDS"))

# graphs <- readRDS(file = file.path(output_folder, "graphs_comms.RDS"))


############### AGGREGATION ###############

# filter graphs for 0 modularity
graphs_01 <- graphs[(lapply(graphs, get.graph.attribute, name = "proteoform_modularity") > 0) %>% unlist()]

proteoforms_intensities <- aggregate_peptides_to_proteoforms(e_set = peptides_raw,
                                                             graphs = graphs_01,
                                                             aggregation_fun = sum,
                                                             BPPARAM = BPPARAM)

# filter out small number of peptides and ambiguous only proteoforms
proteoforms_intensities_filtered <- proteoforms_intensities %>%
  .[fData(.)$ambiguous_peptides_only == FALSE, ] %>%
  .[fData(.)$num_peptides > 2] %>%
  .[fData(.)$ambiguity_ratio < 0.5]

# VSN normalisation
proteoforms_vsn_norm <- vsn_normalize_by_temperature(e_set = proteoforms_intensities_filtered)

# ratios to lowest temperature
proteoforms <- build_ratios_to_lowest_temperature(e_set = proteoforms_vsn_norm, 
                                                  sample_col = "sample_name")

saveRDS(object = proteoforms, file = file.path(output_folder, "proteoforms_narrow_range_focused.RDS"))

# proteoforms <- readRDS(file = file.path(output_folder, "proteoforms.RDS"))

iois <- proteoforms %>%
  fData() %>%
  filter(membership != 0) %>%
  .$ioi %>%
  gsub("_[0-9]+$", "", .) %>%
  unique()

# plot
o <- file.path(output_folder, "proteoforms")
if (!dir.exists(o)) dir.create(o, recursive = TRUE)
for (ioi in iois) {
  pdf(file = file.path(o, paste0(ioi, ".pdf")))
  plot_proteoform(ioi = ioi,
                  e_set = proteoforms,
                  x_col = "temperature",
                  facet_col = "sample_name",
                  colors = membership_colors,
                  add_splines = TRUE)
  dev.off()
}


############### EVALUATION ###############

cat("=== Evaluate parameters\n")

# evaluate graph parameters
graph_parameters <- calculate_graph_parameters(graphs = graphs,
                                               BPPARAM = BPPARAM)
write_delim(x = graph_parameters,
            file = file.path(output_folder, "graph_parameters.txt"),
            delim = "\t",
            col_names = TRUE)

# evaluate community parameters
community_parameters <- calculate_community_parameters(graphs = graphs,
                                                       BPPARAM = BPPARAM)
write_delim(x = community_parameters,
            file = file.path(output_folder, "community_parameters.txt"),
            delim = "\t",
            col_names = TRUE)

# combine and get some overviews
parameters_df <- community_parameters %>%
  left_join(graph_parameters, by = "ioi") %>%
  mutate(cutoff = if_else(max_modularity >= 0.2, "0.20",
                          if_else(max_modularity >= 0.1, "0.10",
                                  if_else(max_modularity >= 0.05, "0.05",
                                          if_else(max_modularity >= 0, "0.00", "neg")))))

pdf(file = file.path(output_folder, "community_parameters.pdf"))

lapply(X = colnames(parameters_df)[3:6],
       FUN = function (parameter_name) {
         p <- parameters_df %>%
           ggplot(aes_string(x = parameter_name, group = "cutoff", color = "cutoff")) +
           geom_density() +
           ggtitle(parameter_name)
         
         if (parameter_name == "compactness") {
           p <- p + xlim(0, 100)
         }
         
         p
       })

dev.off()


# evaluate within proteoform similarities
pdf(file = file.path(output_folder, "proteoform_similarities.pdf"))

proteoform_similarities <- get_similarities_per_proteoform(graphs = graphs,
                                                           modularity_cutoffs = c(-Inf, 0, 0.05, 0.1, 0.2),
                                                           do_plot = TRUE,
                                                           BPPARAM = BPPARAM)

dev.off()


############### QC PLOTS ###############

pdf(file = file.path(output_folder, "qc_plots.pdf"))

plot_num_proteoforms_per_id(graphs = graphs, BPPARAM = BPPARAM)

plot_num_peptides_per_proteoform(graphs = graphs, BPPARAM = BPPARAM)

# modularity comparison between protein IDs and memberships
modularities <- compare_modularities(graphs = graphs,
                                     x_col = "membership",
                                     y_col = "protein_ids",
                                     h_line_intercept = 0,
                                     v_line_intercept = 0.2,
                                     BPPARAM = BPPARAM)

dev.off()

write_delim(x = modularities,
            path = file.path(output_folder, "modularity_comparison.txt"),
            delim = "\t",
            col_names = TRUE)


############### VARIOUS PLOTS ###############

# plot graphs and profiles
graphs_to_plot <- graphs %>%
  .[lapply(., get.graph.attribute, "proteoform_modularity") > 0.1] %>%
  .[lapply(., get.graph.attribute, "num_peptides_in_largest_community") > 5]

cat(">>> Plotting", length(graphs_to_plot), "graphs...\n")

plot_combo(graphs = graphs_to_plot,
           e_set = peptides,
           include_ambiguous_ids = TRUE,
           x_col = "sample_name_temp",
           box_group_by = "temperature",
           membership_colors = membership_colors,
           output_folder = file.path(output_folder, "plots"))
