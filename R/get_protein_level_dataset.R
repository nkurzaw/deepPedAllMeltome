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
library(data.table)
library(RColorBrewer)
library(here)

sourceDir <- function(path, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        source(file.path(path, nm), ...)
    }
}
sourceDir(path = file.path("R/pepnet"))

proteins <- import_proteins_from_nf(file = here("data/genes_table.txt"),
                                    sample_meta_file = here("meta/meltome_sample_meta.txt"))

# # VSN normalisation
proteins_vsn_norm <- vsn_normalize_by_temperature(e_set = proteins)

saveRDS(proteins_vsn_norm, here("proteoform_detection/output/standard/proteins.RDS"))

## build ratios to lowest temperature
protein_ratios <- build_ratios_to_lowest_temperature(e_set = proteins_vsn_norm, 
                                                     sample_col = "sample_name")

saveRDS(protein_ratios, here("proteoform_detection/output/standard/protein_ratios.RDS"))

proteins_df <- biobroom::tidy.ExpressionSet(protein_ratios, addPheno = TRUE)
