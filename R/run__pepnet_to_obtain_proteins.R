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

# set up parallelisation
BPPARAM <- BiocParallel::MulticoreParam(workers = 2, progressbar = TRUE)

sourceDir <- function(path, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        source(file.path(path, nm), ...)
    }
}
sourceDir(path = file.path("R/pepnet"))

set2one <- function(x) 1

graphs <- readRDS(file = file.path(output_folder, "graphs_comms.RDS"))

graphs_only_one <- lapply(graphs, function(gr){
    community <- get.graph.attribute(gr, "communities")
    community$membership <- sapply(community$membership, set2one)
    gr <- set.graph.attribute(gr, "communities", community)
    return(gr)
})

peptides_raw <- readRDS(file = file.path(output_folder, "peptides_raw.RDS"))

proteoforms_intensities <- aggregate_peptides_to_proteoforms(e_set = peptides_raw,
                                                             graphs = graphs_only_one,
                                                             aggregation_fun = sum,
                                                             BPPARAM = BPPARAM)