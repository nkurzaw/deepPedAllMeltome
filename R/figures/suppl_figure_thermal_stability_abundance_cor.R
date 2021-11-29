library(data.table)
library(tidyverse)
library(Biobase)
library(biobroom)
library(matrixStats)
library(Hmisc)
library(BiocParallel)
library(NPARC)
library(here)

sourceDir <- function(path, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        source(file.path(path, nm), ...)
    }
}
sourceDir(path = file.path("nparc"))

# get sample meta data
sample_meta_raw <- read_tsv(here("meta/sample_meta.txt"))

# set up parallelisation
BPPARAM <- BiocParallel::MulticoreParam(workers = 2, progressbar = TRUE)

# read in raw peptides
peptides_raw <- readRDS(here("proteoform_detection/output/standard/peptides_raw.RDS"))

# make tidy data frame
peptides_raw_df <- biobroom::tidy.ExpressionSet(peptides_raw, addPheno = FALSE)

# fully annotate data frame
pdata_peptides_df <- pData(peptides_raw) %>% 
    as_tibble %>% 
    dplyr::select(sample_id, temperature, sample_name, sample_name_machine)

fdata_peptides_df <- fData(peptides_raw) %>% 
    as_tibble() %>% 
    dplyr::select(peptide, id)

peptides_anno_raw_df <- peptides_raw_df %>% 
    left_join(pdata_peptides_df, by = c("sample" = "sample_id")) %>% 
    left_join(fdata_peptides_df, by = c("gene" = "peptide"))

protein_df <- peptides_anno_raw_df %>% 
    mutate(id = sub(";.+", "", id)) %>% 
    group_by(id, sample, temperature, sample_name_machine) %>% 
    summarize(value = sum(value, na.rm = TRUE)) %>% 
    ungroup

# make eset out of protein data frame
## make assay may
protein_assay_mat <- protein_df %>% 
    dplyr::select(id, sample, value) %>% 
    spread(sample, value) %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "id") %>% 
    as.matrix

## make pData
protein_pdata_df <- protein_df %>% 
    dplyr::select(sample_id = sample, temperature, sample_name_machine) %>% 
    distinct() %>% 
    mutate(rowname = sample_id) %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "rowname")

## assemble eset
protein_eset <- ExpressionSet(
    assayData = protein_assay_mat,
    phenoData = AnnotatedDataFrame(protein_pdata_df))

# VSN normalisation
proteins_vsn_norm <- vsn_normalize_by_temperature(
    e_set = protein_eset)

# # get ratios
# protein_ratios <- build_ratios_to_lowest_temperature(
#     e_set = proteins_vsn_norm, sample_col = "sample_name_machine")

# get protein ratios df
protein_ratios_df <- tidy.ExpressionSet(proteins_vsn_norm)

# fully annotate data frame
pdata_protein_ratios_df <- pData(proteins_vsn_norm) %>%
    as_tibble

protein_ratios_anno_df <- protein_ratios_df %>%
    left_join(pdata_protein_ratios_df, by = "sample") %>%
    filter(value > 0) %>%
    group_by(gene, sample_name_machine) %>%
    filter(n() == 8) %>%
    mutate(rel_value = value / value[temperature == 41]) %>%
    ungroup

# define unique saome
unique_sample_ids <- unique(pData(protein_ratios_anno_df)$sample_name_machine)
    
control = NPARC:::getParams()

nparc_fit_res <- NPARC:::invokeParallelFits(
    x = protein_ratios_anno_df$temperature, 
    y = protein_ratios_anno_df$rel_value, 
    id = protein_ratios_anno_df$gene, 
    groups = protein_ratios_anno_df$sample_name_machine,
    BPPARAM = BPPARAM,
    maxAttempts = control$maxAttempts,
    returnModels = FALSE,
    start = control$start)

