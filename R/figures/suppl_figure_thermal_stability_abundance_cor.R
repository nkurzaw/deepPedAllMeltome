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
sourceDir(path = here("R/nparc"))
sourceDir(path = here("R/pepnet"))

output_dir <- here("nparc/thermal_stability_abundance")
if(!exists(output_dir))
    dir.create(output_dir)

# get sample meta data
sample_meta_raw <- read_tsv(here("meta/sample_meta.txt"))

# set up parallelisation
BPPARAM <- BiocParallel::MulticoreParam(workers = 5, progressbar = TRUE)

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
    dplyr::summarize(value = sum(value, na.rm = TRUE)) %>% 
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
    left_join(pdata_protein_ratios_df, by = c("sample" = "sample_id")) %>%
    filter(value > 0) %>%
    group_by(gene, sample_name_machine) %>%
    filter(n() == 8) %>%
    mutate(rel_value = value / value[temperature == 41]) %>%
    ungroup

# # define unique saome
unique_sample_ids <- unique(protein_ratios_anno_df$sample_name_machine)

control = NPARC:::getParams()

reRun <- FALSE

if(reRun){
    lapply(unique_sample_ids, function(samp){
        print(samp)
        protein_ratios_sub_df <- protein_ratios_anno_df %>% 
            filter(sample_name_machine == samp)
        
        nparc_fit_res <- NPARC:::invokeParallelFits(
            x = protein_ratios_sub_df$temperature, 
            y = protein_ratios_sub_df$rel_value, 
            id = protein_ratios_sub_df$gene, 
            groups = protein_ratios_sub_df$sample_name_machine,
            BPPARAM = BPPARAM,
            maxAttempts = control$maxAttempts,
            returnModels = FALSE,
            start = control$start)
        
        save(nparc_fit_res, file =
                 paste0(output_dir, "/",
                        "nparc_fit_res_03_12", samp, ".RData"))
    })
}

protein_aumc_df <- bind_rows(lapply(list.files(output_dir, full.names = TRUE), function(fpath){
    sample_name <- sub("\\.RData", "", sub("nparc_fit_res_03_12", "", basename(fpath)))
    load(fpath)
    nparc_fit_df <- 
        nparc_fit_res$modelMetrics %>% 
        dplyr::select(id, sample_name_machine = group,
                      tm, aumc, resid_sd, conv) %>% 
        filter(conv)
    return(nparc_fit_df)
})) %>% 
    group_by(id) %>% 
    mutate(norm_aumc = aumc / mean(aumc, na.rm = TRUE)) %>% 
    ungroup

# read in qMS data and make abundance boxplot for FBP1
qms_eset <- readRDS(here("data/proteins.B.RDS"))

# make proteoform data frame
qms_df <- biobroom::tidy.ExpressionSet(
    qms_eset) 

# get qms sample annotation
qms_pdata_df <- pData(qms_eset) %>% 
    as_tibble() %>% 
    dplyr::select(sample = proteomics_id,
                  sample_name = Cell_Line_Name_Paper) %>% 
    mutate(sample_name_machine = sub("h", "", sub("_LL", "", gsub("-", "_", sample_name)))) %>% 
    filter(!grepl("_BR2", sample_name_machine)) %>% 
    dplyr::select(-sample_name)

# join anno with qms data
qms_anno_df <- left_join(qms_df, qms_pdata_df, by = "sample")

# join 
tpp_qms_df <- protein_aumc_df %>% 
    left_join(qms_anno_df, by = c("id" = "gene", "sample_name_machine"))

ggplot(tpp_qms_df, aes(value, log2(norm_aumc))) + 
    ggpointdensity::geom_pointdensity() + 
    ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE) + 
    labs(x = bquote('log'[2]*'relative abundance fold change'), 
         y = bquote('log'[2]*'relative thermal stability fold change')) +
    theme_paper
    
    
