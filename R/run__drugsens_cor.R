library(tidyverse)
library(Biobase)
library(preprocessCore)
library(limma)
library(here)

sourceDir <- function(path, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        source(file.path(path, nm), ...)
    }
}
sourceDir(path = file.path(here("R/drugsens_cor")))

# define ouput folder
output_folder <- here("drugsens_cor", "output")
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# read in meta data file
sample_meta_file <- here("meta/sample_meta.txt")
sample_meta_raw <- read_tsv(file = sample_meta_file)

# read in and preprocess drugsens data
dss_df <- read_delim(
    here("data/2020-09-09_sDSS_2D_matrix_20_with_annotation_Meltome.txt"), 
    delim = "\t") %>% 
    dplyr::select(-Mechanism.Targets, -Putative.Target.Protein,
                  -Class.explained, -High.phase.Approval.status,
                  -Res_code, -Alias, -Activity.modifier,
                  -Active_inactive.in.clinic, -Solvent,
                  -High.conc..nM., -InChI, -ChEMBL.ID) %>% 
    gather(sample, dss, -FIMM.ID, -DRUG.NAME) %>% 
    dplyr::select(id = FIMM.ID, drug_name = DRUG.NAME, sample, dss) %>% 
    mutate(sample = gsub("\\.", "_", sample)) %>% 
    mutate(sample = gsub("X_", "", sample)) %>% 
    mutate(sample = gsub("Kasumi", "KASUMI", sample)) %>% 
    # filter out COG_319 due to bad quality
    filter(sample != "COG_319")

# read in NPARC result tables
nparc_res_hq_df <- readRDS(here("nparc/output/standard/nparc_res_hq_df.RDS"))
nparc_fstat_df <- readRDS(here("nparc/output/standard/nparc_fstat_df.RDS"))

# create auc df
auc_df <- nparc_res_hq_df %>% 
    filter(id %in% filter(nparc_fstat_df, F_statistic > 
                              quantile(nparc_fstat_df$F_statistic, .9))$id) %>% 
    dplyr::select(sample = sample_name, id, aumc) 

# convert to matrix
auc_mat <- auc_df %>% 
    dplyr::select(id, aumc, sample) %>% 
    spread(sample, aumc) %>% 
    column_to_rownames("id") %>% 
    as.matrix()

# normalize AUC values
auc_mat_norm <- preprocessCore::normalize.quantiles(auc_mat)
dimnames(auc_mat_norm) <- dimnames(auc_mat)

# filter for minimal dss effect
dss_df <- filter_drugsens_for_minimal_effect(dss_df, min_effect = 6)

# convert dss table to matrix
dss_mat <- dss_df %>% 
    dplyr::select(sample, drug_name, dss) %>% 
    spread(drug_name, dss) %>% 
    column_to_rownames("sample") %>% 
    as.matrix()

# match ids of auc and dss matices and adapt auc matrix rownames
auc_matched_ids <- match(rownames(dss_mat), colnames(auc_mat_norm))
rownames(auc_mat_norm) <- gsub("-", "_", rownames(auc_mat_norm))

# loop across all proteoforms and perform limma analysis
all_limma_out_df <- bind_rows(lapply(rownames(auc_mat_norm), function(prot){
    oneProtDesignMatrix <- model.matrix(as.formula(paste("~1 +", prot)), t(auc_mat_norm[,auc_matched_ids]) %>% as.data.frame())
    fit <- lmFit(t(dss_mat[rownames(dss_mat) %in% rownames(oneProtDesignMatrix),]), 
                 design = oneProtDesignMatrix)
    fit2 <- eBayes(fit)
    resTab <- topTable(fit2, number = Inf, sort.by = "p") %>% 
        rownames_to_column() %>% 
        tbl_df %>% 
        mutate(testedProt = prot)
    return(resTab)
})) %>% 
    as_tibble() %>% 
    arrange(P.Value) %>% 
    mutate(p_adj = p.adjust(P.Value, "BH"))

# print out significant hits
print(all_limma_out_df %>% filter(p_adj < 0.1), n = 30)

# volcano plot
ggplot(filter(all_limma_out_df, p_adj >= 0.1), 
       aes(logFC, -log10(P.Value))) +
    stat_binhex(aes(alpha = log10(..count..)), bins = 100, fill = "black") +
    geom_point(alpha = 0.5, color = "black", 
               data = filter(all_limma_out_df, p_adj < 0.1)) +
    labs(x = bquote('log'[2]*'(fold change)'),
         y = bquote('-log'[10]*'('*italic(p)*'-value)')) +
    theme(legend.position = "none")
