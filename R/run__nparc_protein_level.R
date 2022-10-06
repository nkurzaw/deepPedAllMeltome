library(dplyr)
library(tidyr)
library(readr)
library(Biobase)
library(NPARC)
library(here)

sourceDir <- function(path, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        source(file.path(path, nm), ...)
    }
}
sourceDir(path = file.path(here("R/nparc")))

# set up parallelisation
BPPARAM <- BiocParallel::MulticoreParam(workers = 5, progressbar = TRUE)
# BPPARAM <- BiocParallel::MulticoreParam(workers = 20, progressbar = TRUE)

# the output folder
output_folder <- file.path(here(), "nparc", "output", "protein_level")
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set path for meta data
sample_meta_file <- here("meta/sample_meta.txt")
sample_meta_raw <- read_tsv(file = sample_meta_file)

# read in proteoform ratios
protein_ratios <- readRDS("/g/savitski/01_users/Nils/deepmeltome_output_folder/protein_ratios.RDS")
unique_sample_ids <- unique(pData(protein_ratios)$sample_name_machine)

# loop over samples and compute per sample sigmoid fit statistics using NPARC 
nparc_res_list <- lapply(unique_sample_ids, function(samp_id){
    print(samp_id)
    file_name <- paste0("nparc_res_", gsub(" ", "_", samp_id))
    
    nparc_res <- fit_and_eval_melting_curves(
        proteins = protein_ratios, 
        sample_meta = sample_meta_raw,
        by_var = "cellline", 
        subset_var = samp_id,
        BPPARAM = BPPARAM)
    
    saveRDS(nparc_res, file = file.path(output_folder, paste0(file_name, "_protein_level.RDS")))
    return(nparc_res)
})
names(nparc_res_list) <- unique_sample_ids

# get table of alternative model results per sample
nparc_res_df <- bind_rows(nparc_res_list, .id = "sample_name")

# filter out high noise single-sample fits to get appropriate null models
alt_hq_nparc_df <- nparc_res_df %>% 
    filter(conv, !grepl("_BR2", sample_name), resid_sd < 0.1)
saveRDS(alt_hq_nparc_df, file = file.path(output_folder, "nparc_res_hq_protein_level_df.RDS"))

# remove biological replicates from analysis 
# (they are only used to assess alternative model replicability)
proteins_ratios_wo_rep <- 
    protein_ratios[,-grep("_BR2", pData(protein_ratios)$sample_name_machine)]

# compute null model
nparc_null_res_df <- fit_and_eval_melting_curves(
    proteins = proteins_ratios_wo_rep, 
    sample_meta = sample_meta_raw,
    by_var = "none",
    filter_based_on_alt_model = TRUE,
    alt_model_df = alt_hq_nparc_df,
    BPPARAM = BPPARAM)

saveRDS(nparc_null_res_df, file = file.path(output_folder, "nparc_res_hq_only_null_model_protein_level.RDS"))

# summarize alternative model fit results
nparc_res_df_summarized <- alt_hq_nparc_df %>% 
    dplyr::select(id, rss, nCoeffs, nFitted) %>% 
    group_by(id) %>% 
    dplyr::summarize(rssAlternative = sum(rss),
                     nCoeffsAlternative = sum(nCoeffs),
                     nFittedAlternative = sum(nFitted),
                     .groups = "keep") %>% 
    ungroup

# join full and null model table and compute F statistic
nparc_fstat_df <- left_join(
    nparc_null_res_df %>% 
        filter(conv), 
    nparc_res_df_summarized, 
    by = "id") %>% 
    filter(nFitted == nFittedAlternative,
           nFitted >= 80) %>% 
    # only consider proteoforms that were found in at 
    # least half of the cell lines
    dplyr::select(-nFittedAlternative) %>% 
    mutate(d1 = nCoeffsAlternative - nCoeffs,
           d2 = nFitted - nCoeffs) %>% 
    mutate(F_statistic = ((rss - rssAlternative)/rssAlternative) * (d2/d1)) 

saveRDS(nparc_fstat_df, file = file.path(output_folder, "nparc_fstat_protein_level_df.RDS"))

# volcano plot
ggplot(nparc_fstat_df, aes(rss-rssAlternative, F_statistic)) +
    geom_point(alpha = 0.25, color = "gray") +
    geom_point(color = "black", alpha = 0.5,
               data = filter(nparc_fstat_df, 
                             F_statistic >= quantile(nparc_fstat_df$F_statistic, 0.9))) +
    # geom_label_repel(
    #     label = "NEK2",  
    #     nudge_x = 1,
    #     direction = "x",
    #     segment.size = 0.25,
    #     color = "black", 
    #     data = filter(nparc_fstat_df,id == "NEK2_0")) +
    scale_x_log10() +
    labs(x = bquote('RSS'^0~' - RSS'^1~''),
         y = expression(''*italic(F)*'-statistic')) +
    coord_cartesian(xlim = c(0.05, 17.5)) 
