library(tidyverse)
library(Biobase)
library(ggthemes)
library(NPARC)
library(BiocParallel)
library(cowplot)
library(GGally)
library(here)

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

naprc_result_path <- here("nparc/output/standard/")

nparc_res_files <- list.files(naprc_result_path, pattern = "nparc_res_[A-Z,6]", full.names = TRUE)
nparc_res_list <- lapply(nparc_res_files, function(file){
    readRDS(file)
})
names(nparc_res_list) <- sub("nparc_res_", "", sub(".RDS", "", basename(nparc_res_files)))

nparc_res_df <- bind_rows(nparc_res_list, .id = "sample_name") %>% 
    filter(conv, resid_sd < 0.1)

nparc_res_auc_spread_df <- nparc_res_df %>% 
    dplyr::select(sample_name, id, aumc) %>% 
    spread(sample_name, aumc)

cor_mat <- cor(nparc_res_auc_spread_df[,-1], use = "complete.obs")

lowerFn <- function(data, mapping, dot.size) {
    p <- ggplot(data = data, mapping = mapping) +
        geom_point(alpha = 0.15, size = 0.5) +
        geom_line(aes(x,y), data = tibble(x = 0:35, y = 0:35),
                  color = "darkgray") +
        ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE, size = 2.5) +
        theme_bw() +
        theme(axis.text = element_text(size = 6))
    p
}

GGally::ggpairs(nparc_res_auc_spread_df[,c(2, 5, 18, 19)],
                upper = "blank", 
                lower = list(continuous = GGally::wrap(lowerFn)))

ggsave("R/figures/suppl_figure_AUC_scatterplots_replicates_other.pdf", 
       width = 12, height = 12, units = "cm")

# compare to human meltome atlas data analysis
human_meltome_atlas_df <- readRDS("~/repos/meltomeAtlasApp/data/humanCellsUpdated.rds")
human_meltome_tm_df <- human_meltome_atlas_df %>% 
    filter(temperature == 37) %>% 
    dplyr::select(gene_name, cell_line_or_type, meltPoint)  %>% 
    na.omit() %>% 
    group_by(gene_name) %>% 
    mutate(n = n(), n_cell_lines = length(unique(cell_line_or_type))) %>% 
    ungroup() %>% 
    filter(n_cell_lines > 4, n > 11)

## ANOVA
anova_pval_list <- lapply(unique(human_meltome_tm_df$gene_name), function(gene){
    tmp_df <- filter(human_meltome_tm_df, gene_name == gene)
    anova_test <- aov(tmp_df$meltPoint ~ tmp_df$cell_line_or_type)
    test_summary <- summary(anova_test)
    return(test_summary[[1]]$`Pr(>F)`[1])
})
names(anova_pval_list) <- unique(human_meltome_tm_df$gene_name)

anova_pval_df <- bind_cols(anova_pval_list, .id = "gene_name") %>% 
    gather(gene_name, p_value) %>% 
    filter(gene_name != ".id") %>% 
    mutate(p_value = as.numeric(p_value)) %>% 
    mutate(p_adj = p.adjust(p_value))

## NPARC
human_meltome_fc_df <- human_meltome_atlas_df %>% 
    filter(gene_name %in% human_meltome_tm_df$gene_name) %>% 
    group_by(gene_name, cell_line_or_type, temperature) %>% 
    dplyr::summarize(rel_value = mean(fold_change, na.rm = TRUE)) %>% 
    ungroup()

# set up parallelisation
BPPARAM <- BiocParallel::MulticoreParam(workers = 2, progressbar = TRUE)
# BPPARAM <- BiocParallel::MulticoreParam(workers = 20, progressbar = TRUE)

# the output folder
#output_folder <- file.path(here(), "nparc", "output", "standard")
output_folder <- file.path(here(), "nparc", "output", "meltome_atlas")
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

unique_sample_ids <- unique(human_meltome_fc_df$cell_line_or_type)

# loop over samples and compute per sample sigmoid fit statistics using NPARC 
nparc_meltome_res_list <- lapply(unique_sample_ids, function(samp_id){
    print(samp_id)
    file_name <- paste0("nparc_res_", gsub(" ", "_", samp_id))
    
    proteins_df <- filter(human_meltome_fc_df, cell_line_or_type == samp_id)
    
    control = NPARC:::getParams()
    
    fit_res <- NPARC:::invokeParallelFits(
        x = proteins_df$temperature, 
        y = proteins_df$rel_value, 
        id = proteins_df$gene_name, 
        groups = NULL,
        BPPARAM = BPPARAM,
        maxAttempts = control$maxAttempts,
        returnModels = FALSE,
        start = control$start)
    
    nparc_res <- fit_res$modelMetrics
    
    saveRDS(nparc_res, file = file.path(output_folder, paste0(file_name, ".RDS")))
    return(nparc_res)
})
names(nparc_meltome_res_list) <- unique_sample_ids

# get table of alternative model results per sample
nparc_meltome_res_df <- bind_rows(nparc_meltome_res_list, .id = "sample_name")

# filter out high noise single-sample fits to get appropriate null models
alt_hq_nparc_meltome_df <- nparc_meltome_res_df %>% 
    filter(conv, resid_sd < 0.1)
saveRDS(alt_hq_nparc_meltome_df, file = file.path(output_folder, "nparc_meltome_res_hq_df.RDS"))

# compute null model
control = NPARC:::getParams()

null_fit_res <- NPARC:::invokeParallelFits(
    x = human_meltome_fc_df$temperature, 
    y = human_meltome_fc_df$rel_value, 
    id = human_meltome_fc_df$gene_name, 
    groups = NULL,
    BPPARAM = BPPARAM,
    maxAttempts = control$maxAttempts,
    returnModels = FALSE,
    start = control$start)

nparc_null_meltome_res_df <- null_fit_res$modelMetrics

saveRDS(nparc_null_meltome_res_df, 
        file = file.path(output_folder, "nparc_meltome_res_hq_only_null_model.RDS"))


# summarize alternative model fit results
nparc_meltome_res_df_summarized <- alt_hq_nparc_meltome_df %>% 
    dplyr::select(id, rss, nCoeffs, nFitted) %>% 
    group_by(id) %>% 
    dplyr::summarize(rssAlternative = sum(rss),
                     nCoeffsAlternative = sum(nCoeffs),
                     nFittedAlternative = sum(nFitted),
                     .groups = "keep") %>% 
    ungroup

# join full and null model table and compute F statistic
nparc_meltome_fstat_df <- left_join(
    nparc_null_meltome_res_df %>% 
        filter(conv), 
    nparc_meltome_res_df_summarized, 
    by = "id") %>% 
    filter(nFitted == nFittedAlternative,
           nFitted >= 40) %>% 
    # only consider proteoforms that were found in at 
    # least half of the cell lines
    dplyr::select(-nFittedAlternative) %>% 
    mutate(d1 = nCoeffsAlternative - nCoeffs,
           d2 = nFitted - nCoeffs) %>% 
    mutate(F_statistic = ((rss - rssAlternative)/rssAlternative) * (d2/d1)) 

saveRDS(nparc_meltome_fstat_df, file = file.path(output_folder, "nparc_meltome_fstat_df.RDS"))


fstat_90_quan <- quantile(nparc_meltome_fstat_df$F_statistic, 0.9)

nparc_volcano <- 
    ggplot(nparc_meltome_fstat_df, aes(rss-rssAlternative, F_statistic)) +
    geom_point(alpha = 0.25, color = "gray", size = 0.5) +
    geom_point(color = "black", alpha = 0.5, size = 0.5,
               data = filter(nparc_meltome_fstat_df, F_statistic >= fstat_90_quan)) +
    geom_point(color = "red", alpha = 0.5, size = 0.5,
               data = filter(nparc_meltome_fstat_df, F_statistic >= fstat_90_quan,
                             id %in% filter(anova_pval_df, p_adj < 0.1)$gene_name)) +
    scale_x_log10() +
    labs(x = bquote('RSS'^0~' - RSS'^1~''),
         y = expression(''*italic(F)*'-statistic')) +
    coord_cartesian(xlim = c(0.05, 17.5)) +
    theme_paper

nparc_anova_combo_df <- anova_pval_df %>% 
    left_join(nparc_meltome_fstat_df %>% dplyr::select(gene_name = id, F_statistic),
              by = "gene_name")

ggplot(nparc_anova_combo_df, aes(F_statistic, -log10(p_value))) + 
    geom_point(alpha = 0.2, size = 0.75) + 
    geom_hline(yintercept = -log10(0.0000213), color = "darkorange") + 
    geom_vline(xintercept = fstat_90_quan, color = "darkorange") + 
    theme_paper

ggsave("R/figures/suppl_figure_NPARC_vs_ANOVA_meltome_atlas.pdf", width = 9, height = 9, units = "cm")
