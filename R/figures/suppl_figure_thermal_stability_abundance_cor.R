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

# load proteoform data
proteoforms <- readRDS(
    here("proteoform_detection/output/standard/proteoforms_narrow_range_focused.RDS"))

# make proteoform data frame
proteoform_df <- biobroom::tidy.ExpressionSet(
    proteoforms, 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

cl_colors <-  c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
                "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", 
                "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", 
                "#771122", "#AA4455")
# Peter Carl - R-Bloggers 02/2013 "The Paul Tol 21-color salute"
names(cl_colors) <- sort(
    unique(gsub("_BR2", "", proteoform_df$sample_name_machine)))

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

global_aumc_abund_cor <- ggplot(tpp_qms_df, aes(value, log2(norm_aumc))) + 
    #ggpointdensity::geom_pointdensity() + 
    #geom_point(alpha = 0.1) +
    stat_binhex(aes(alpha = log10(..count..)), bins = 100, fill = "black") +
    ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE) + 
    labs(x = bquote('log'[2]*'relative abundance fold change'), 
         y = bquote('log'[2]*'relative thermal stability fold change')) +
    theme_paper +
    theme(legend.position = "bottom")

ggsave(global_aumc_abund_cor, 
       filename = here("R/figures/suppl_fig_thermal_stab_abund.pdf"), 
       width = 7, height = 8, units = "cm")
    
# combine with FBP1 figures
### FBP1_1
fbp1_1_df <- filter(proteoform_df, gene == "FBP1_1") %>% 
    filter(!grepl("_BR2", sample_name_machine)) %>% 
    filter(sample_name_machine %in% filter(nparc_res_hq_df, id == "FBP1_1")$sample_name) %>% 
    na.omit()

fbp1_1_alt_fit_param <- NPARC:::invokeParallelFits(
    x = fbp1_1_df$temperature, 
    y = fbp1_1_df$rel_value, 
    id = fbp1_1_df$gene, 
    groups = fbp1_1_df$sample_name_machine,
    BPPARAM = SerialParam(progressbar = TRUE),
    maxAttempts = control$maxAttempts,
    returnModels = FALSE,
    start = control$start)

fbp1_1_alt_fit_df <- 
    tibble(temperature = rep(temp_range, length(unique(fbp1_1_df$sample_name_machine))),
           group = rep(
               unique(fbp1_1_df$sample_name_machine),
               each = length(temp_range)
           )) %>% 
    left_join(fbp1_1_alt_fit_param$modelMetrics, 
              by = "group") %>% 
    rowwise() %>% 
    mutate(y_hat = (1 - pl)  / (1 + exp((b - a/temperature))) + pl) %>% 
    ungroup

fbp1_1_plot_df <- fbp1_1_alt_fit_param$modelPredictions %>% 
    left_join(sample_meta_raw %>% 
                  dplyr::select(group = sample_name_machine, subtype),
              by = "group") 

fbp1_p <- ggplot(fbp1_1_df %>% mutate(group = sample_name_machine), 
                 aes(temperature, rel_value)) +
    geom_line(aes(temperature, y_hat, color = group, group = group), #, linetype = Assigned_stage), 
              data = fbp1_1_alt_fit_df)+#,
    #color = "darkgray") +
    geom_point(aes(color = group)) +
    # geom_segment(aes(xend = x, yend = .fitted), 
    #              linetype = "dashed") +
    #scale_color_manual("", values = cl_colors) +
    scale_color_manual("Sample", values = cl_colors) +
    #scale_x_continuous(breaks = c(40, 55)) +
    labs(x = x_label,
         y = y_label) +
    #facet_wrap(~ subtype, ncol = 11) +
    ggtitle("FBP1 proteoform 1") +
    theme_paper +
    theme(legend.position = "none")


# prepare FBP1-specific data set with thermal stability annotation
fbp1_qms_df <- qms_anno_df %>% 
    filter(gene == "FBP1") %>% 
    mutate(thermal_stability_group = case_when(
        sample_name_machine %in% c("RCH_ACV", "LC4_1", "P30_OHKUBO", "KASUMI_2", 
                                   "MHH_CALL_3", "KOPN_8","COG_319", "RCH_ACV_BR2", 
                                   "LC4_1_BR2", "P30_OHKUBO_BR2", "KASUMI_2_BR2", 
                                   "MHH_CALL_3_BR2", "KOPN_8_BR2","COG_319_BR2") ~ 
            "high",
        TRUE ~ "low")) %>% 
    mutate(thermal_stability_group = factor(
        thermal_stability_group, levels = c("low",
                                            "high"))) %>% 
    filter(cell_line %in% fbp1_1_df$sample_name_machine)

# make boxplot
fbp1_qms_bp <- ggplot(fbp1_qms_df, aes(thermal_stability_group, value)) +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(aes(color = cell_line), width = 0.1) +
    scale_color_manual("Cell line", values = cl_colors) +
    ggsignif::geom_signif(comparisons = list(c("low", "high")),
                          test = t.test) +
    labs(x = "FBP1_1 thermal stability", 
         y = bquote('FBP1 log'[2]*'fold change to mean')) +
    coord_cartesian(ylim = c(-2, 3.5)) +
    theme_paper +
    theme(legend.position = "none")

# check thermal stability of G6PD
g6pd_auc_df <- auc_full_hq_df %>% 
    filter(grepl("G6PD_", id)) %>% 
    filter(!sample %in% c("697","COG_355", "COG_394", "MHH_CALL_2")) %>% 
    mutate(group = ifelse(sample %in%  c("COG_319", "LC4_1", "KASUMI_9", "KASUMI_2", 
                                         "MHH_CALL_3", "P30_OHKUBO", "RCH_ACV", "KOPN_8"), 
                          "high", "low")) %>% 
    mutate(group = factor(group, levels = c("low", "high"))) %>% 
    mutate(sample = gsub("_", "-", sample))

g6pd_boxplot <- ggplot(g6pd_auc_df, aes(group, aumc)) +
    geom_boxplot() +
    geom_jitter(aes(color = sample), width = 0.1) +
    ggsignif::geom_signif(comparisons = list(c("high", "low")), test = t.test) +
    scale_color_manual(values = cl_colors) +
    facet_wrap(~id) +
    labs(x = "FBP1_1 thermal stability",
         y = "Area under the melting curve") +
    coord_cartesian(ylim = c(11, 15.25)) +
    theme_paper +
    theme(legend.position = "none")

plot_grid(plot_grid(fbp1_p, fbp1_qms_bp, 
                    global_aumc_abund_cor + theme(legend.position = "none"),
                    labels = letters[1:3], ncol = 3, nrow = 1),
          plot_grid(g6pd_boxplot, NULL, labels = letters[4],
                    rel_widths = c(2, 1), ncol = 2, nrow = 1),
          ncol = 1, nrow = 2)

ggsave(filename = here("R/figures/suppl_fig_nparc_fbp1_combo.pdf"), 
              width = 21, height = 14, units = "cm")

# check for correlation with FBP1 AUC
fbp1_1_auc_df <- auc_full_hq_df %>% 
    filter(grepl("^FBP1_1", id)) %>% 
    filter(!sample %in% c("697","COG_355", "COG_394", "MHH_CALL_2"))
