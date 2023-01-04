library(tidyverse)
library(Biobase)
library(cowplot)
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

# get sample meta data
sample_meta_raw <- read_tsv(here("meta/sample_meta.txt"))

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

# set axis labels
x_label <- expression("Temperature"* " " * "("*degree*C*")")
y_label <- "Fraction non-denatured"

# define color scheme
cl_colors <-  c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
                "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", 
                "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", 
                "#771122", "#AA4455")
# Peter Carl - R-Bloggers 02/2013 "The Paul Tol 21-color salute"
names(cl_colors) <- sort(
    unique(gsub("_BR2", "", proteoform_df$sample_name_machine)))

# get unfiltered NPAC alternative model results
nparc_alt_result_df <- bind_rows(lapply(
    list.files(here("nparc/output/standard/"), pattern = "res_[^(hq)]", full.names = TRUE),
    function(file){
        cell_line <- sub(".RDS", "", sub("nparc_res_", "", basename(file)))
        df <- readRDS(file)
        df$sample <- cell_line
        return(df)
    }))


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


# filter nparc data to only _0 proteoforms
nparc_sd_df <- nparc_alt_result_df %>% 
    filter(grepl("_0", id)) %>% 
    dplyr::select(id, resid_sd, sample) %>% 
    mutate(gene = sub("_0", "", id)) %>% 
    left_join(qms_anno_df %>% 
                  dplyr::select(gene, sample = sample_name_machine, rel_abundance = value), 
              by = c("gene", "sample")) %>% 
    filter(!is.na(resid_sd), !grepl("BR2", sample)) 


ggplot(nparc_sd_df, aes(resid_sd < 0.1, rel_abundance)) +
    geom_violin(fill = "gray", alpha = 0.5, color = NA) +
    geom_boxplot(outlier.colour = NA, width = 0.25) +
    geom_hline(yintercept = 0, alpha = 0.5, color = "orange") +
    facet_wrap(~sample) +
    theme_paper


# examples of profiles with high noise
## AATK_0
AATK_0_df <- filter(proteoform_df, gene == "AATK_0") %>% 
    filter(!grepl("_BR2", sample_name_machine)) %>% 
    filter(sample_name_machine %in% filter(nparc_res_hq_df, id == "AATK_0")$sample_name) %>% 
    na.omit()

control <-  NPARC:::getParams()
AATK_0_null_fit_param <- NPARC:::invokeParallelFits(
    x = AATK_0_df$temperature, 
    y = AATK_0_df$rel_value, 
    id = AATK_0_df$gene, 
    groups = NULL,
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
    maxAttempts = control$maxAttempts,
    returnModels = FALSE,
    start = control$start)

temp_range <- seq(from = 40, to = 65, by = 0.1)

AATK_0_alt_fit_param <- NPARC:::invokeParallelFits(
    x = AATK_0_df$temperature, 
    y = AATK_0_df$rel_value, 
    id = AATK_0_df$gene, 
    groups = AATK_0_df$sample_name_machine,
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
    maxAttempts = control$maxAttempts,
    returnModels = FALSE,
    start = control$start)

AATK_0_alt_fit_df <- 
    tibble(temperature = rep(temp_range, 6),
           group = rep(
               unique(AATK_0_df$sample_name_machine),
               each = length(temp_range)
           )) %>% 
    left_join(AATK_0_alt_fit_param$modelMetrics, 
              by = "group") %>% 
    rowwise() %>% 
    mutate(y_hat = (1 - pl)  / (1 + exp((b - a/temperature))) + pl) %>% 
    ungroup

AATK_0_plot_df <- AATK_0_alt_fit_param$modelPredictions %>% 
    left_join(sample_meta_raw %>% 
                  dplyr::select(group = sample_name_machine, subtype),
              by = "group") 

AATK_0_plot_df <- AATK_0_plot_df %>% 
    mutate(resid_sd_lt_01 = case_when(group %in% c("LC4_1", "RCH_ACV") ~ TRUE,
                                       TRUE ~ FALSE))

ggplot(AATK_0_plot_df, 
       aes(x, y)) +
    geom_rect(data = AATK_0_plot_df, aes(fill = resid_sd_lt_01), xmin = -Inf, xmax = Inf,
              ymin = -Inf, ymax = Inf, alpha = 0.2, show.legend = FALSE) +
    scale_fill_manual(values = c("white", "gray")) +
    geom_line(aes(temperature, y_hat, color = group, group = group), #, linetype = Assigned_stage), 
              data = AATK_0_alt_fit_df)+#,
    #color = "darkgray") +
    geom_point(aes(color = group, group = group)) +
    geom_text(aes(label = as.character(round(resid_sd, 4))), 
              data = filter(nparc_alt_result_df, id == "AATK_0") %>% 
                  mutate(x = 55, y = 1.5, group = sample) %>% 
                  filter(!grepl("BR", group))) +
    # geom_segment(aes(xend = x, yend = .fitted), 
    #              linetype = "dashed") +
    scale_color_manual("", values = cl_colors) +
    #scale_x_continuous(breaks = c(40, 55)) +
    labs(x = x_label,
         y = y_label) +
    facet_wrap(~ group, ncol = 11) +
    ggtitle("AATK_0") +
    theme_paper +
    theme(legend.position = "none")

ggsave("~/Downloads/AATK_0_profile_resid_sd_filter.pdf", width = 18, height = 7.5, units = "cm")

AATK_0_rss_0 <- AATK_0_null_fit_param$modelMetrics$rss
AATK_0_rss_1 <- (AATK_0_alt_fit_param$modelMetrics %>% 
                     dplyr::summarize(rss = sum(rss)))$rss

AATK_0_nCoeffsNull <- AATK_0_null_fit_param$modelMetrics$nCoeffs

AATK_0_nCoeffsAlternative <- (AATK_0_alt_fit_param$modelMetrics %>% 
                                  dplyr::summarize(nCoeffs = sum(nCoeffs)))$nCoeffs

AATK_0_nFittedAlternative <- (AATK_0_alt_fit_param$modelMetrics %>% 
                                  dplyr::summarize(nFitted = sum(nFitted)))$nFitted

AATK_0_F_statistic <- ((AATK_0_rss_0 - AATK_0_rss_1)/AATK_0_rss_1) * 
    ((AATK_0_nFittedAlternative - AATK_0_nCoeffsAlternative)/(AATK_0_nCoeffsAlternative - AATK_0_nCoeffsNull))

## CD96_2
CD96_2_df <- filter(proteoform_df, gene == "CD96_2") %>% 
    filter(!grepl("_BR2", sample_name_machine)) %>% 
    filter(sample_name_machine %in% filter(nparc_res_hq_df, id == "TP53_1")$sample_name) %>% 
    na.omit()

control <-  NPARC:::getParams()

CD96_2_null_fit_param <- NPARC:::invokeParallelFits(
    x = CD96_2_df$temperature, 
    y = CD96_2_df$rel_value, 
    id = CD96_2_df$gene, 
    groups = NULL,
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
    maxAttempts = control$maxAttempts,
    returnModels = FALSE,
    start = control$start)

temp_range <- seq(from = 40, to = 65, by = 0.1)

CD96_2_alt_fit_param <- NPARC:::invokeParallelFits(
    x = CD96_2_df$temperature, 
    y = CD96_2_df$rel_value, 
    id = CD96_2_df$gene, 
    groups = CD96_2_df$sample_name_machine,
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
    maxAttempts = control$maxAttempts,
    returnModels = FALSE,
    start = control$start)

CD96_2_alt_fit_df <- 
    tibble(temperature = rep(temp_range, 14),
           group = rep(
               unique(CD96_2_df$sample_name_machine),
               each = length(temp_range)
           )) %>% 
    left_join(CD96_2_alt_fit_param$modelMetrics, 
              by = "group") %>% 
    rowwise() %>% 
    mutate(y_hat = (1 - pl)  / (1 + exp((b - a/temperature))) + pl) %>% 
    ungroup

CD96_2_plot_df <- CD96_2_alt_fit_param$modelPredictions %>% 
    left_join(sample_meta_raw %>% 
                  dplyr::select(group = sample_name_machine, subtype),
              by = "group") 

CD96_2_plot_df <- CD96_2_plot_df %>% 
    mutate(resid_sd_lt_01 = case_when(group %in% c("MHH_CALL_3") ~ TRUE,
                                      TRUE ~ FALSE))

ggplot(CD96_2_plot_df, 
       aes(x, y)) +
    geom_rect(data = CD96_2_plot_df, aes(fill = resid_sd_lt_01), xmin = -Inf, xmax = Inf,
              ymin = -Inf, ymax = Inf, alpha = 0.2, show.legend = FALSE) +
    scale_fill_manual(values = c("white", "gray")) +
    geom_line(aes(temperature, y_hat, color = group, group = group), #, linetype = Assigned_stage), 
              data = CD96_2_alt_fit_df)+#,
    #color = "darkgray") +
    geom_point(aes(color = group, group = group)) +
    # geom_segment(aes(xend = x, yend = .fitted), 
    #              linetype = "dashed") +
    geom_text(aes(label = as.character(round(resid_sd, 4))), 
              data = filter(nparc_alt_result_df, id == "CD96_2") %>% 
                  mutate(x = 55, y = 0.1, group = sample) %>% 
                  filter(!grepl("BR", group))) +
    scale_color_manual("", values = cl_colors) +
    scale_y_continuous(limits = c(0, 1.3)) +
    labs(x = x_label,
         y = y_label) +
    facet_wrap(~ group, ncol = 7) +
    ggtitle("CD96_2") +
    theme_paper +
    theme(legend.position = "none")

ggsave("~/Downloads/CD96_2_profile_resid_sd_filter.pdf", width = 18, height = 14, units = "cm")

CD96_2_rss_0 <- CD96_2_null_fit_param$modelMetrics$rss
CD96_2_rss_1 <- (CD96_2_alt_fit_param$modelMetrics %>% 
                     dplyr::summarize(rss = sum(rss)))$rss

CD96_2_nCoeffsNull <- CD96_2_null_fit_param$modelMetrics$nCoeffs

CD96_2_nCoeffsAlternative <- (CD96_2_alt_fit_param$modelMetrics %>% 
                                  dplyr::summarize(nCoeffs = sum(nCoeffs)))$nCoeffs

CD96_2_nFittedAlternative <- (CD96_2_alt_fit_param$modelMetrics %>% 
                                  dplyr::summarize(nFitted = sum(nFitted)))$nFitted

CD96_2_F_statistic <- ((CD96_2_rss_0 - CD96_2_rss_1)/CD96_2_rss_1) * 
    ((CD96_2_nFittedAlternative - CD96_2_nCoeffsAlternative)/(CD96_2_nCoeffsAlternative - CD96_2_nCoeffsNull))

# make combined supplementary figure
