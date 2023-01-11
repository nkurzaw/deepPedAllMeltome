library(tidyverse)
library(Biobase)
library(ggthemes)
library(NPARC)
library(BiocParallel)
library(cowplot)
library(ggrepel)
library(ggExtra)
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

# source all files in the functions directory
sourceDir <- function(path, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        source(file.path(path, nm), ...)
    }
}
sourceDir(path = file.path(here("R/nparc")))

# get sample meta data
sample_meta_raw <- read_tsv(here("meta/sample_meta.txt"))

# get sample meta data on b-cell development
bcell_maturation_anno <- read_tsv(here("meta/sample_meta_stages.txt"))

# define output folder for figures
figure_output_folder <- here("R/figures")

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

# load nparc fstat data frame
nparc_fstat_df <- readRDS(here("nparc/output/standard/nparc_fstat_df.RDS"))

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

# NEK2 example plots
nek2_0_df <- filter(proteoform_df, gene == "NEK2_0") %>% 
    filter(!grepl("_BR2", sample_name_machine)) %>% 
    filter(sample_name_machine %in% filter(nparc_res_hq_df, id == "NEK2_0")$sample_name) %>% 
    na.omit()

control <-  NPARC:::getParams()
nek2_0_null_fit_param <- NPARC:::invokeParallelFits(
    x = nek2_0_df$temperature, 
    y = nek2_0_df$rel_value, 
    id = nek2_0_df$gene, 
    groups = NULL,
    BPPARAM = SerialParam(progressbar = TRUE),
    maxAttempts = control$maxAttempts,
    returnModels = FALSE,
    start = control$start)


temp_range <- seq(from = 40, to = 65, by = 0.1)

nek2_0_null_fit_df <- 
    tibble(temperature = rep(temp_range, length(unique(nek2_0_df$sample_name_machine))),
           group = rep(
               unique(nek2_0_df$sample_name_machine),
               each = length(temp_range)),
           id = "NEK2_0") %>% 
    left_join(nek2_0_null_fit_param$modelMetrics, 
              by = "id") %>% 
    rowwise() %>% 
    mutate(y_hat = (1 - pl)  / (1 + exp((b - a/temperature))) + pl) %>% 
    ungroup

nek2_0_alt_fit_param <- NPARC:::invokeParallelFits(
    x = nek2_0_df$temperature, 
    y = nek2_0_df$rel_value, 
    id = nek2_0_df$gene, 
    groups = nek2_0_df$sample_name_machine,
    BPPARAM = SerialParam(progressbar = TRUE),
    maxAttempts = control$maxAttempts,
    returnModels = FALSE,
    start = control$start)

nek2_0_alt_fit_df <- 
    tibble(temperature = rep(temp_range, length(unique(nek2_0_df$sample_name_machine))),
           group = rep(
               unique(nek2_0_df$sample_name_machine),
               each = length(temp_range)
           )) %>% 
    left_join(nek2_0_alt_fit_param$modelMetrics, 
              by = "group") %>% 
    rowwise() %>% 
    mutate(y_hat = (1 - pl)  / (1 + exp((b - a/temperature))) + pl) %>% 
    ungroup

nek2_0_rss_0 <- nek2_0_null_fit_df$rss[1]
nek2_0_rss_1 <- (nek2_0_alt_fit_param$modelMetrics %>% 
                     dplyr::summarize(rss = sum(rss)))$rss

nek2_0_null_plot <- 
    ggplot(nek2_0_null_fit_param$modelPredictions %>% 
               mutate(group = nek2_0_df$sample_name_machine), 
           aes(x, y)) +
    geom_line(aes(temperature, y_hat), 
              data = nek2_0_null_fit_df,
              color = "darkgray") +
    geom_point(aes(color = group)) +
    geom_segment(aes(xend = x, yend = .fitted), 
                 linetype = "dashed") +
    # geom_text(data = tibble(),
    #           aes(x = 60, y = 0.95),
    #               label = as.expression(
    #                 bquote('RSS'^0 == .(round(tp53_2_rss_0, 3))))) +
    scale_color_manual("", values = cl_colors) +
    scale_x_continuous(breaks = c(40, 55)) +
    labs(x = x_label,
         y = y_label) +
    facet_wrap(~ group, ncol = 6) +
    ggtitle(bquote('NEK2 null model; RSS'^0 == .(round(nek2_0_rss_0, 4)))) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    theme_paper

nek2_0_alt_plot <- 
    ggplot(nek2_0_alt_fit_param$modelPredictions , 
           aes(x, y)) +
    geom_line(aes(temperature, y_hat, group = group, color = group), 
              data = nek2_0_alt_fit_df)+
    geom_point(aes(color = group)) +
    geom_segment(aes(xend = x, yend = .fitted), 
                 linetype = "dashed") +
    # geom_text(data = tibble(),
    #           aes(x = 60, y = 0.95),
    #                 label = as.expression(
    #                   bquote('RSS'^1 == .(round(abl1_1_rss_1, 3))))) +
    scale_color_manual("", values = cl_colors) +
    scale_x_continuous(breaks = c(40, 55)) +
    labs(x = x_label,
         y = y_label) +
    facet_wrap(~ group, ncol = 6) +
    ggtitle(bquote('NEK2 alternative model; RSS'^1 == .(round(nek2_0_rss_1, 4)))) +
    theme(legend.position = "bottom",
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    theme_paper

nek2_0_legend <- get_legend(nek2_0_alt_plot)

nparc_fit_p <- 
    plot_grid(nek2_0_null_plot + theme(legend.position = "none"), 
              nek2_0_alt_plot + theme(legend.position = "none"),
              ncol = 1, rel_heights = c(45, 45))

ggsave(nparc_fit_p, filename = here("R/figures/figure_nparc_fit_schematic_nek2.pdf"), 
       width = 21, height = 14, units = "cm")

fstat_90_quan <- quantile(nparc_fstat_df$F_statistic, 0.9)

# inset <- ggplot(nparc_fstat_df, aes(F_statistic)) +
#     geom_density(fill = "black", alpha = 0.7) +
#     geom_rug(alpha = 0.7) +
#     geom_text(
#         data = tibble(x = 100, y = 0.065),
#         label = expression(''*italic(F)*'' == frac(('RSS'^0 - 'RSS'^1) * d[2], ('RSS'^1) * d[1])),
#         aes(x = x, y = y),
#         size = 2.5) +
#     # geom_label_repel(
#     #   y = 0,
#     #   label = "TP53_2",  
#     #   nudge_y = 0.02,
#     #   direction = "x",
#     #   segment.size = 0.25,
#     #   color = "cadetblue3", 
#     #   data = filter(nparc_fstat_df,id == "TP53_2")) +
#     xlab(expression(''*italic(F)*'-statistic')) +
#     theme_paper

nparc_volcano <- 
    ggplot(nparc_fstat_df, aes(rss-rssAlternative, F_statistic)) +
    geom_point(alpha = 0.25, color = "gray", size = 0.5) +
    geom_point(color = "black", alpha = 0.5, size = 0.5,
               data = filter(nparc_fstat_df, F_statistic >= fstat_90_quan)) +
    geom_label_repel(
        aes(label = id),  
        nudge_x = 1,
        direction = "x",
        segment.size = 0.25,
        color = "darkgray", 
        data = filter(nparc_fstat_df,id %in% 
                          c("NEK2_0", "DNTT_1", "INPP4B_1", "TP53_1", "FBP1_1"))) +
    scale_x_log10() +
    labs(x = bquote('RSS'^0~' - RSS'^1~''),
         y = expression(''*italic(F)*'-statistic')) +
    coord_cartesian(xlim = c(0.05, 17.5)) +
    theme_paper


# nparc_volcano_p <- ggdraw(nparc_volcano) +
#     draw_plot(inset, .15, .6, .5, .35) 

ggsave(nparc_volcano, 
       filename = here("R/figures/figure_nparc_volcano_all_highlighted_black_gray.pdf"), 
       width = 7, height = 8, units = "cm")

# DNTT1 
dntt1_1_df <- filter(proteoform_df, gene == "DNTT_1") %>% 
    filter(!grepl("_BR2", sample_name_machine)) %>% 
    filter(sample_name_machine %in% filter(nparc_res_hq_df, id == "DNTT_1")$sample_name) %>% 
    na.omit()

control <-  NPARC:::getParams()
temp_range <- seq(from = 40, to = 65, by = 0.1)

dntt1_1_alt_fit_param <- NPARC:::invokeParallelFits(
    x = dntt1_1_df$temperature, 
    y = dntt1_1_df$rel_value, 
    id = dntt1_1_df$gene, 
    groups = dntt1_1_df$sample_name_machine,
    BPPARAM = SerialParam(progressbar = TRUE),
    maxAttempts = control$maxAttempts,
    returnModels = FALSE,
    start = control$start)

dntt1_1_alt_fit_df <- 
    tibble(temperature = rep(temp_range, 20),
           group = rep(
               unique(dntt1_1_df$sample_name_machine),
               each = length(temp_range)
           )) %>% 
    left_join(dntt1_1_alt_fit_param$modelMetrics, 
              by = "group") %>% 
    rowwise() %>% 
    mutate(y_hat = (1 - pl)  / (1 + exp((b - a/temperature))) + pl) %>% 
    ungroup

dntt1_1_plot_df <- dntt1_1_alt_fit_param$modelPredictions %>% 
    left_join(sample_meta_raw %>% 
                  dplyr::select(group = sample_name_machine, subtype),
              by = "group") %>% 
    left_join(bcell_maturation_anno %>% 
                  dplyr::select(group = sample, Assigned_stage),
              by = "group") %>% 
    mutate(simple_stage = factor(
        case_when(grepl("other", Assigned_stage) ~ "other",
                  grepl("pre-pro", Assigned_stage) ~ "pre-pro-B",
                  grepl("pro", Assigned_stage) ~ "pro-B",
                  grepl("pre", Assigned_stage) ~ "pre-B"),
        levels = c("pre-pro-B", "pro-B", "pre-B", "other")))

dntt1_1_alt_fit_df <- dntt1_1_alt_fit_df %>% 
    left_join(bcell_maturation_anno %>% 
                  dplyr::select(group = sample, Assigned_stage),
              by = "group") %>% 
    mutate(simple_stage = factor(
        case_when(grepl("other", Assigned_stage) ~ "other",
                  grepl("pre-pro", Assigned_stage) ~ "pre-pro-B",
                  grepl("pro", Assigned_stage) ~ "pro-B",
                  grepl("pre", Assigned_stage) ~ "pre-B"),
        levels = c("pre-pro-B", "pro-B", "pre-B", "other")))

dntt1_p <- ggplot(dntt1_1_plot_df, 
                  aes(x, y)) +
    geom_line(aes(temperature, y_hat, color = simple_stage, group = group), #, linetype = Assigned_stage), 
              data = dntt1_1_alt_fit_df)+#,
    #color = "darkgray") +
    geom_point(aes(color = simple_stage, group = group)) +
    # geom_segment(aes(xend = x, yend = .fitted), 
    #              linetype = "dashed") +
    #scale_color_manual("", values = cl_colors) +
    scale_color_manual("Developmental\nstage", 
                       values = c("pre-pro-B" = "darkslateblue", "pre-B" = "purple", "pro-B" = "orange", "other" = "gray")) +
    #scale_x_continuous(breaks = c(40, 55)) +
    labs(x = x_label,
         y = y_label) +
    #facet_wrap(~ subtype, ncol = 11) +
    ggtitle("DNTT proteoform 1") +
    theme_paper #+
#theme(legend.position = "bottom")

ggsave(dntt1_p,
       filename = here("R/figures/figure_dntt_1_stage_example_color_optimized.pdf"), 
       width = 10, height = 6, units = "cm")

### INPP4B
inpp4b_1_df <- filter(proteoform_df, gene == "INPP4B_1") %>% 
    filter(!grepl("_BR2", sample_name_machine)) %>% 
    filter(sample_name_machine %in% filter(nparc_res_hq_df, id == "INPP4B_1")$sample_name) %>% 
    na.omit()

control <-  NPARC:::getParams()
temp_range <- seq(from = 40, to = 65, by = 0.1)

inpp4b_1_alt_fit_param <- NPARC:::invokeParallelFits(
    x = inpp4b_1_df$temperature, 
    y = inpp4b_1_df$rel_value, 
    id = inpp4b_1_df$gene, 
    groups = inpp4b_1_df$sample_name_machine,
    BPPARAM = SerialParam(progressbar = TRUE),
    maxAttempts = control$maxAttempts,
    returnModels = FALSE,
    start = control$start)

inpp4b_1_alt_fit_df <- 
    tibble(temperature = rep(temp_range, length(unique(inpp4b_1_df$sample_name_machine))),
           group = rep(
               unique(inpp4b_1_df$sample_name_machine),
               each = length(temp_range)
           )) %>% 
    left_join(inpp4b_1_alt_fit_param$modelMetrics, 
              by = "group") %>% 
    rowwise() %>% 
    mutate(y_hat = (1 - pl)  / (1 + exp((b - a/temperature))) + pl) %>% 
    ungroup

inpp4b_1_plot_df <- inpp4b_1_alt_fit_param$modelPredictions %>% 
    left_join(sample_meta_raw %>% 
                  dplyr::select(group = sample_name_machine, subtype),
              by = "group") %>% 
    left_join(bcell_maturation_anno %>% 
                  dplyr::select(group = sample, Assigned_stage),
              by = "group") %>% 
    mutate(subtype_tcfpbx = factor(
        case_when(subtype == "TCF3.PBX1" ~ "TCF3-PBX1",
                  TRUE ~ "other"),
        levels = c("TCF3-PBX1", "other")))

inpp4b_1_alt_fit_df <- inpp4b_1_alt_fit_df %>% 
    left_join(bcell_maturation_anno %>% 
                  dplyr::select(group = sample, Assigned_stage, subtype),
              by = "group") %>% 
    mutate(subtype_tcfpbx = factor(
        case_when(subtype == "TCF3.PBX1" ~ "TCF3-PBX1",
                  TRUE ~ "other"),
        levels = c("TCF3-PBX1", "other")))

inpp4b_p <- ggplot(inpp4b_1_plot_df, 
                   aes(x, y)) +
    geom_line(aes(temperature, y_hat, color = subtype_tcfpbx, group = group), #, linetype = Assigned_stage), 
              data = inpp4b_1_alt_fit_df)+#,
    #color = "darkgray") +
    geom_point(aes(color = subtype_tcfpbx)) +
    # geom_segment(aes(xend = x, yend = .fitted), 
    #              linetype = "dashed") +
    #scale_color_manual("", values = cl_colors) +
    scale_color_manual("Subtype", values = c("TCF3-PBX1" = "deeppink", "other" = "gray")) +
    #scale_x_continuous(breaks = c(40, 55)) +
    labs(x = x_label,
         y = y_label) +
    #facet_wrap(~ subtype, ncol = 11) +
    ggtitle("INPP4B proteoform 1") +
    theme_paper


ggsave(inpp4b_p, filename = here("R/figures/figure_inpp4b_1_stage_example_opt_colors.pdf"), 
       width = 10, height = 6, units = "cm")


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

ggsave(fbp1_p, 
       filename = here("R/figures/figure_nparc_fbp1_1_melting_profile.pdf"), 
       width = 7, height = 7, units = "cm")

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
    mutate(cell_line = sub("_BR2", "", sample_name_machine))

# join anno with qms data
qms_anno_df <- left_join(qms_df, qms_pdata_df, by = "sample")

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

ggsave(fbp1_qms_bp, 
       filename = here("R/figures/figure_fbp1_qms_boxplot.pdf"), 
       width = 7, height = 7, units = "cm")

# prepare TBC1D10C-specific data set with thermal stability annotation
tbc1_qms_df <- qms_anno_df %>% 
    filter(gene == "TBC1D10C") %>% 
    mutate(eps8l2_thermal_stability_group = case_when(
        sample_name %in% c("TMD5_BR2", "SEM", "NALL-1", "KOPN-8", "LC4-1", 
                           "COG-LL-319h", "REH") ~ "high",
        sample_name %in% c("KASUMI-2", "MHH-CALL-2", "KASUMI-9_BR2",
                           "SUP-B15", "HAL-01") ~ "low",
        TRUE ~ "none")) %>% 
    filter(eps8l2_thermal_stability_group != "none") %>% 
    mutate(eps8l2_thermal_stability_group = factor(
        eps8l2_thermal_stability_group, levels = c("low",
                                            "high"))) 

# make boxplot
tbc1_qms_bp <- ggplot(tbc1_qms_df, aes(eps8l2_thermal_stability_group, value)) +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(aes(color = cell_line), width = 0.1) +
    scale_color_manual("Cell line", values = cl_colors) +
    ggsignif::geom_signif(comparisons = list(c("low", "high")),
                          test = t.test) +
    labs(x = "EPS8L2_2 thermal stability", 
         y = bquote('TBC1D10C log'[2]*'fold change to mean')) +
    #coord_cartesian(ylim = c(-2, 3.5)) +
    theme_paper +
    theme(legend.position = "none")

# Source data
source_data_fig3a_1 <- nek2_0_df %>% 
    dplyr::select(proteoform = gene, channel, temperature, 
                  sample_name = sample_name_machine, rel_value)

write_csv(source_data_fig3a_1, file = here("R/tables/source_data_fig3a_1.csv"))

source_data_fig3a_2 <- nek2_0_null_fit_param$modelMetrics

write_csv(source_data_fig3a_2, file = here("R/tables/source_data_fig3a_2.csv"))

source_data_fig3a_3 <- nek2_0_df %>% 
    dplyr::select(proteoform = gene, channel, temperature, 
                  sample_name = sample_name_machine, rel_value)

write_csv(source_data_fig3a_3, file = here("R/tables/source_data_fig3a_3.csv"))

source_data_fig3a_4 <- nek2_0_alt_fit_param$modelMetrics

write_csv(source_data_fig3a_4, file = here("R/tables/source_data_fig3a_4.csv"))

write_csv(nparc_fstat_df, file = here("R/tables/source_data_fig3b.csv"))

source_data_fig3c <- dntt1_1_df %>% 
    dplyr::select(proteoform = gene, channel, temperature, 
                  sample_name = sample_name_machine, rel_value,
                  state)

write_csv(source_data_fig3c, file = here("R/tables/source_data_fig3c.csv"))

source_data_fig3d <- inpp4b_1_df %>% 
    dplyr::select(proteoform = gene, channel, temperature, 
                  sample_name = sample_name_machine, rel_value,
                  subtype)

write_csv(source_data_fig3d, file = here("R/tables/source_data_fig3d.csv"))
