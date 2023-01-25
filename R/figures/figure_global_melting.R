library(tidyverse)
library(Biobase)
library(cowplot)
library(ggrepel)
library(ggExtra)
library(here)

# Plotting theme
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

# define proteoform detection output folder
proteoform_detection_folder <- here("proteoform_detection/output/standard/")

# read in meta data
sample_meta_raw <- read_tsv(file = here("meta/sample_meta.txt"))

# read in peptide data
peptides <- readRDS(file.path(proteoform_detection_folder, "peptides.RDS"))

# read in proteoform data
proteoforms <- readRDS(file.path(proteoform_detection_folder, "proteoforms.RDS"))

# get all proteoform data frame
proteoform_df <- biobroom::tidy.ExpressionSet(
    proteoforms, 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

# define colors
cl_colors <-  c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455")
# Peter Carl - R-Bloggers 02/2013 "The Paul Tol 21-color salute"
names(cl_colors) <- sort(
    unique(gsub("_BR2", "", proteoform_df$sample_name_machine)))

# get all peptides data frame
all_peptides_df <- biobroom::tidy.ExpressionSet(
    peptides, 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

# get all 697 peptides data frame
all_697_peptides_df <- biobroom::tidy.ExpressionSet(
    peptides[,grep("697", phenoData(peptides)$sample_name_machine)], 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

# get all LC4-1 peptides data frame
all_lc41_peptides_df <- biobroom::tidy.ExpressionSet(
    peptides[,which(phenoData(peptides)$sample_name_machine == "LC4_1")], 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup


# get all Kansumi 9 peptides data frame
all_kans9_peptides_df <- biobroom::tidy.ExpressionSet(
    peptides[,which(phenoData(peptides)$sample_name_machine == "KASUMI_9")], 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

# get all Kansumi 9 peptides data frame
all_allpo_peptides_df <- biobroom::tidy.ExpressionSet(
    peptides[,which(phenoData(peptides)$sample_name_machine == "ALL_PO")], 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

# get all REH peptides data frame
all_reh_peptides_df <- biobroom::tidy.ExpressionSet(
    peptides[,which(phenoData(peptides)$sample_name_machine == "REH")], 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

missing_feat_data <- tibble(
    id = featureData(peptides)$id,
    peptide = featureData(peptides)$peptide
)

all_peptides_df <-  left_join(
    all_peptides_df,
    missing_feat_data,
    by = c("gene" = "peptide")
)

all_697_peptides_df <- left_join(
    all_697_peptides_df,
    missing_feat_data,
    by = c("gene" = "peptide")
)

all_lc41_peptides_df <- left_join(
    all_lc41_peptides_df,
    missing_feat_data,
    by = c("gene" = "peptide")
)

all_kans9_peptides_df <- left_join(
    all_kans9_peptides_df,
    missing_feat_data,
    by = c("gene" = "peptide")
)


all_allpo_peptides_df <- left_join(
    all_allpo_peptides_df,
    missing_feat_data,
    by = c("gene" = "peptide")
)

all_reh_peptides_df <- left_join(
    all_reh_peptides_df,
    missing_feat_data,
    by = c("gene" = "peptide")
)

# median protein 697 melting curve
all_697_median_df <- all_697_peptides_df %>% 
    group_by(id, temperature) %>% 
    summarize(median_rel_value = median(rel_value, na.rm = TRUE)) %>% 
    ungroup()

# median protein LC4_1 melting curve
all_lc41_median_df <- all_lc41_peptides_df %>% 
    group_by(id, temperature) %>% 
    summarize(median_rel_value = median(rel_value, na.rm = TRUE)) %>% 
    ungroup()

# median protein Kansumi 9 melting curve
all_kans9_median_df <- all_kans9_peptides_df %>% 
    group_by(id, temperature) %>% 
    summarize(median_rel_value = median(rel_value, na.rm = TRUE)) %>% 
    ungroup()

# median protein ALL PO melting curve
all_allpo_median_df <- all_allpo_peptides_df %>% 
    group_by(id, temperature) %>% 
    summarize(median_rel_value = median(rel_value, na.rm = TRUE)) %>% 
    ungroup()

# median protein ALL PO melting curve
all_reh_median_df <- all_reh_peptides_df %>% 
    group_by(id, temperature) %>% 
    summarize(median_rel_value = median(rel_value, na.rm = TRUE)) %>% 
    ungroup()

# median protein melting curve
all_cell_line_median_df <- all_peptides_df %>% 
    group_by(sample_name_machine, id, temperature) %>% 
    dplyr::summarize(median_rel_value = median(rel_value, na.rm = TRUE)) %>% 
    ungroup()

## supplementary figure 1a
ggplot(all_cell_line_median_df, aes(as.factor(temperature), median_rel_value)) +
    geom_violin(color = NA, fill = "gray", alpha = 0.5) +
    geom_smooth(method = "lm",
                formula = y ~ splines::bs(x, df = 4)) +
    facet_wrap(~sample_name_machine) +
    coord_cartesian(ylim = c(0, 1.5)) +
    theme_paper

## source data supplementary figure 1a
source_data_suppl_fig1a <- all_cell_line_median_df %>% 
    filter(!grepl("BR", sample_name_machine),
           !sample_name_machine %in% c("697", "KASUMI_9", "REH")) %>% 
    dplyr::select(sample = sample_name_machine, gene = id, 
                  temperature, rel_value = median_rel_value) %>% 
    na.omit()

write_csv(source_data_suppl_fig1a, file = here("R/tables/source_data_suppl_fig_1a.csv"))


# median protein melting curve
all_median_df <- all_peptides_df %>% 
    group_by(id, temperature) %>% 
    summarize(median_rel_value = median(rel_value, na.rm = TRUE)) %>% 
    ungroup() 

all_median_rank_df <- all_median_df %>% 
    filter(temperature == 53) %>% 
    mutate(rank = dense_rank(median_rel_value)) %>% 
    dplyr::select(-temperature, -median_rel_value) 

all_median_df <- left_join(
    all_median_df, all_median_rank_df, 
    by = "id")

# define axis labels
x_label <- expression("Temperature"* " " * "("*degree*C*")")
y_label <- "Fraction non-denatured"

ggplot(bind_rows(all_697_median_df,
                 tibble(id = "pseudo", temperature = 70, median_rel_value = 0)), 
       aes(temperature, median_rel_value)) +
    stat_density2d(aes(alpha=..level..), size = 2, 
                   fill = cl_colors[["697"]],#"chartreuse4",
                   bins = 10, geom = "polygon") + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle("697") +
    theme_paper +
    labs(x = x_label, y = y_label) +
    theme(legend.position = "bottom")

ggplot(all_697_median_df, 
       aes(temperature, median_rel_value)) +
    geom_pointdensity() +
    geom_smooth(method = "lm", formula = y ~ splines::bs(x, 4), se = TRUE,
                color = "gray60", alpha = 0.8) +
    coord_cartesian(ylim = c(0, 2.5)) +
    scale_colour_viridis_b() +
    ggtitle("697") +
    theme_paper +
    labs(x = x_label, y = y_label) 
#theme(legend.position = "bottom")

ggplot(all_697_median_df, 
       aes(as.factor(temperature), median_rel_value)) +
    geom_violin(color = NA, fill = "gray60", alpha = 0.5) +
    #geom_boxplot(width = 0.05, alpha = 0.75) +
    geom_smooth(aes(group = 1),
                method = "lm", formula = y ~ splines::bs(x, 4), se = TRUE,
                color = cl_colors[["697"]], alpha = 0.8) +
    coord_cartesian(ylim = c(0, 1.5)) +
    ggtitle("697") +
    theme_paper +
    labs(x = x_label, y = y_label) 

ggsave(filename = "R/figures/meltome_overview_697.pdf", width = 7, height = 4, units = "cm")

ggplot(all_lc41_median_df, 
       aes(as.factor(temperature), median_rel_value)) +
    geom_violin(color = NA, fill = "gray60", alpha = 0.5) +
    #geom_boxplot(width = 0.05, alpha = 0.75) +
    geom_smooth(aes(group = 1),
                method = "lm", formula = y ~ splines::bs(x, 4), se = TRUE,
                color = cl_colors[["LC4_1"]], alpha = 0.8) +
    coord_cartesian(ylim = c(0, 1.5)) +
    ggtitle("LC4-1") +
    theme_paper +
    labs(x = x_label, y = y_label) 

ggsave(filename = "R/figures/meltome_overview_lc41.pdf", width = 7, height = 4, units = "cm")

ggplot(all_kans9_median_df, 
       aes(as.factor(temperature), median_rel_value)) +
    geom_violin(color = NA, fill = "gray60", alpha = 0.5) +
    #geom_boxplot(width = 0.05, alpha = 0.75) +
    geom_smooth(aes(group = 1),
                method = "lm", formula = y ~ splines::bs(x, 4), se = TRUE,
                color = cl_colors[["KASUMI_9"]], alpha = 0.8) +
    coord_cartesian(ylim = c(0, 1.5)) +
    ggtitle("KASUMI-9") +
    theme_paper +
    labs(x = x_label, y = y_label) 

ggsave(filename = "R/figures/meltome_overview_kas9.pdf", width = 7, height = 4, units = "cm")

ggplot(all_allpo_median_df, 
       aes(as.factor(temperature), median_rel_value)) +
    geom_violin(color = NA, fill = "gray60", alpha = 0.5) +
    geom_boxplot(width = 0.05, alpha = 0.75) +
    geom_smooth(aes(group = 1),
                method = "lm", formula = y ~ splines::bs(x, 4), se = TRUE,
                color = cl_colors[["ALL_PO"]], alpha = 0.8) +
    coord_cartesian(ylim = c(0, 1.5)) +
    ggtitle("ALL-PO") +
    theme_paper +
    labs(x = x_label, y = y_label) 

ggsave(filename = "R/figures/meltome_overview_allpo.pdf", width = 7, height = 4, units = "cm")


ggplot(all_reh_median_df, 
       aes(as.factor(temperature), median_rel_value)) +
    geom_violin(color = NA, fill = "gray60", alpha = 0.5) +
    #geom_boxplot(width = 0.05, alpha = 0.75) +
    geom_smooth(aes(group = 1),
                method = "lm", formula = y ~ splines::bs(x, 4), se = TRUE,
                color = cl_colors[["REH"]], alpha = 0.8) +
    coord_cartesian(ylim = c(0, 1.5)) +
    ggtitle("REH") +
    theme_paper +
    labs(x = x_label, y = y_label) 

ggsave(filename = "R/figures/meltome_overview_reh.pdf", width = 7, height = 4, units = "cm")

source_data_fig1b <- bind_rows(
    all_697_median_df %>% mutate(sample = "697"),
    all_kans9_median_df %>% mutate(sample = "KASUMI-9"),
    all_reh_median_df %>%  mutate(sample = "REH")
)

write_csv(source_data_fig1b,
          file = here("R/tables/source_data_figure1b.csv"))
