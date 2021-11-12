library(data.table)
library(tidyverse)
library(cowplot)

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

sourceDir <- function(path, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        source(file.path(path, nm), ...)
    }
}
sourceDir(path = here("R/pepnet"))

# define ppi coaggregation analysis output folder
proteoform_detection_folder <- here("proteoform_detection/output/standard")

# load proteoform ratios
proteoforms <- readRDS(
    file.path(proteoform_detection_folder, "proteoforms_narrow_range_focused.RDS"))

# create proteoform data frame
proteoform_df <- biobroom::tidy.ExpressionSet(
    proteoforms, 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

# load psms 
psms <- import_psms_from_nf(file = "/Users/kurzawa/Downloads/all_elac2_ptm_psmtable.txt",
                            id_col = "Gene Name",
                            protein_id_col = "Protein",
                            peptide_col = "Peptide",
                            quan_regex = "tmt16plex")

fData(psms)$gene <- paste(fData(psms)$first_id, fData(psms)$peptide, sep = "_")

psms_df <- 
    bind_cols(fData(psms) %>% 
                  as_tibble() %>% 
                  dplyr::select(gene, set, matches("tmt16plex")),
              exprs(psms) %>% 
                  as_tibble()) %>% 
    gather(tmt_channel, value, matches("tmt16plex")) %>% 
    unite(sample_id, set, tmt_channel, sep = "_")

sample_meta_df = read_tsv(here("meta", "meltome_sample_meta.txt")) %>% 
    dplyr::select(sample_id, sample_name_machine, temperature)

psms_anno_df <- left_join(
    psms_df, sample_meta_df, by = "sample_id") %>% 
    group_by(sample_name_machine, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

filter(psms_anno_df) %>% 
    ggplot(aes(temperature, rel_value)) + 
    geom_point(aes(color = gene == "ELAC2_HQPWQSPERPLSR_Phospho:S6")) +
    geom_smooth(aes(color = gene == "ELAC2_HQPWQSPERPLSR_Phospho:S6"), method = "lm", formula = y ~ splines::ns(x, 4),
                se = FALSE) +
    facet_wrap(~sample_name_machine) +
    theme(legend.position = "bottom")

# make ELAC2 pS199 comparison plot
x_label <- expression("Temperature"* " " * "("*degree*C*")")
y_label <- "Fraction non-denatured"

sem <- function(x, na.rm = FALSE){
    if(na.rm){
        x <- x[!is.na(x)]
    }
    sem_x <- sd(x) / sqrt(length(x))
    return(sem_x)
}

phospho_peptide_plot <- psms_anno_df %>% 
    mutate(group = case_when(gene == "ELAC2_HQPWQSPERPLSR_Phospho:S6" ~ "ELAC2_pS199_peptide",
                             TRUE ~ "all other ELAC2 phospho peptides")) %>% 
    group_by(group, temperature) %>% 
    summarize(mean_rel_value = mean(rel_value, na.rm = TRUE), 
              sem_rel_value = sem(rel_value, na.rm = TRUE)) %>% 
    ungroup %>% 
    ggplot(aes(temperature, mean_rel_value)) + 
    geom_point(aes(color = group)) +
    geom_errorbar(aes(ymin = mean_rel_value - sem_rel_value, 
                      ymax = mean_rel_value + sem_rel_value,
                      color = group), 
                  width=.1, size = 0.5) +
    geom_smooth(aes(color = group), method = "lm", formula = y ~ splines::ns(x, 4),
                se = FALSE) +
    coord_cartesian(y = c(0, 1.25)) +
    scale_color_manual("Phosphopeptides", values = c("all other ELAC2 phospho peptides" = "plum2", "ELAC2_pS199_peptide" = "purple")) +
    #facet_wrap(~sample_name_machine) +
    labs(x = x_label, y = y_label) + 
    theme_paper +
    theme(legend.position = "bottom")

proteoform_plot <- proteoform_df %>% 
    filter(grepl("^ELAC2_", gene)) %>% 
    group_by(gene, temperature) %>% 
    summarize(mean_rel_value = mean(rel_value, na.rm = TRUE), 
              sem_rel_value = sem(rel_value, na.rm = TRUE)) %>% 
    ungroup %>% 
    ggplot(aes(temperature, mean_rel_value)) + 
    geom_point(aes(color = gene)) +
    geom_errorbar(aes(ymin = mean_rel_value - sem_rel_value, 
                      ymax = mean_rel_value + sem_rel_value,
                      color = gene), 
                  width=.1, size = 0.5) +
    geom_smooth(aes(color = gene), method = "lm", formula = y ~ splines::ns(x, 4),
                se = FALSE) +
    coord_cartesian(y = c(0, 1.25)) +
    scale_color_manual("Proteoform", values = c("#FF5376", "#72AFD9")) +
    #facet_wrap(~sample_name_machine) +
    labs(x = x_label, y = y_label) + 
    theme_paper +
    theme(legend.position = "bottom")

plot_grid(phospho_peptide_plot, proteoform_plot, ncol = 2)
ggsave(here("R/figures/figure_elac2_phospho_peptide_suppl.pdf")), 
            width = 12.5, height = 6, unit = "cm")
