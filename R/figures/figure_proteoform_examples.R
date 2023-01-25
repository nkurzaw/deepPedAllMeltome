library(here)
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

# define ppi coaggregation analysis output folder
proteoform_detection_folder <- here("proteoform_detection/output/standard")

# read in peptide data
peptides <- readRDS(file.path(proteoform_detection_folder, "peptides.RDS"))

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

# define colors
cl_colors <-  c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455")
# Peter Carl - R-Bloggers 02/2013 "The Paul Tol 21-color salute"
names(cl_colors) <- sort(
    unique(gsub("_BR2", "", proteoform_df$sample_name_machine)))

# get TMPO proteoform data frame
tmpo_proteoform_df <- biobroom::tidy.ExpressionSet(
    proteoforms[grep("TMPO_", featureNames(proteoforms)),], 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

# get TMPO proteoform data frame
tmpo_peptides_df <- biobroom::tidy.ExpressionSet(
    peptides[grep("TMPO", featureData(peptides)$id),], 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

# median protein TMPO RCH-ACV melting curve
tmpo_RCH_ACV_median_df <- tmpo_peptides_df %>% 
    filter(sample_name_machine == "RCH_ACV") %>% 
    group_by(temperature) %>% 
    summarize(median_rel_value = median(rel_value, na.rm = TRUE)) %>% 
    ungroup()

# median protein TMPO RCH-ACV melting curve
tmpo_COG_355_median_df <- tmpo_peptides_df %>% 
    filter(sample_name_machine == "COG_355") %>% 
    group_by(temperature) %>% 
    summarize(median_rel_value = median(rel_value, na.rm = TRUE)) %>% 
    ungroup()


x_label <- expression("Temperature"* " " * "("*degree*C*")")
y_label <- "Fraction non-denatured"

median_peptide_gg <- 
    ggplot(tmpo_COG_355_median_df, aes(temperature, median_rel_value)) +
    geom_point() +
    geom_smooth(method = "nls", se = FALSE, color = "black",
                formula = y ~ (1-a)/(1 + exp(-(b/x - c))) + a,
                method.args = list(start = c(a = 0, b = 550, c = 10),
                                   algorithm = 'port'),
                alpha = 0.25, size = 0.5) +
    ggtitle("TMPO median peptide signal, COG-355") +
    labs(x = x_label, y = y_label) +
    theme_paper +
    coord_cartesian(ylim = c(0, 1.2))

# ggsave

proteoform_gg <- 
    ggplot(filter(tmpo_proteoform_df, sample_name_machine == "COG_355"), 
           aes(temperature, rel_value, group = gene, color = gene)) +
    geom_point() +
    geom_smooth(method = "nls", se = FALSE,
                formula = y ~ (1-a)/(1 + exp(-(b/x - c))) + a,
                method.args = list(start = c(a = 0, b = 550, c = 10),
                                   algorithm = 'port'),
                alpha = 0.25, size = 0.5) +
    ggtitle("TMPO proteoform signal, COG-355") +
    scale_color_manual("Proteoform", values = c("#FF5376", "#72AFD9")) +
    labs(x = x_label, y = y_label) +
    theme_paper +
    theme(legend.position = "bottom") +
    coord_cartesian(ylim = c(0, 1.2))

proteoform_legend <- get_legend(proteoform_gg)

plot_grid(median_peptide_gg, 
          proteoform_gg + theme(legend.position = "none"),
          NULL, proteoform_legend,
          ncol = 2, rel_heights = c(.9, .1),
          align = "v")

ggsave(filename = "R/figures/tmpo_cog_355_proteoform_vs_average.pdf", 
       width = 12.5, height = 6, units = "cm")

# get ELAC2 proteoform data frame
elac2_proteoform_df <- biobroom::tidy.ExpressionSet(
    proteoforms[grep("ELAC2_", featureNames(proteoforms)),], 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

# get DELAC2peptides data frame
elac2_peptides_df <- biobroom::tidy.ExpressionSet(
    peptides[grep("ELAC2", featureData(peptides)$id),], 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

# median protein DLGAP4 MHH_CALL_3 melting curve
elac2_NALL_1_median_df <- elac2_peptides_df %>% 
    filter(sample_name_machine == "NALL_1") %>% 
    group_by(temperature) %>% 
    summarize(median_rel_value = median(rel_value, na.rm = TRUE)) %>% 
    ungroup()

elac2_median_peptide_gg <- 
    ggplot(elac2_NALL_1_median_df, aes(temperature, median_rel_value)) +
    geom_point() +
    geom_smooth(method = "nls", se = FALSE, color = "black",
                formula = y ~ (1-a)/(1 + exp(-(b/x - c))) + a,
                method.args = list(start = c(a = 0, b = 550, c = 10),
                                   algorithm = 'port'),
                alpha = 0.25, size = 0.5) +
    ggtitle("ELAC2 median peptide signal, NALL-1") +
    labs(x = x_label, y = y_label) +
    theme_paper +
    coord_cartesian(ylim = c(0, 1.2))

elac2_proteoform_gg <- 
    ggplot(filter(elac2_proteoform_df, sample_name_machine == "NALL_1"), 
           aes(temperature, rel_value, group = gene, color = gene)) +
    geom_point() +
    # geom_smooth(method = "nls", se = FALSE,
    #             formula = y ~ (1-a)/(1 + exp(-(b/x - c))) + a,
    #             method.args = list(start = c(a = 0, b = 550, c = 10),
    #                                algorithm = 'port'),
    #             alpha = 0.25, size = 0.5) +
    geom_smooth(method = "lm",
                formula = y ~ splines::ns(x, df = 4),
                se = FALSE, size = 0.5) +
    ggtitle("EALC2 proteoform signal, NALL-1") +
    scale_color_manual("Proteoform", values = c("#FF5376", "#72AFD9")) +
    labs(x = x_label, y = y_label) +
    theme_paper +
    theme(legend.position = "bottom") +
    coord_cartesian(ylim = c(0, 1.2))

elac2_proteoform_legend <- get_legend(elac2_proteoform_gg)

plot_grid(elac2_median_peptide_gg, 
          elac2_proteoform_gg + theme(legend.position = "none"),
          NULL, elac2_proteoform_legend,
          ncol = 2, rel_heights = c(.9, .1),
          align = "v")

ggsave(filename = "R/figures/elac2_nall_1_proteoform_vs_average.pdf", 
       width = 12.5, height = 6, units = "cm")

# define melting curve fit
fitMeltcurve <- function(df, curve_func = drc::LL.4(fixed = c(NA, NA, 1, NA))){
    h1_model <- try(drc::drm(df$rel_value ~ df$temperature,
                             data = df, fct = curve_func),
                    silent = TRUE)
    return(h1_model)
}

# define phosphoTPP plotting function
plotDiffPhosphoMeltcurve <- function(g_name, p_seq, 
                                     predict_vec = seq(35, 70, by = 0.1),
                                     nbf_df, phospho_df, Gene_pSite,
                                     in_ylim = c(0, 1.25)){
    prot_df <- bind_rows(
        filter(phospho_df, gene_name == g_name,
               grepl(p_seq, mod_sequence)) %>% 
            mutate(group = "phosphopeptide"),
        filter(nbf_df, gene_name == g_name)%>% 
            mutate(group = "all unmodified")) 
    
    prot_summarized_df <- prot_df %>% 
        group_by(group, temperature) %>% 
        summarize(mean_rel_value = mean(rel_value, na.rm = TRUE),
                  #sd_rel_value  = sd(rel_value, na.rm = TRUE),
                  se_rel_value  = sd(rel_value, na.rm = TRUE)/
                      sqrt(n())) %>% 
        ungroup()
    
    prot_nbf_fit <- fitMeltcurve(df = filter(prot_df, group == "all unmodified"))
    prot_phospho_fit <- fitMeltcurve(df = filter(prot_df, group == "phosphopeptide"))
    
    prot_fit_df <- tibble(
        temperature = rep(predict_vec, 2),
        mean_rel_value = c(
            predict(prot_nbf_fit, newdata = data.frame(
                temperature = predict_vec)),
            predict(prot_phospho_fit, newdata = data.frame(
                temperature = predict_vec))
        ),
        group = c(
            rep("all unmodified", length(predict_vec)),
            rep("phosphopeptide", length(predict_vec)))
    )
    
    ggplot(prot_summarized_df, aes(temperature, mean_rel_value, color = group)) +
        geom_point(shape = 16) +
        geom_errorbar(aes(ymin = mean_rel_value - se_rel_value, 
                          ymax = mean_rel_value + se_rel_value), 
                      width=.1, size = 0.5) +
        geom_line(data = prot_fit_df, size = 0.5) +
        scale_color_manual("", values = c("all unmodified" = "chartreuse3", "phosphopeptide" = "purple")) +
        theme_bw() +
        labs(x = expression('Temperature ('*~degree*C*')'),
             y = "fraction non-denatured") +
        theme(legend.position = "bottom") +
        ggtitle(gsub("_", " ", Gene_pSite)) +
        coord_cartesian(ylim = in_ylim)
}

plotDiffPhosphoSpline <- function(g_name, p_seq, spline_df = 4,
                                  predict_vec = seq(35, 70, by = 0.1),
                                  nbf_df, phospho_df, Gene_pSite){
    prot_df <- bind_rows(
        filter(phospho_df, gene_name == g_name,
               grepl(p_seq, mod_sequence)) %>% 
            mutate(group = "phosphopeptide"),
        filter(nbf_df, gene_name == g_name)%>% 
            mutate(group = "all unmodified")) 
    
    prot_summarized_df <- prot_df %>% 
        group_by(group, temperature) %>% 
        summarize(mean_rel_value = mean(rel_value, na.rm = TRUE),
                  #sd_rel_value  = sd(rel_value, na.rm = TRUE),
                  se_rel_value  = sd(rel_value, na.rm = TRUE)/
                      sqrt(n())) %>% 
        ungroup()
    
    ggplot(prot_summarized_df, aes(temperature, mean_rel_value, color = group)) +
        geom_point() +
        geom_errorbar(aes(ymin = mean_rel_value - se_rel_value, 
                          ymax = mean_rel_value + se_rel_value), 
                      width=.1) +
        geom_smooth(method = "lm",
                    formula = y ~ splines::ns(x, df = spline_df),
                    se = FALSE, size = 0.5) +
        scale_color_manual("", values = c("all unmodified" = "chartreuse3", "phosphopeptide" = "purple")) +
        theme_bw() +
        labs(x = expression('Temperature ('*~degree*C*')'),
             y = "fraction non-denatured") +
        theme(legend.position = "bottom") +
        ggtitle(gsub("_", " ", Gene_pSite))
}

# read in phosphoTPP data on ELAC2
elac2_phospho_peptide_tab_filtered <-  read_tsv(
    here("data/elac2_phospho_peptide_tab_filtered.txt"))
elac2_nbf_peptide_tab_filtered <- read_tsv(
    here("data/elac2_nbf_peptide_tab_filtered.txt"))

elac2_phosphoTPP_fit <- plotDiffPhosphoSpline(
    predict_vec = seq(41, 63, by = 0.1),
    g_name = "ELAC2", p_seq = "_HQPWQpSPERPLSR_", 
    nbf_df = elac2_nbf_peptide_tab_filtered,
    phospho_df = elac2_phospho_peptide_tab_filtered,
    Gene_pSite = "ELAC2_pS199") +
    coord_cartesian(xlim = c(40, 63),
                    ylim = c(0, 1.25)) +
    theme_paper +
    theme(legend.position = "bottom")

phosphoTPP_legend <- get_legend(elac2_phosphoTPP_fit)

plot_grid(elac2_phosphoTPP_fit + theme(legend.position = "none"),
          phosphoTPP_legend,
          ncol = 1, rel_heights = c(.9, .1),
          align = "v")

ggsave(filename = "R/figures/elac2_phosphoTPP.pdf", 
       width = 6.5, height = 6, units = "cm")

# Source data
source_data_fig2c_1 <- tmpo_peptides_df %>% 
    filter(sample_name_machine == "COG_355") %>% 
    mutate(gene_name = "TMPO") %>% 
    dplyr::select(gene_name, peptide = gene, channel, temperature, sample_name, rel_value) %>% 
    na.omit()

write_csv(source_data_fig2c_1,
          file = here("R/tables/source_data_figure2c_1.csv"))

source_data_fig2c_2 <- tmpo_proteoform_df %>% 
    filter(sample_name_machine == "COG_355") %>% 
    mutate(gene_name = "TMPO") %>% 
    dplyr::select(gene_name, proteoform_group = gene, channel, temperature, 
                  sample_name = sample_name_machine, rel_value = value)

write_csv(source_data_fig2c_2,
          file = here("R/tables/source_data_figure2c_2.csv"))

suppl_table2 <- read_delim(here("R/tables/suppl_table_2_proteoform_detection.txt"), delim = "\t")

source_data_fig2e <- suppl_table2 %>% 
    filter(gene == "TMPO") %>% 
    dplyr::select(gene, proteoform, proteoform_id, peptides, group)

write_csv(source_data_fig2e,
          file = here("R/tables/source_data_figure2e.csv"))

source_data_fig2g_1 <- elac2_peptides_df %>% 
    filter(sample_name_machine == "NALL_1") %>% 
    mutate(gene_name = "ELAC2") %>% 
    dplyr::select(gene_name, peptide = gene, channel, temperature, sample_name, rel_value) %>% 
    na.omit()

write_csv(source_data_fig2g_1,
          file = here("R/tables/source_data_figure2g_1.csv"))

source_data_fig2g_2 <- elac2_proteoform_df %>% 
    filter(sample_name_machine == "NALL_1") %>% 
    mutate(gene_name = "ELAC2") %>% 
    dplyr::select(gene_name, proteoform_group = gene, channel, temperature, 
                  sample_name = sample_name_machine, rel_value = value)

write_csv(source_data_fig2g_2,
          file = here("R/tables/source_data_figure2g_2.csv"))

source_data_fig2h <- suppl_table2 %>% 
    filter(gene == "ELAC2") %>% 
    dplyr::select(gene, proteoform, proteoform_id, peptides, group)

write_csv(source_data_fig2h,
          file = here("R/tables/source_data_figure2h.csv"))

source_data_fig2i <- bind_rows(
    elac2_phospho_peptide_tab_filtered %>% 
        filter(mod_sequence == "_HQPWQpSPERPLSR_") %>% 
        dplyr::select(gene_name, sequence, mod_sequence, replicate, temperature, rel_value),
    elac2_nbf_peptide_tab_filtered %>% 
        mutate(mod_sequence = "unmodified")  %>% 
        dplyr::select(gene_name, sequence, mod_sequence, replicate, temperature, rel_value)
)

write_csv(source_data_fig2i,
          file = here("R/tables/source_data_figure2i.csv"))
