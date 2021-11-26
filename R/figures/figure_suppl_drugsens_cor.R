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

sourceDir <- function(path, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
        source(file.path(path, nm), ...)
    }
}
sourceDir(path = file.path(here("R/drugsens_cor")))

# define ppi coaggregation analysis output folder
drugsens_cor_folder <- here("drugsens_cor/output")

# define output folder for figures
figure_output_folder <- here("R/figures")

# load nparc hq results
nparc_res_hq_df <- readRDS(here("nparc/output/standard/nparc_res_hq_df.RDS"))

# load proteoform ratios
proteoforms <- readRDS(
    here("proteoform_detection/output/standard/proteoforms_narrow_range_focused.RDS"))

# create proteoform data frame
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
    unique(gsub("_", "-", sub("_BR2", "", proteoform_df$sample_name_machine))))

# define axis lavels
x_label <- expression("Temperature"* " " * "("*degree*C*")")
y_label <- "Fraction non-denatured"

# get fignl specific proteoform dataset
fignl_proteoform_df <- proteoform_df %>% 
    filter(gene == "FIGNL1_1", !grepl("BR", sample_name_machine)) %>% 
    left_join(nparc_res_hq_df %>% dplyr::select(id, sample_name, conv),
              by = c("gene" = "id", "sample_name_machine" = "sample_name")) %>% 
    filter(conv == TRUE)  %>% 
    mutate(sample_name_machine = gsub("_", "-", sample_name_machine))

# plot FIGNL_1 proteoforms across cell lines
fignl_proteoform_profile_plot <- 
    ggplot(fignl_proteoform_df, 
           aes(temperature, rel_value, color = sample_name_machine)) +
    geom_smooth(aes(group = sample_name_machine), method = "lm",
                formula = 'y ~ splines::ns(x, df = 4)',
                se = FALSE, alpha = 0.5, size = 0.75) +
    # stat_smooth(aes(group = sample_name_machine),
    #             method = "nls", se = FALSE,
    #             formula = y ~ (1-a)/(1 + exp(-(b/x - c))) + a,
    #             method.args = list(start = c(a = 0.1, b = 1550, c = 40),
    #                                algorithm = 'port'),
    #             geom = "line",
    #             alpha = 0.5) +
    scale_color_manual("Subtype", 
                       values = cl_colors) +
    theme_paper +
    theme(legend.position = "none") +
    labs(x = x_label,
         y = y_label) +
    ggtitle("FIGNL proteoform 1")

ggsave(filename = here("R/figures/fig_fignl1_spline_fit.pdf"), 
       width = 7, height = 7, units = "cm")

# read in results from limma analysis
limma_out_df <- readRDS(file.path(drugsens_cor_folder, "limma_out_df.RDS"))

limma_out_anno_df <- limma_out_df %>% 
    mutate(id = 1:n()) 

# loop over significant hits and make scatter plots of thermal stability and drug sens.
pList <- lapply(seq(nrow(filter(limma_out_anno_df, p_adj < 0.1))), function(i) {
    protn <- limma_out_anno_df$testedProt[i] 
    drugn <- limma_out_anno_df$rowname[i]
    p_adj <- limma_out_anno_df$p_adj[i]
    
    title_string <- paste(drugn, protn, sep = " vs. ")
    
    auc_fil_df <- auc_mat_norm %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "id") %>% 
        filter(id == protn) %>% 
        gather(sample, aumc, -id) %>% 
        filter(!grepl("_BR", sample)) 
    
    
    dss_fil_df <- dss_mat %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "sample") %>% 
        gather(drug, dss, -sample) %>% 
        filter(drug == drugn) 
    
    plot_df <- left_join(
        auc_fil_df,
        dss_fil_df,
        by = "sample") 
    
    ggplot(plot_df, aes(aumc, dss)) +
        geom_point(aes(color = sample)) +
        geom_smooth(method = "lm", color = "black") +
        scale_color_manual(values = cl_colors) +
        labs(x = "Area under the melting curve",
             y = "Drug sensitivity") +
        ggpubr::stat_cor(method = "pearson",
                         label.x.npc = "left",
                         cor.coef.name = "rho") +
        # geom_text(aes(x, y), label = round(p_adj, 3), 
        #           data = tibble(x = 0.9 * max(plot_df$aumc, na.rm = TRUE),
        #                         y = 0.9 * max(plot_df$dss, na.rm = TRUE))) +
        ggtitle(title_string) +
        theme_paper +
        theme(legend.position = "none")
})

# make combo plot of FIGNL1_1 profile and drug associations
plot_grid(fignl_proteoform_profile_plot, pList[[1]], pList[[3]],
          pList[[4]], pList[[6]], pList[[7]],
          ncol = 3, nrow = 2)

ggsave(filename = here("R/figures/suppl_fig_drugsens_cor_fignl.pdf"), 
       width = 21, height = 14, units = "cm")

# get EPS8L2 specific proteoform dataset
eps8l2_proteoform_df <- proteoform_df %>% 
    filter(grepl("EPS8L2_", gene), !grepl("BR", sample_name_machine)) %>% 
    left_join(nparc_res_hq_df %>% dplyr::select(id, sample_name, conv),
          by = c("gene" = "id", "sample_name_machine" = "sample_name")) %>% 
    filter(conv == TRUE) %>% 
    mutate(sample_name_machine = gsub("_", "-", sample_name_machine))

# plot EPS8L2 proteoforms across cell lines
eps8l2_proteoform_profile_plot <- 
    ggplot(eps8l2_proteoform_df, 
           aes(temperature, rel_value, color = sample_name_machine)) +
    #geom_point() +
    geom_smooth(aes(group = gene, linetype = gene), method = "lm",
                formula = 'y ~ splines::ns(x, df = 4)',
                se = FALSE, alpha = 0.5, size = 0.75) +
    # stat_smooth(aes(group = sample_name_machine),
    #             method = "nls", se = FALSE,
    #             formula = y ~ (1-a)/(1 + exp(-(b/x - c))) + a,
    #             method.args = list(start = c(a = 0, b = 550, c = 10),
    #                                algorithm = 'port'),
    #             geom = "line",
    #             alpha = 0.5) +
    scale_color_manual("Cell line", 
                       values = cl_colors) +
    facet_wrap(~sample_name_machine) +
    theme_paper +
    theme(legend.position = "bottom") +
    labs(x = x_label,
         y = y_label) +
    ggtitle("EPS8L2 proteoforms 1 and 2") 
