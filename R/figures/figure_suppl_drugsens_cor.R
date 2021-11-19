library(here)
library(tidyverse)

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
    unique(gsub("_BR2", "", proteoform_df$sample_name_machine)))

fignl_proteoform_df <- proteoform_df %>% 
    filter(gene == "FIGNL1_1", !grepl("BR", sample_name_machine))
    

# define axis lavels
x_label <- expression("Temperature"* " " * "("*degree*C*")")
y_label <- "Fraction non-denatured"

# plot FIGNL_1 proteoforms highlighting cell lines with fusion variants
fignl_proteoform_profile_plot <- 
    ggplot(fignl_proteoform_df, 
           aes(temperature, rel_value, color = sample_name_machine)) +
    geom_smooth(aes(group = sample_name_machine), method = "lm",
                formula = 'y ~ splines::ns(x, df = 4)',
                se = FALSE, alpha = 0.5) +
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
    theme(legend.position = "bottom") +
    labs(x = x_label,
         y = y_label) +
    ggtitle("FIGNL proteoform 1") 
