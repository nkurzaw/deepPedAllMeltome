library(tidyverse)
library(ggalluvial)

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

# read in sample meta on stages
sample_meta_stages_raw <- read_tsv(file = here("meta/sample_meta_stages.txt")) 

# read in proteoform data
proteoforms <- readRDS(file.path(proteoform_detection_folder, 
                                 "proteoforms_narrow_range_focused.RDS.RDS"))

# get all proteoform data frame
proteoform_df <- biobroom::tidy.ExpressionSet(
    proteoforms, 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    ungroup

# define colors
cl_colors <-  c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", 
                "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", 
                "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", 
                "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455")
# Peter Carl - R-Bloggers 02/2013 "The Paul Tol 21-color salute"
names(cl_colors) <- sort(
    unique(gsub("_", "-", sub("_BR2", "", proteoform_df$sample_name_machine))))

# create full sample meta data frame
sample_meta_df <-  sample_meta_stages_raw %>% 
    left_join(sample_meta_raw %>% 
                  dplyr::select(sample_name, sample_name_machine),
              bt = "sample_name") %>% 
    filter(!grepl("BR2", sample_name)) %>% 
    dplyr::select(sample_name_machine, subtype, 
                  stage = Assigned_stage, sex = gender) %>% 
    mutate(subtype = sub("\\.", "-", subtype)) %>% 
    group_by(subtype) %>% 
    mutate(freq = n()) %>% 
    ungroup()

# remove false positive gene fusion and make stage to factor
sample_meta_df <- sample_meta_df %>% 
    within(subtype[subtype == "PTMA-TMSB4X"] <- "other") %>%
    within(stage[stage == "pre-B other"] <- "pre-B") %>%
    mutate(stage = factor(stage, levels = c("pre-pro-B", "pro-B", "early pre-B", 
                                            "pre-B", "late pre-B"))) %>% 
    mutate(sample_name_machine = gsub("_", "-", sample_name_machine))
    
# make alluvial plot
ggplot(sample_meta_df, aes(y = 1, axis1 = subtype, 
                           axis2 = sample_name_machine,
                           axis3 = stage)) +
    geom_alluvium(aes(fill = sample_name_machine), width = 1/20) +
    geom_stratum(width = 1/20, fill = "gray", color = "gray20") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum)),
               size = 2.75) +
    scale_fill_manual(values = cl_colors) +
    theme_void(base_size = 6) +
    theme(legend.position = "none")

ggsave(filename = here("R/figures/figure_sample_overview.pdf"), 
       width = 18, height = 9, units = "cm")
