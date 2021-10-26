library(here)
library(tidyverse)

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

ppi_coaggregation_folder <- here("ppi_coaggregation/output")

tpca_result_list <- readRDS(
    file.path(ppi_coaggregation_folder,
              "tpca_result_list_narrow_range_focused_hq_filtered.RDS"))

multi_cell_line_rtpca_df <- readRDS(
    file.path(ppi_coaggregation_folder,
              "multi_cell_line_rtpca_robust_df.RDS"))


ppi_volcano <- ggplot(multi_cell_line_rtpca_df, aes(max_rss - min_rss, f_stat)) + 
    geom_point(alpha = 0.25) + 
    geom_point(color = "red", alpha = 0.5,
               data = filter(multi_cell_line_rtpca_df,
                             f_stat > quantile(multi_cell_line_rtpca_df$f_stat, 0.9))) + 
    geom_label_repel(
        label = "CXXC1_2:SETD1A_3",  
        nudge_y = 200,
        direction = "y",
        segment.size = 0.25,
        color = "black", 
        data = filter(multi_cell_line_rtpca_df, pair == "CXXC1_2:SETD1A_3")) +
    scale_x_log10() +
    coord_cartesian(xlim = c(0.005, 2.5)) +
    labs(x =  bquote('RSS'['(n-1)']~' - RSS'['(2)']~''),
         y = expression(''*italic(F)*'-statistic')) +
    theme_paper
