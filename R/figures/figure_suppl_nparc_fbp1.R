library(here)
library(tidyverse)
library(DESeq2)
library(fgsea)
library(reactome.db)
library(clusterProfiler)

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

# read in sample meta info
sample_meta_stages <- read_tsv(here("meta/sample_meta_stages.txt"))

# load nparc hq results
nparc_res_hq_df <- readRDS(here("nparc/output/standard/nparc_res_hq_df.RDS"))

# load proteoform ratios
proteoforms <- readRDS(
    here("proteoform_detection/output/standard/proteoforms_narrow_range_focused.RDS"))

# define proteoform detection output folder
proteoform_detection_folder <- here("proteoform_detection/output/standard/")

# read in meta data
sample_meta_raw <- read_tsv(file = here("meta/sample_meta.txt"))

# read in peptide data
peptides <- readRDS(file.path(proteoform_detection_folder, "peptides.RDS"))

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

# create full hq auc df
auc_full_hq_df <- nparc_res_hq_df %>% 
    dplyr::select(sample = sample_name, id, aumc) 

# RNA-seq analysis
# read in RNA-seq counts
rna_seq_df <- read_tsv(here("data/ALL_raw_counts_with_samples_info.txt")) %>% 
    column_to_rownames(var = "Gene.Name") %>% 
    dplyr::select(-matches("BR"), -matches("697|COG\\.355|COG\\.394|MHH\\.CALL\\.2")) %>% 
    as.matrix()
colnames(rna_seq_df) <- sub("h$", "", sub("X_", "", gsub("\\.", "_", colnames(rna_seq_df))))

# define coldata
coldata <- data.frame(
    sample_name = colnames(rna_seq_df),
    condition = ifelse(colnames(rna_seq_df) %in% 
                           c("COG_319", "LC4_1", "KASUMI_9", "KASUMI_2", 
                             "MHH_CALL_3", "P30_OHKUBO", "RCH_ACV", "KOPN_8"), 
                       "high", "low")
) %>% 
    left_join(sample_meta_stages %>% 
                  mutate(sample = sub("h$", "", sub("_LL", "", gsub("-", "_", sample_name)))) %>% 
                  dplyr::select(sample_name=sample, sex=gender, stage=Assigned_stage),
              by = "sample_name")

rownames(coldata) <- colnames(rna_seq_df)

# perform DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = rna_seq_df,
                              colData = coldata,
                              design= ~sex + condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "high", "low"))

# get DESeq2 results
plot_df <- res %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    as_tibble() %>% 
    arrange(padj)

# GSEA
# map gene symbols to entrez ids
entrez_map_df <- bitr(geneID = plot_df$rowname, fromType = "SYMBOL", 
                      toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db)

plot_df <- left_join(plot_df, entrez_map_df, by = c("rowname" = "SYMBOL"))
deseq_ranks <- filter(plot_df, !is.na(ENTREZID), !is.na(stat))$stat
names(deseq_ranks) <- filter(plot_df, !is.na(ENTREZID), !is.na(stat))$ENTREZID

# get reactome pathways
pathways <- reactomePathways(names(deseq_ranks))

# perform gene set enrichment
fgseaRes <- fgsea(pathways = pathways, 
                  stats    = deseq_ranks,
                  minSize  = 13,
                  maxSize  = 500)

# check thermal stability of G6PD
g6pd_auc_df <- auc_full_hq_df %>% 
    filter(grepl("G6PD_", id)) %>% 
    filter(!sample %in% c("697","COG_355", "COG_394", "MHH_CALL_2")) %>% 
    mutate(group = ifelse(sample %in%  c("COG_319", "LC4_1", "KASUMI_9", "KASUMI_2", 
                                         "MHH_CALL_3", "P30_OHKUBO", "RCH_ACV", "KOPN_8"), 
                                                         "high", "low")) %>% 
    mutate(group = factor(group, levels = c("low", "high"))) %>% 
    mutate(sample = gsub("_", "-", sample))

ggplot(g6pd_auc_df, aes(group, aumc)) +
    geom_boxplot() +
    geom_jitter(aes(color = sample), width = 0.1) +
    ggsignif::geom_signif(comparisons = list(c("high", "low")), test = t.test) +
    scale_color_manual(values = cl_colors) +
    facet_wrap(~id) +
    labs(x = "FBP1_1 thermal stability",
         y = "Area under the melting curve") +
    theme_paper +
    theme(legend.position = "none")

# check for correlation with FBP1 AUC
fbp1_1_auc_df <- auc_full_hq_df %>% 
    filter(grepl("^FBP1_1", id)) %>% 
    filter(!sample %in% c("697","COG_355", "COG_394", "MHH_CALL_2"))

# check G6PD median peptide melting curves
g6pd_median_peptides_df <- biobroom::tidy.ExpressionSet(
    peptides[which(peptides@featureData$first_protein_id %in% c("ENSP00000377192.3", "ENSP00000377194.2")),], 
    addPheno = TRUE) %>% 
    mutate(temperature = as.numeric(temperature)) %>% 
    group_by(sample_name, gene) %>% 
    mutate(rel_value = value / value[temperature == 41]) %>% 
    group_by(sample_name, temperature) %>% 
    summarize(median_rel_value = median(rel_value, na.rm = TRUE)) %>% 
    ungroup %>% 
    filter(!grepl("_BR", sample_name)) %>% 
    mutate(group = ifelse(sample_name %in%  c("COG-319", "LC4-1", "KASUMI-9", "KASUMI-2", 
                                         "MHH_CALL-3", "P30-OHKUBO", "RCH-ACV", "KOPN-8"), 
                          "high", "low")) 

