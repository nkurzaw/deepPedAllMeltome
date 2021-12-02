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
