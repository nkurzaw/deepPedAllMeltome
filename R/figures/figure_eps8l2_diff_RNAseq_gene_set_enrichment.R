library(DESeq2)
library(tidyverse)
library(readxl)
library(vsn)
library(limma)
library(here)
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

# define output folder for figures
figure_output_folder <- here("R/figures")

# read in sample meta info
sample_meta_stages <- read_tsv(here("meta/sample_meta_stages.txt"))

# read in RNA-seq counts
rna_seq_df <- read_tsv(here("data/ALL_raw_counts_with_samples_info.txt")) %>% 
    column_to_rownames(var = "Gene.Name") %>% 
    dplyr::select(-matches("BR")) %>% 
    as.matrix()
colnames(rna_seq_df) <- sub("h$", "", sub("X_", "", gsub("\\.", "_", colnames(rna_seq_df))))

# subset dataset to relevant cell lines
relevant_cell_lines <- c("SEM", "TMD5", "NALL_1", "KOPN_8", "HAL_01", "LC4_1",
                         "KASUMI_2", "MHH_CALL_2", "KASUMI_9", "MHH_CALL_3")
rna_seq_df <- rna_seq_df[,colnames(rna_seq_df) %in% relevant_cell_lines]

# define coldata
coldata <- data.frame(
    sample_name = colnames(rna_seq_df),
    condition = case_when(
        colnames(rna_seq_df) %in% c("SEM", "TMD5", "NALL_1", "HAL_01") ~ "low", 
        colnames(rna_seq_df) %in% c("KOPN_8", "KASUMI_2", "MHH_CALL_3",
                                    "MHH_CALL_2", "KASUMI_9", "LC4_1") ~"high",
        TRUE ~ "none")
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
entrez_map_df <- clusterProfiler::bitr(geneID = plot_df$rowname, fromType = "SYMBOL", 
                                       toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db)

plot_df <- left_join(plot_df, entrez_map_df, by = c("rowname" = "SYMBOL"))
deseq_ranks <- filter(plot_df, !is.na(ENTREZID), !is.na(stat))$stat
names(deseq_ranks) <- filter(plot_df, !is.na(ENTREZID), !is.na(stat))$ENTREZID

# get reactome pathways
pathways <- reactomePathways(names(deseq_ranks))

# perform gene set enrichment
fgseaRes <- fgsea(pathways = pathways, 
                  stats    = deseq_ranks,
                  minSize  = 15,
                  maxSize  = 500)

# plot NF-kappaB enrichment plot
plotEnrichment(pathways[["Activation of NF-kappaB in B cells"]],
               deseq_ranks) + labs(title="Activation of NF-kappaB in B cells") +
    theme_paper

# plot Immunoregulatory interaction enrichment plot
plotEnrichment(pathways[[944]],
               -deseq_ranks) + 
    labs(title="Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell") +
    theme_paper

ggsave(file.path(figure_output_folder, "fig6_immunregulatory_interaction_enrichment.pdf"), 
       width = 8, height = 5, units = "cm")


# plot Interleukin-10 enrichment plot
plotEnrichment(pathways[["Signaling by Interleukins"]],
               -deseq_ranks) + 
    labs(title="Signaling by Interleukins") +
    theme_paper

ggsave(file.path(figure_output_folder, "fig6_signaling_interleukins_enrichment.pdf"), 
       width = 8, height = 5, units = "cm")

# plot Interleukin-10 enrichment plot
plotEnrichment(pathways[["Interleukin-10 signaling"]],
               -deseq_ranks) + 
    labs(title="Interleukin-10 signaling") +
    theme_paper

ggsave(file.path(figure_output_folder, "fig6_signaling_interleukin10_enrichment.pdf"), 
       width = 8, height = 5, units = "cm")

# make Volcano plot with colored dots of Interleukin-10 signaling
plot_df <- plot_df %>% 
    mutate(interleukin10 = ENTREZID %in% pathways[["Interleukin-10 signaling"]])

ggplot(plot_df, aes(log2FoldChange, -log10(pvalue))) +
    geom_point(color = "gray", alpha = 0.2) +
    geom_point(color = "green", data = filter(plot_df, interleukin10)) +
    theme_paper

# plot Interleukin-10 enrichment plot
plotEnrichment(pathways[["Signaling by Interleukins"]],
               deseq_ranks) + 
    labs(title="Signaling by Interleukins") +
    theme_paper

ggsave(file.path(figure_output_folder, "fig6_signaling_interleukins_depletion.pdf"), 
       width = 8, height = 5, units = "cm")

# plot Interleukin-10 enrichment plot
plotEnrichment(pathways[["Interleukin-10 signaling"]],
               deseq_ranks) + 
    labs(title="Interleukin-10 signaling") +
    theme_paper

ggsave(file.path(figure_output_folder, "fig6_signaling_interleukin10_depletion.pdf"), 
       width = 8, height = 5, units = "cm")


high_eps8l2_genes <- filter(plot_df, padj < 0.1, log2FoldChange > 0, !is.na(ENTREZID))$ENTREZID
universe_genes <- filter(plot_df, !is.na(ENTREZID))$ENTREZID

# GO analysis
ego <- enrichGO(gene          = high_eps8l2_genes,
                universe      = universe_genes,
                OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1,
                readable      = TRUE)

dotplot(ego)
