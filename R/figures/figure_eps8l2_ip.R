library(tidyverse)
library(readxl)
library(vsn)
library(limma)
library(here)

# define data path
data_path <- here("data/eps8l2_ip/")

# read in metadata
metadata_df <- read_xlsx(list.files(data_path, full.names = TRUE, pattern = "metadata")) %>% 
    dplyr::select(label, lysate, sample, data_file = `text file`)

# read in MS data
all_df <- bind_rows(lapply(list.files(data_path, full.names = TRUE, pattern = "MSdata"), 
                            function(file_path){
    readin_df <- read_tsv(file_path) %>% 
        dplyr::select(protein_id, gene_name, matches("signal_sum")) %>% 
        filter(!grepl("##", protein_id)) %>% 
        gather(key, value, -protein_id, -gene_name) %>% 
        mutate(label = sub("signal_sum_", "", key),
               data_file = basename(file_path))
}))

# join meta data
full_df <- left_join(all_df, metadata_df, by = c("data_file", "label"))

# create data matrix
full_mat_df <- full_df %>% 
    filter(value > 0, !is.na(value), !is.na(lysate), !is.na(sample)) %>% 
    dplyr::select(-key, -label, -data_file) %>% 
    unite("sample_name", lysate, sample, sep = "_") %>% 
    spread(sample_name, value)

# check if normalization is required
full_mat <- as.matrix(full_mat_df[,-c(1,2)])
boxplot(log2(full_mat))

# vsn norm
vsn_fit <- vsn2(full_mat)
vsn_norm_mat <- predict(vsn_fit, full_mat)
boxplot(vsn_norm_mat)

# create expression set
col_dat <- data.frame(
    cell_line = sub("_.+", "", colnames(vsn_norm_mat)),
    antibody = gsub(".+_", "", colnames(vsn_norm_mat)),
    group = factor(c(rep("high", 12), rep("low", 9))))
rownames(col_dat) <- colnames(vsn_norm_mat)

rownames(vsn_norm_mat) <- paste(full_mat_df$protein_id, full_mat_df$gene_name, sep = "_")

eset <- ExpressionSet(vsn_norm_mat, phenoData = AnnotatedDataFrame(col_dat))

# perform limma analysis
X <- model.matrix(~0 + col_dat$antibody + col_dat$group)
# colnames(X) <- c("ab_abcam", "ab_bethyl", "ab_ig_control", "cl_KASUMI9", 
#                  "cl_KOPN8", "cl_MHH_CALL2", "cl_NALL1", "cl_SEM", "cl_TMD5")
colnames(X) <- c("ab_abcam", "ab_bethyl", "ab_ig_control", "gr_low")
fit <- lmFit(eset, X)

# check abcam 
contrast.matrix.abcam <- makeContrasts(ab_abcam - ab_ig_control, levels = X)
fit_abcam <- contrasts.fit(fit, contrast.matrix.abcam)
fit_abcam <- eBayes(fit_abcam)
volcanoplot(fit_abcam)

abcam_out_df <- topTable(fit_abcam, n=Inf, adjust="BH")
abcam_out_df <- as_tibble(abcam_out_df, rownames = "gene_name") %>% 
    mutate(gene_name_only = sub(".+_", "", gene_name))

ggplot(abcam_out_df, aes(logFC, -log10(P.Value))) +
    geom_point(color = "gray") +
    geom_point(color = "black", data = filter(abcam_out_df, adj.P.Val < 0.10)) +
    ggrepel::geom_text_repel(aes(label = gene_name_only),
                             data = filter(abcam_out_df, adj.P.Val < 0.10 | grepl("ABI1|EPS8", gene_name))) +
    theme_bw()

# check bethyl
contrast.matrix.bethyl <- makeContrasts(ab_bethyl - ab_ig_control, levels = X)
fit_bethyl <- contrasts.fit(fit, contrast.matrix.bethyl)
fit_bethyl <- eBayes(fit_bethyl)
volcanoplot(fit_bethyl)

bethyl_out_df <- topTable(fit_bethyl, n=Inf, adjust="BH")
bethyl_out_df <- as_tibble(bethyl_out_df, rownames = "gene_name") %>% 
    mutate(gene_name_only = sub(".+_", "", gene_name))

ggplot(bethyl_out_df, aes(logFC, -log10(P.Value))) +
    geom_point(color = "gray") +
    geom_point(color = "black", data = filter(bethyl_out_df, adj.P.Val < 0.10)) +
    ggrepel::geom_text_repel(aes(label = gene_name_only),
                             data = filter(bethyl_out_df, adj.P.Val < 0.10 | grepl("ABI1|EPS8", gene_name))) +
    theme_bw()

# turn input data into log fold changes to compare thermal stability groups
abcam_rel_fc_mat_df <- tibble(
    KASUMI2_EPS8L2 = vsn_norm_mat[,1] - vsn_norm_mat[,3],
    KASUMI9_EPS8L2 = vsn_norm_mat[,4] - vsn_norm_mat[,6],
    KOPN8_EPS8L2 = vsn_norm_mat[,7] - vsn_norm_mat[,9],
    MHH_CALL2_EPS8L2 = vsn_norm_mat[,10] - vsn_norm_mat[,12],
    NALL1_EPS8L2 = vsn_norm_mat[,13] - vsn_norm_mat[,15],
    SEM_EPS8L2 = vsn_norm_mat[,16] - vsn_norm_mat[,18],
    TMD5_EPS8L2 = vsn_norm_mat[,19] - vsn_norm_mat[,21],
)

abcam_rel_fc_mat <- as.matrix(abcam_rel_fc_mat_df)
rownames(abcam_rel_fc_mat) <- rownames(vsn_norm_mat)
    
# create expression set
fc_col_dat <- data.frame(
    cell_line = sub("_.+", "", colnames(abcam_rel_fc_mat)),
    group = factor(c(rep("high", 4), rep("low", 3))))
rownames(fc_col_dat) <- colnames(abcam_rel_fc_mat)

fc_eset <- ExpressionSet(abcam_rel_fc_mat, 
                         phenoData = AnnotatedDataFrame(fc_col_dat))

# perform limma analysis
X <- model.matrix(~0 + fc_col_dat$group)
colnames(X) <- c("high", "low")

fit <- lmFit(fc_eset, X)    

# check contrast
contrast.matrix.fc <- makeContrasts(high - low, levels = X)
fit_fc <- contrasts.fit(fit, contrast.matrix.fc)
fit_fc <- eBayes(fit_fc)
volcanoplot(fit_fc)

fc_out_df <- topTable(fit_fc, n=Inf, adjust="BH")
fc_out_df <- as_tibble(fc_out_df, rownames = "gene_name") %>% 
    mutate(gene_name_only = sub(".+_", "", gene_name))

ggplot(fc_out_df, aes(logFC, -log10(P.Value))) +
    geom_point(color = "gray") +
    geom_point(color = "black", data = filter(fc_out_df, adj.P.Val < 0.10)) +
    ggrepel::geom_text_repel(aes(label = gene_name_only),
                             data = filter(fc_out_df, adj.P.Val < 0.10 | grepl("ABI1|EPS8", gene_name))) +
    theme_bw()
