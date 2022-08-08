library(tidyverse)
library(readxl)
library(Rtpca)
library(here)
library(ggpointdensity)

data("string_ppi_df")

# string PPIs
string_ppi_975_df <- string_ppi_df %>% 
    filter(combined_score >= 975)

# read in Becher et al. dataset
supp_tab_becher_s4 <- read_xlsx(
    here("data", "1-s2.0-S0092867418303854-mmc4.xlsx"),
    sheet = "TableS4_TPP-TR")

temperature_anno <- 
    as.numeric(
        gsub("T", "", gsub("_.+", "", colnames(
            supp_tab_becher_s4 %>% 
                dplyr::select(matches("mean\\.fc"))))))

# M
m_df <- supp_tab_becher_s4 %>% 
    filter(cell.cycle == "M") %>% 
    dplyr::select(
        gene_name,
        replicates = found.in.reps,
        max_qupm = max.qupm,
        min_qupm = min.qupm,
        matches("mean\\.fc")) %>% 
    filter(min_qupm > 3,
           replicates == 3)

m_mat <- as.matrix(
    m_df %>% dplyr::select(dplyr::matches("mean\\.fc"))
)
rownames(m_mat) <- m_df$gene_name
attributes(m_mat)$temperature <- temperature_anno


becher_m_ppi_tpca <- runTPCA(
    objList = list(m_mat),
    ppiAnno = string_ppi_975_df,
    nSamp = 10^6)

#saveRDS(becher_m_ppi_tpca, here("R/benchmark/becher_m_ppi_tpca.RDS"))

plotPPiRoc(becher_m_ppi_tpca, computeAUC = TRUE)

# read in Heusel et al. dataset
heusel_df <- read_xlsx(here("data", "1-s2.0-S2405471220300016-mmc2.xlsx"), sheet = 2)

m_heusel_mat <- heusel_df %>% 
    dplyr::select(gene_name = Gene_names, matches("mitosis.intensity.mean")) %>% 
    mutate(gene_name = gsub("\\s.+", "", gene_name)) %>% 
    filter(!is.na(gene_name) & !duplicated(gene_name)) %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "gene_name") %>% 
    as.matrix()

m_heusel_mat_norm <- m_heusel_mat / apply(m_heusel_mat, 1, max)
attributes(m_heusel_mat_norm)$temperature <- 
    seq_len(ncol(m_heusel_mat_norm))

heusel_m_ppi_tpca <- runTPCA(
    objList = list(m_heusel_mat_norm),
    ppiAnno = string_ppi_975_df,
    nSamp = 10^6)

saveRDS(heusel_m_ppi_tpca, here("R/benchmark/heusel_m_ppi_tpca.RDS"))

plotPPiRoc(heusel_m_ppi_tpca, computeAUC = TRUE)

# select 10 random columns to make data size comparable to TPP
set.seed(1)
m_heusel_mat_norm_sub <- m_heusel_mat_norm[,sort(sample(seq(ncol(m_heusel_mat_norm)), 10))]

heusel_sub_m_ppi_tpca <- runTPCA(
    objList = list(m_heusel_mat_norm_sub),
    ppiAnno = string_ppi_975_df,
    nSamp = 10^6)

saveRDS(heusel_sub_m_ppi_tpca, here("R/benchmark/heusel_sub_m_ppi_tpca.RDS"))

plotPPiRoc(heusel_sub_m_ppi_tpca, computeAUC = TRUE)

# intra-complex PPIs
data("ori_et_al_complex_ppis")

becher_m_ori_ppi_tpca <- runTPCA(
    objList = list(m_mat),
    ppiAnno = ori_et_al_complex_ppis,
    nSamp = 10^6)

plotPPiRoc(becher_m_ori_ppi_tpca, computeAUC = TRUE)

saveRDS(becher_m_ori_ppi_tpca, here("R/benchmark/becher_m_ori_ppi_tpca.RDS"))


heusel_sub_m_ori_ppi_tpca <- runTPCA(
    objList = list(m_heusel_mat_norm_sub),
    ppiAnno = ori_et_al_complex_ppis,
    nSamp = 10^6)

saveRDS(heusel_sub_m_ori_ppi_tpca, here("R/benchmark/heusel_sub_m_ori_ppi_tpca.RDS"))

plotPPiRoc(heusel_sub_m_ori_ppi_tpca, computeAUC = TRUE)

# get protein table for one of the deepmeltome cell lines and compare its performance
reh_protein_table <- read_tsv(here("R/benchmark/proteins_table.txt")) %>% 
    dplyr::select(protein_id = `Protein ID`, gene_name = `Gene Name`, 
                  qupm = `Set2_tmt16plex_127N - Amount quanted PSMs`, matches("REH")) %>% 
    filter(qupm > 4) %>% 
    group_by(gene_name) %>% 
    filter(qupm == max(qupm)) %>% 
    ungroup() %>% 
    filter(!duplicated(gene_name)) %>% 
    gather(key, value, matches("X__POOL")) %>% 
    group_by(gene_name) %>% 
    mutate(rel_value = value/value[key == "X__POOL_REH_41C_Set2_tmt16plex_126"]) %>% 
    ungroup() %>% 
    dplyr::select(gene_name, key, rel_value) %>% 
    spread(key, rel_value)

reh_mat <- reh_protein_table %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "gene_name") %>% 
    as.matrix()

reh_ori_ppi_tpca <- runTPCA(
    objList = list(reh_mat),
    ppiAnno = ori_et_al_complex_ppis,
    nSamp = 10^6)

saveRDS(reh_ori_ppi_tpca, here("R/benchmark/reh_ori_ppi_tpca.RDS"))

reh_ori_ppi_tpca_roc_df <- bind_cols(PPiRocTableAnno(reh_ori_ppi_tpca),
                                     PPiRocTable(reh_ori_ppi_tpca))

saveRDS(reh_ori_ppi_tpca_roc_df, here("R/benchmark/reh_ori_ppi_tpca_roc_df.RDS"))


plotPPiRoc(reh_ori_ppi_tpca, computeAUC = TRUE)

inter_heusel_mat <- heusel_df %>% 
    dplyr::select(gene_name = Gene_names, matches("interphase.intensity.mean")) %>% 
    mutate(gene_name = gsub("\\s.+", "", gene_name)) %>% 
    filter(!is.na(gene_name) & !duplicated(gene_name)) %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "gene_name") %>% 
    as.matrix()

inter_heusel_mat_norm <- inter_heusel_mat / apply(inter_heusel_mat, 1, max)
attributes(inter_heusel_mat_norm)$temperature <- 
    seq_len(ncol(inter_heusel_mat_norm))

# select 8 random columns to make data size comparable to deepmeltome TPP
set.seed(8)
inter_heusel_mat_norm_sub <- inter_heusel_mat_norm[,sort(sample(seq(ncol(inter_heusel_mat_norm)), 8))]

heusel_sub_inter_ppi_tpca <- runTPCA(
    objList = list(inter_heusel_mat_norm_sub),
    ppiAnno = ori_et_al_complex_ppis,
    nSamp = 10^6)

plotPPiRoc(heusel_sub_inter_ppi_tpca, computeAUC = TRUE)

heusel_sub_inter_ori_ppi_tpca_roc_df <- bind_cols(
    PPiRocTableAnno(heusel_sub_inter_ppi_tpca),
    PPiRocTable(heusel_sub_inter_ppi_tpca))

saveRDS(heusel_sub_inter_ori_ppi_tpca_roc_df, 
        here("R/benchmark/heusel_sub_inter_ori_ppi_tpca_roc_df.RDS"))

# per PPI scatterplot
heusel_reh_roc_df <- left_join(heusel_sub_inter_ori_ppi_tpca_roc_df, 
                               reh_ori_ppi_tpca_roc_df,
                               by = "pair") %>% 
    na.omit()

ggplot(heusel_reh_roc_df, aes(eucl_dist.x, eucl_dist.y)) +
    geom_point(alpha = 0.2) +
    facet_wrap(~annotated.x) +
    labs(x = "Heusel et al. Euclidean distance",
         y = "REH Deepmeltome Euclidean distance") +
    geom_abline(slope = 1, color = "gray", linetype = "dashed") +
    coord_fixed(xlim = c(0, 15))

ggplot(heusel_reh_roc_df %>% filter(annotated.x), aes(eucl_dist.x, eucl_dist.y)) +
    geom_pointdensity() +
    #facet_wrap(~annotated.x) +
    labs(x = "Heusel et al. Euclidean distance",
         y = "REH Deepmeltome Euclidean distance") +
    geom_abline(slope = 1, color = "gray", linetype = "dashed") +
    coord_fixed(xlim = c(0, 15)) +
    viridis::scale_fill_viridis()


ggplot(heusel_reh_roc_df %>% filter(annotated.x), aes(log10(eucl_dist.x), log10(eucl_dist.y))) +
    geom_point(alpha = 0.2) +
    geom_point(alpha = 0.2, color = "blue", 
               data = heusel_reh_roc_df %>% filter(annotated.x) %>% 
                   filter(grepl("PSM[A,B].+PSM[A,B]", pair))) +
    geom_point(alpha = 0.2, color = "orange",
               data = heusel_reh_roc_df %>% filter(annotated.x) %>%
                   filter(grepl("PSM[C,D].+PSM[C,D]", pair))) +
    facet_wrap(~annotated.x) +
    labs(x = "Heusel et al. Euclidean distance, log10",
         y = "REH Deepmeltome Euclidean distance, log10") +
    geom_abline(slope = 1, color = "gray", linetype = "dashed") +
    coord_fixed(ylim = c(-3, 1), xlim = c(-3, 1)) +
    viridis::scale_fill_viridis()

# membrane protein focused analysis
uniprot_membrane_anno_df <- read_tsv(here("R/benchmark/uniprot_membrane_annotation.tsv")) %>% 
    filter(!is.na(Transmembrane) | !is.na(Intramembrane)) %>% 
    mutate(gene_name = sub(" .+", "", `Gene Names`))

# reh membrane protein filtered
reh_membrane_mat <- reh_protein_table %>% 
    filter(gene_name %in% uniprot_membrane_anno_df$gene_name) %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "gene_name") %>% 
    as.matrix()

reh_membrane_ori_ppi_tpca <- runTPCA(
    objList = list(reh_membrane_mat),
    ppiAnno = ori_et_al_complex_ppis,
    nSamp = 10^6)

plotPPiRoc(reh_membrane_ori_ppi_tpca, computeAUC = TRUE)

inter_membrane_heusel_mat <- heusel_df %>% 
    dplyr::select(gene_name = Gene_names, matches("interphase.intensity.mean")) %>% 
    mutate(gene_name = gsub("\\s.+", "", gene_name)) %>% 
    filter(gene_name %in% uniprot_membrane_anno_df$gene_name) %>% 
    filter(!is.na(gene_name) & !duplicated(gene_name)) %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "gene_name") %>% 
    as.matrix()

inter_membrane_heusel_mat_norm <- inter_membrane_heusel_mat / apply(inter_membrane_heusel_mat, 1, max)
attributes(inter_membrane_heusel_mat_norm)$temperature <- 
    seq_len(ncol(inter_membrane_heusel_mat_norm))

# select 8 random columns to make data size comparable to deepmeltome TPP
set.seed(8)
inter_membrane_heusel_mat_norm_sub <- inter_membrane_heusel_mat_norm[,sort(sample(seq(ncol(inter_membrane_heusel_mat_norm)), 8))]

heusel_membrane_sub_inter_ppi_tpca <- runTPCA(
    objList = list(inter_membrane_heusel_mat_norm_sub),
    ppiAnno = ori_et_al_complex_ppis,
    nSamp = 10^6)

plotPPiRoc(heusel_membrane_sub_inter_ppi_tpca, computeAUC = TRUE)


# correlate distances
heusel_membrane_roc_df <- bind_cols(PPiRocTableAnno(heusel_membrane_sub_inter_ppi_tpca),
                                    PPiRocTable(heusel_membrane_sub_inter_ppi_tpca) )
reh_membrane_roc_df <- bind_cols(PPiRocTableAnno(reh_membrane_ori_ppi_tpca),
                                    PPiRocTable(reh_membrane_ori_ppi_tpca) )

heusel_reh_membrane_roc_df <- left_join(heusel_membrane_roc_df, 
                                        reh_membrane_roc_df,
                                        by = "pair") %>% 
    na.omit()

ggplot(heusel_reh_membrane_roc_df, aes(eucl_dist.x, eucl_dist.y)) +
    geom_point(alpha = 0.2) +
    facet_wrap(~annotated.x) +
    labs(x = "Heusel et al. Euclidean distance",
         y = "REH Deepmeltome Euclidean distance") +
    geom_abline(slope = 1, color = "gray", linetype = "dashed") +
    coord_fixed(xlim = c(0, 15))



