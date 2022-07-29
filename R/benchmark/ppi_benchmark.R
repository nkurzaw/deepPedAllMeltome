library(tidyverse)
library(readxl)
library(Rtpca)
library(here)

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