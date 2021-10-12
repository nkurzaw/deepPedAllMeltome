library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(Biobase)
library(Rtpca) # install from github: devtools::install_github("nkurzaw/Rtpca")
library(BiocParallel)
library(here)

# source all files in the functions directory
sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = here("R/ppi_coaggregation"))

output_folder <- file.path(here(), "ppi_coaggregation", "output")
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# read in meta file locations
sample_meta_file <- file.path("meta", "sample_meta.txt")
sample_meta_raw <- read_tsv(file = sample_meta_file)

# read in proteoforms
proteoforms <- 
  readRDS(here("proteoform_detection/output/standard/proteoforms_narrow_range_focused.RDS"))

# read in proteoforms with hq nparc fits
alt_hq_nparc_df <- readRDS(here("nparc/output/standard/nparc_res_hq_df.RDS"))

# remove biological replicates from analysis
proteoforms_wo_rep <- 
  proteoforms[,-grep("_BR2", pData(proteoforms)$sample_name_machine)]

# filter proteoforms to only include hq proteoform profiles
proteoform_by_sample_list <- lapply(unique(pData(proteoforms_wo_rep)$sample_name_machine), function(sn){
  hq_by_sample_df <- alt_hq_nparc_df %>% 
    filter(sample_name == sn)
  sample_meta_fil <- sample_meta_raw %>% 
    filter(sample_name_machine == sn)
  p_eset <- proteoforms_wo_rep[which(featureNames(proteoforms_wo_rep) %in% hq_by_sample_df$id),
                    pData(proteoforms_wo_rep)$sample_name_machine == sn]
  return(p_eset)
})
names(proteoform_by_sample_list) <- unique(pData(proteoforms_wo_rep)$sample_name_machine)


# extend string annotation to capture all possible proteoform interactions 
data("string_ppi_df")
string_ppi_975_df <- string_ppi_df %>% 
    filter(combined_score >= 975)

string_ppi_975_proteoform_df <- string_ppi_975_df
  
for(proteo_feature in featureNames(proteoforms)){
  prot_feature <- str_remove(proteo_feature, "_.+")
  isoform_id <- str_remove(proteo_feature, ".+_")
  # find annotations for this protein
  x_matches <- which(string_ppi_975_proteoform_df$x == prot_feature)
  y_matches <- which(string_ppi_975_proteoform_df$y == prot_feature)
  if(length(x_matches) != 0 | length(x_matches) != 0){
    if(isoform_id == "0"){
      if(length(x_matches) != 0){
        string_ppi_975_proteoform_df$x[x_matches] <- proteo_feature
      }
      if(length(y_matches) != 0){
        string_ppi_975_proteoform_df$y[y_matches] <- proteo_feature
      }
      string_ppi_975_proteoform_df$pair[union(x_matches, y_matches)] <- 
      str_replace(string_ppi_975_proteoform_df$pair[union(x_matches, y_matches)],
                  pattern = prot_feature, replacement = proteo_feature)
  }else{
    new_rows <- string_ppi_975_proteoform_df[union(x_matches, y_matches),]
    if(any(new_rows$x == prot_feature)){
      new_rows$x[which(new_rows$x == prot_feature)] <- proteo_feature
    }
    if(any(new_rows$y == prot_feature)){
      new_rows$y[which(new_rows$y == prot_feature)] <- proteo_feature
    }
    new_rows$pair <- str_replace(new_rows$pair, 
                                 pattern = prot_feature, 
                                 replacement = proteo_feature)
    string_ppi_975_proteoform_df <- bind_rows(
      string_ppi_975_proteoform_df, new_rows
    )
  }
  }
}

saveRDS(string_ppi_975_proteoform_df, 
        file.path(output_folder, "string_ppi_975_proteoform_df.RDS"))

# define parallelization
BBPARAM <-  BiocParallel::MulticoreParam(workers = 5, progressbar = TRUE)
# BBPARAM <-  BiocParallel::MulticoreParam(workers = 5, progressbar = TRUE)

# run standard `Rtpca` on all individual cell lines
tpca_result_list <- bplapply(proteoform_by_sample_list, 
                             function(proteoform_by_sample){
  fc_mat <- exprs(proteoform_by_sample) %>% na.omit()
  attributes(fc_mat)$temperature <- pData(proteoforms)$temperature[1:8]

  tpca_obj <- runTPCA(
    objList = list(fc_mat),
    ppiAnno = string_ppi_975_proteoform_df,
    doRocAnalysis = FALSE)

  return(tpcaResultTable(tpca_obj))
})

saveRDS(tpca_result_list, 
        file.path(output_folder, "tpca_result_list_narrow_range_focused_hq_filtered"))

# tpca_result_list <- 
#   readRDS(file.path(output_folder, "tpca_result_list_narrow_range_focused_hq_filtered.RDS"))

# inspect results for 697
print(tpca_result_list[[10]] %>% arrange(p_value) %>% filter(p_adj < 0.1), n = 80)


# run differential `Rtpca` across cell lines: which PPIs show differences across the cell lines?
fc_mat_list <- lapply(proteoform_by_sample_list, 
                             function(proteoform_by_sample){
  fc_mat <- exprs(proteoform_by_sample) %>% na.omit()
  attributes(fc_mat)$temperature <- pData(proteoforms)$temperature[1:8]
  return(fc_mat)
})

# get all interactions which co-melt in at least one of the samples
significantly_coaggregating_ppis_at_least_one_cl <- unique(
  unlist(lapply(tpca_result_list, function(df){
    filter(df, p_adj < 0.1)$complex_name
  }))
)
  
# multi_cell_line_rtpca_df <- 
#   compute_Fstat_PPI_coaggregation_differences(
#     fc_mat_list = fc_mat_list,
#     ppi_anno = filter(string_ppi_975_proteoform_df, 
#                       pair %in% significantly_coaggregating_ppis_at_least_one_cl),
#     BPPARAM = BPPARAM
# )

# run mutli-cell lines comparison
multi_cell_line_rtpca_robust_df <- 
  compute_robust_Fstat_PPI_coaggregation_differences(
    fc_mat_list = fc_mat_list,
    ppi_anno = filter(string_ppi_975_proteoform_df, 
                      pair %in% significantly_coaggregating_ppis_at_least_one_cl),
    BPPARAM = BPPARAM
)

saveRDS(multi_cell_line_rtpca_robust_df, 
        file.path(output_folder, "multi_cell_line_rtpca_robust_df.RDS"))
