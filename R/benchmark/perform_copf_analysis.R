library(dplyr)
library(tidyr)
library(readr)
library(data.table)
library(CCprofiler)

df <- read_csv("simulated_pep_cov_15_15_intra_noise_more_intermediate_hard_cases_spread_df.csv") %>%
	mutate(intensity = rel_value * 1000) %>%
	dplyr::select(protein_id = protein_name, peptide_id = peptide,
		      filename = temperature_sample, intensity)

fraction_annotation <- tibble(filename = unique(df$filename))  %>%
	mutate(fraction_number = 1:n())

pepTraces <- importPCPdata(input_data = as.data.table(df),
                           fraction_annotation = as.data.table(fraction_annotation),
                           rm_decoys = FALSE)

traces_corr <- CCprofiler::calculateGeneCorrMatrices(pepTraces)

traces_clustered <- clusterPeptides(
    	traces_corr,
    	method = "average", plot = F, PDF=F,
    	name="ProteoformClusters_deepmeltome_benchmark_")

traces_clusteredInN <- cutClustersInNreal(traces_clustered, clusterN = 2,
                                          min_peptides_per_cluster = 5)

traces_scored <- calculateProteoformScore(traces_clusteredInN)

scores <- unique(traces_scored$trace_annotation[,.(proteoform_score,proteoform_score_z,proteoform_score_pval_adj), by=.(protein_id)])

write_csv(scores, file = "cc_profiler_proteoform_scores.csv")
