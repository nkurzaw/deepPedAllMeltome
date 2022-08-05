library(dplyr)
library(tidyr)
library(ggplot2)
library(here)
library(BiocParallel)

simulate_peptide_profiles <- function(protein_name = "test",
                                      proteoform_name = "test_0",
                                      pep_cov = 10, tm = 50, 
                                      within_pf_noise = 0.1, 
                                      noise_sd = 0.1,
                                      slope = 0.75,
                                      temperature_range = 
                                          c(41, 44, 47, 50, 53, 56, 59, 63),
                                      n_cell_lines = 4){
    peptide_tms <- lapply(seq(n_cell_lines), function(i) 
        rnorm(mean = tm, sd = within_pf_noise, n = pep_cov))
    peptide_df <- bind_rows(lapply(seq(n_cell_lines), function(ncell){
        bind_rows(lapply(seq(pep_cov), function(i){
            rounded_tm <- round(peptide_tms[[ncell]][i], 3)
            pep_df <- tibble(protein_name = protein_name,
                   proteoform_name = proteoform_name,
                   peptide = paste(proteoform_name, i, sep = "_"),
                   temperature = seq((41 - rounded_tm), 25, 1) + rounded_tm,
                   rel_value = sapply(1/(1 + exp(slope*seq((41 - rounded_tm), 25, 1))), function(x)
                       rnorm(n = 1, mean = x, sd = noise_sd)))
            return(pep_df)
        })) %>% 
            filter(temperature %in% temperature_range) %>% 
            within(rel_value[temperature == 41] <- 1) %>% 
            # make sure values don't get negative
            within(rel_value[rel_value < 0.1] <- rel_value[rel_value < 0.1] + 0.1) %>% 
            within(rel_value[rel_value == 0] <- 0.01) %>% 
            within(rel_value[rel_value < 0] <-  0.001) %>% 
            mutate(cell_line = paste0("cellline_", ncell))
    }))
       
    return(peptide_df)
}

simulate_protein_with_2_proteoforms <- function(protein_name = "test",
                                                peptide_coverage = 15, 
                                                tm_diff = 1,
                                                melting_point_range = seq(50, 60, 0.1)){
    pf1_pep_cov <- sample(5:(peptide_coverage - 5), 1)
    pf2_pep_cov <- peptide_coverage - pf1_pep_cov
    pf1_tm <- sample(melting_point_range, 1)
    # decide whether pf2_tm is smaller or larger than pf1_tm
    coin_flip <- sample(c(-1,1), 1)
    pf2_tm <- pf1_tm + (coin_flip * tm_diff)
    pf1_peptide_df <- simulate_peptide_profiles(
        protein_name = protein_name,
        proteoform_name = paste(protein_name, "1", sep = "_"),
        pep_cov = pf1_pep_cov, tm = pf1_tm)
    
    pf2_peptide_df <- simulate_peptide_profiles(
        protein_name = protein_name,
        proteoform_name = paste(protein_name, "2", sep = "_"),
        pep_cov = pf2_pep_cov, tm = pf2_tm)
    combo_df <- 
        bind_rows(pf1_peptide_df, pf2_peptide_df)
    
    return(combo_df)
}

simulate_dataset <- function(n_tns = 1000,
                             n_tps = list("1" = 50,
                                          "2" = 50,
                                          "3" = 50,
                                          "4" = 50),
                             peptide_coverage = 15,
                             melting_point_range = seq(50, 60, 0.1),
                             BPPARAM = BiocParallel::MulticoreParam(workers = 4)
                             ){
    # simulate true negatives 
    tn_tms <- sample(melting_point_range, n_tns, replace = TRUE)
    tn_df <- bind_rows(bplapply(seq(n_tns), function(i){
        simulate_peptide_profiles(
            protein_name = paste(c("tn_protein", as.character(i)), collapse = ""),
            proteoform_name = paste(c("tn_protein", as.character(i), "_0"), collapse = ""),
            pep_cov = peptide_coverage,
            tm = tn_tms[i])
    }, BPPARAM = BPPARAM))
    
    # simulate true positives 
    tp_df <- bind_rows(bplapply(names(n_tps), function(tm_diff){
        bind_rows(lapply(seq(n_tps[[tm_diff]]), function(i){
            simulate_protein_with_2_proteoforms(
                protein_name = paste(c("tp_protein_", tm_diff, "_", i), collapse = ""),
                peptide_coverage = peptide_coverage,
                tm_diff = as.numeric(tm_diff)
            )
        }))
    }, BPPARAM = BPPARAM))
    
    tn_tp_df <- bind_rows(tn_df, tp_df)
    return(tn_tp_df)
} 

# convert_simulated_dataset_2_eset <- function(df){
#     
# }

set.seed(123)
test_peptide_df <- simulate_peptide_profiles()
ggplot(test_peptide_df, aes(temperature, rel_value)) + 
    geom_line(aes(color = peptide)) +
    facet_wrap(~cell_line)

# try out simulating a single protein with two proteoforms
test_combo_df <- simulate_protein_with_2_proteoforms(
    peptide_coverage = 15, tm_diff = 1
)

ggplot(test_combo_df, aes(temperature, rel_value,
                          group = interaction(peptide, proteoform_name), 
                          color = proteoform_name)) +
    geom_line() +
    facet_wrap(~cell_line)

test_combo_df <- simulate_protein_with_2_proteoforms(
    peptide_coverage = 15, tm_diff = 2
)

ggplot(test_combo_df, aes(temperature, rel_value,
                          group = interaction(peptide, proteoform_name), 
                          color = proteoform_name)) +
    geom_line() +
    facet_wrap(~cell_line)


# simulate full dataset
set.seed(123)
full_simulated_pep_cov_15_df <- simulate_dataset()

saveRDS(full_simulated_pep_cov_15_df, 
        file = here("R/benchmark/full_simulated_pep_cov_15_df.RDS"))

