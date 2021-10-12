filter_drugsens_for_minimal_effect <- function(dss_df, min_effect = 6){
    dss_fil_df <- dss_df %>% 
        group_by(drug_name) %>% 
        filter(any(dss > min_effect)) %>% 
        ungroup
    
    return(dss_fil_df)
}
