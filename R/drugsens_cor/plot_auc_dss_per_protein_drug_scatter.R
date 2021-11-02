plot_auc_dss_per_protein_drug_scatter <- 
  function(auc_mat, dss_mat, sample_anno, protein_id, drug_name){
    
    stopifnot(protein_id %in% rownames(auc_mat) &
                drug_name %in% colnames(dss_mat))
    
    title_string <- paste(drug_name, protein_id, sep = "; ")
    
    auc_fil_df <- auc_mat %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "id") %>% 
      filter(id == protein_id) %>% 
      gather(sample, aumc, -id)
    
    dss_fil_df <- dss_mat %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "sample") %>% 
      gather(drug, dss, -sample) %>% 
      filter(drug == drug_name) 
    
    plot_df <- left_join(
      auc_fil_df,
      dss_fil_df,
      by = "sample") %>% 
      left_join(sample_anno %>% 
                  dplyr::select(sample = sample_name_machine,
                                subtype), 
                by = "sample")
    
    p <- ggplot(plot_df, aes(aumc, dss)) +
      geom_point(aes(color = subtype)) +
      geom_smooth(method = "lm") +
      ggpubr::stat_cor(method = "pearson") +
      ggrepel::geom_text_repel(aes(label = sample)) +
      ggtitle(title_string)
    
    return(p)
  }
