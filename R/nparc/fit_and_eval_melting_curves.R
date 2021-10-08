#' Fit and evaluate melting curves per celline
#' 
#' @param proteins eset of protein fold changes
#' @param sample_meta data frame of sample meta data
#' @param by_var character string indicating the
#' whether the data should be subsetted by a certain
#' variable. Can be either "cellline", "subtype" or
#' "none"
#' @param subset_var character string of indicating
#' group by which the variable indicated by by_var
#' should be subsetted
#' @param BPPARAM BiocParallel variable defining whether
#' to parallelize NPARC fits
#' 
#' @import NPARC
#' 
#' @return data frame containing model metrics per protein
#' measured within a specific cell line

fit_and_eval_melting_curves <- function(
  proteins, sample_meta, 
  by_var = "cellline",
  subset_var = "697 (EU-3)",
  filter_based_on_alt_model = FALSE,
  alt_model_df = NULL,
  BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)){
  
  stopifnot(by_var %in% c("cellline", "subtype", "none"))
  
  if(by_var == "cellline"){
    # create tidy protein data frame per cell lines
    proteins_df <- biobroom::tidy.ExpressionSet(
      proteins, addPheno = TRUE) %>% 
      filter(sample_name_machine == subset_var) %>% 
      group_by(gene) %>% 
      # re-compute fold changes
      mutate(value = value / value[temperature == "41"]) %>% 
      ungroup %>% 
      # left_join(sample_meta %>% 
      #             dplyr::select(sample_name, subtype),
      #           by = "sample_name") %>% 
      mutate(temperature = as.numeric(temperature)) %>% 
      dplyr::select(uniqueID = gene,
                    relAbundance = value,
                    temperature,
                    channel,
                    set, 
                    sample_name_machine) %>% 
      na.omit() %>% 
      group_by(uniqueID) %>% 
      filter(n() == 8) %>% 
      ungroup
    
  }else if(by_var == "subtype"){
    proteins_df <- biobroom::tidy.ExpressionSet(
      proteins, addPheno = TRUE) %>% 
      group_by(sample_name_machine, gene) %>% 
      # re-compute fold changes
      mutate(value = value / value[temperature == "41"]) %>% 
      ungroup %>% 
      # left_join(sample_meta %>% 
      #             dplyr::select(sample_name, subtype),
      #           by = "sample_name") %>% 
      mutate(temperature = as.numeric(temperature)) %>% 
      dplyr::select(uniqueID = gene,
                    relAbundance = value,
                    temperature,
                    channel,
                    set, 
                    sample_name_machine,
                    subtype) %>% 
      filter(subtype == subset_var) %>% 
      na.omit() %>% 
      group_by(sample_name_machine, uniqueID) %>% 
      filter(n() == 8) %>% 
      ungroup
  }else if(by_var == "none"){
    proteins_df <- biobroom::tidy.ExpressionSet(
      proteins, addPheno = TRUE) %>% 
      group_by(sample_name_machine, gene) %>% 
      # re-compute fold changes
      mutate(value = value / value[temperature == "41"]) %>% 
      ungroup %>% 
      # left_join(sample_meta %>% 
      #             dplyr::select(sample_name, subtype),
      #           by = "sample_name") %>% 
      mutate(temperature = as.numeric(temperature)) %>% 
      dplyr::select(uniqueID = gene,
                    relAbundance = value,
                    temperature,
                    channel,
                    set, 
                    sample_name_machine) %>% 
      na.omit() %>% 
      group_by(sample_name_machine, uniqueID) %>% 
      filter(n() == 8) %>% 
      ungroup
    
    if(filter_based_on_alt_model){
      join_df <- alt_model_df %>% 
        dplyr::select(sample_name_machine = sample_name,
                      uniqueID = id) %>% 
        mutate(hq = 'TRUE')
      
      proteins_df <- proteins_df %>% 
        left_join(join_df, by = c("sample_name_machine", "uniqueID")) %>% 
        filter(hq == 'TRUE')
    }
    
  }
  
  control = NPARC:::getParams()
  
  groups_var <- NULL
  if(by_var == "cellline"){
    groups_var <- proteins_df$sample_name
  }else if(by_var == "subtype"){
    groups_var <- proteins_df$subtype
  }
  
  fit_res <- NPARC:::invokeParallelFits(
    x = proteins_df$temperature, 
    y = proteins_df$relAbundance, 
    id = proteins_df$uniqueID, 
    groups = groups_var,
    BPPARAM = BPPARAM,
    maxAttempts = control$maxAttempts,
    returnModels = FALSE,
    start = control$start)
  
  return(fit_res$modelMetrics)
  
}
