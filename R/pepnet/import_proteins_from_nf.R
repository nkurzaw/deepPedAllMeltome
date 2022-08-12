import_proteins_from_nf <- function (file,
                                     sample_meta_file,
                                     id_col = "Gene Name",
                                     protein_id_col = "Protein ID(s)",
                                     quan_regex = "X__POOL_.+tmt16plex",
                                     sample_id_col = "sample_id") {
    
    # read the protein file
    data_raw <- fread(file = file,
                      sep = "\t",
                      quote = "", 
                      header = TRUE,
                      stringsAsFactors = FALSE,
                      data.table = FALSE)
    
    # filter for relevant columns and some renaming
    data <- data_raw %>%
        dplyr::select(id = !!id_col,
                      protein_ids = !!protein_id_col,
                      matches("q-value"),
                      matches("precursor"),
                      matches(quan_regex)) %>%
        mutate(protein_ids = gsub("\\(pre=[A-Z\\-],post=[A-Z\\-]\\)", "", protein_ids)) %>%
        mutate(id = gsub(";;+", ";", id)) %>%
        mutate(id = gsub("^;", "", id)) %>%
        mutate(id = gsub(";$", "", id))  %>%
        separate(id, into = c("first_id"), sep = ";", remove = FALSE, extra = "drop") %>%
        separate(protein_ids,
                 into = "first_protein_id",
                 sep = ";",
                 extra = "drop",
                 remove = FALSE)
    
    # split off quantitative data
    quan_data <- data %>%
        dplyr::select(matches(quan_regex)) %>%
        as.matrix()
    colnames(quan_data) <- sub("X__POOL_.+Set", "Set",  colnames(quan_data))
        
    
    # split off feature data
    feature_data <- data %>%
        dplyr::select(-matches(quan_regex))
    
    # load sample meta file
    pheno_data <- read_tsv(file = sample_meta_file) %>%
        dplyr::select(sample_id = !!sample_id_col, everything()) %>%
        mutate(rownames = sample_id) %>%
        column_to_rownames("rownames")
    
    # order quan data to match pheno data
    quan_data <- quan_data %>% 
        .[, row.names(pheno_data)]
    
    # build ExpressionSet
    e_set <- ExpressionSet(assayData = quan_data,
                           phenoData = AnnotatedDataFrame(pheno_data),
                           featureData = AnnotatedDataFrame(feature_data))
    
    return(e_set)
}