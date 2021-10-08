export_to_viz <- function (iois,
                           graphs,
                           peptide_features,
                           identifier = "gene_symbols",
                           ensp_column = "protein_ids",
                           membership_colors = c("1" = "#FF5376",
                                                 "2" = "#72AFD9",
                                                 "3" = "#E3D26F",
                                                 "4" = "#A288E3",
                                                 "5" = "#1B5299",
                                                 "6" = "#68D8D6",
                                                 "NA" = "#B78DA3"),
                           output_folder) {
  
  selected_graphs <- graphs[iois]
  
  # --export the nodes
  export_nodes <- function (ioi) {
    g <- graphs[[ioi]]
    
    # get the nodes with attributes
    df_cluster <- data.frame(name = V(g)$name,
                             membership = V(g)$membership,
                             membership_color = membership_colors[as.character(V(g)$membership)],
                             stringsAsFactors = FALSE)
    
    # merge feature meta data
    df_cluster_feature <- df_cluster %>%
      left_join(peptide_features, by = c("name" = "peptide")) %>%
      dplyr::select(peptide = name,
                    !!identifier,
                    protein_id = !!ensp_column,
                    membership,
                    membership_color) %>%
      mutate(protein_id = gsub("\"", "", protein_id)) %>%
      mutate(gene_symbol = ioi)
    
    return(df_cluster_feature)
  }
  
  nodes_list <- sapply(X = names(selected_graphs),
                       FUN = export_nodes,
                       simplify = FALSE)
  
  nodes <- do.call("rbind", nodes_list) %>%
    mutate(id = 0:(nrow(.) - 1)) %>%
    mutate(sequence = gsub("[^A-za-z]", "", peptide))
  
  stopifnot(!duplicated(paste0(nodes$peptide, nodes$gene_symbol)))
  
  
  # --export the edges
  export_edges <- function (ioi) {
    g <- selected_graphs[[ioi]]
    
    ids <- nodes %>%
      filter(gene_symbol == ioi) %>%
      dplyr::select(id, peptide)
    
    edge_df <- g %>%
      as_edgelist() %>%
      as.data.frame() %>%
      mutate_if(is.factor, as.character) %>%
      left_join(ids, by = c("V1" = "peptide")) %>%
      left_join(ids, by = c("V2" = "peptide")) %>%
      mutate(source_id = id.x) %>%
      mutate(target_id = id.y) %>%
      mutate(weight = E(g)$weight) %>%
      group_by(source_id) %>%
      mutate(source_sum_weight = sum(weight, na.rm = TRUE)) %>%
      group_by(target_id) %>%
      mutate(target_sum_weight = sum(weight, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(source_norm_weight = round(weight / source_sum_weight, 3)) %>%
      mutate(target_norm_weight = round(weight / target_sum_weight, 3)) %>%
      mutate(weight = round(weight, 3)) %>%
      dplyr::select(source_id, target_id, weight, source_norm_weight, target_norm_weight)
    
    return(edge_df)
  }
  
  edges_list <- sapply(X = names(selected_graphs),
                       FUN = export_edges,
                       simplify = FALSE)
  
  edges <- do.call("rbind", edges_list)
  
  
  # --get the protein sequences
  # fetch the sequences from ENSEMBL
  # mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  # 
  # sequences <- nodes %>%
  #   group_by(!!!rlang::syms(identifier)) %>%
  #   summarise(protein_id = unique(unlist(str_split(protein_id, ";"))), .groups = "drop") %>%
  #   left_join(getSequence(id = .$protein_id,
  #                         type = "ensembl_peptide_id_version",
  #                         seqType = "peptide",
  #                         mart = mart,
  #                         verbose = TRUE), by = c("protein_id" = "ensembl_peptide_id_version")) %>%
  #   dplyr::rename(sequence = peptide) %>%
  #   mutate(sequence_length = nchar(sequence)) %>%
  #   arrange(!!!rlang::syms(identifier), desc(sequence_length)) %>%
  #   group_by(!!!rlang::syms(identifier)) %>%
  #   mutate(id_rank = 1:n()) %>%
  #   ungroup()
  # 
  # get_sequence_offsets <- function (aligned_sequence) {
  #   aligned_sequence %>%
  #     str_split("") %>%
  #     .[[1]] %>%
  #     set_names(seq_len(length(.))) %>%
  #     .[grepl("[A-Za-z]+", .)] %>%
  #     names() %>%
  #     list()
  # }
  # 
  # # do a multiple sequence alignment
  # aligned_sequences <- sequences %>%
  #   group_by(!!!rlang::syms(identifier)) %>%
  #   summarise(aa_string_set = list(AAStringSet(sequence %>% `names<-`(protein_id))), .groups = "drop") %>%
  #   rowwise() %>%
  #   mutate(msa = list(msa(aa_string_set, method = "ClustalW"))) %>%
  #   ungroup() %>%
  #   .$msa %>%
  #   lapply(FUN = unmasked) %>%
  #   lapply(FUN = as.data.frame) %>%
  #   do.call(rbind, .) %>%
  #   set_names("aligned_sequence") %>%
  #   rownames_to_column("protein_id") %>%
  #   rowwise() %>%
  #   mutate(sequence_offset = get_sequence_offsets(aligned_sequence)) %>%
  #   ungroup() %>%
  #   left_join(sequences, ., by = "protein_id")
  # 
  # # align peptides
  # alignment <- nodes %>%
  #   dplyr::select(id, protein_id, sequence) %>%
  #   mutate(protein_id_split = str_split(protein_id, ";")) %>%
  #   mutate(num_proteins = lapply(protein_id_split, length)) %>%
  #   rowwise() %>%
  #   mutate(df = list(data.frame(id = id, peptide_sequence = sequence, protein_id = unlist(protein_id_split)))) %>%
  #   .$df %>%
  #   do.call(rbind, .) %>%
  #   left_join(aligned_sequences %>% dplyr::select(protein_id, protein_sequence = sequence, sequence_offset), by = "protein_id") %>%
  #   mutate(peptide_start = str_locate_all(protein_sequence, peptide_sequence) %>% lapply(FUN = function (x) paste(x[, 1], collapse = ";")) %>% as.numeric()) %>%
  #   mutate(peptide_length = lapply(peptide_sequence, nchar) %>% as.numeric()) %>%
  #   mutate(peptide_end = peptide_start + peptide_length - 1) %>%
  #   rowwise() %>%
  #   mutate(peptide_aligned = sequence_offset[peptide_start:peptide_end] %>% paste(collapse = ";")) %>%
  #   ungroup() %>%
  #   dplyr::select(id, protein_id, peptide_aligned)
  
  # store
  write_csv(x = nodes,
            file = file.path(output_folder, "nodes.csv"),
            col_names = TRUE,
            quote_escape = FALSE)
  
  write_csv(x = edges,
            file = file.path(output_folder, "edges.csv"),
            col_names = TRUE,
            quote_escape = FALSE)
  
  # write_csv(x = alignment,
  #           path = file.path(output_folder, "alignment.csv"),
  #           col_names = TRUE,
  #           quote_escape = FALSE)
}
