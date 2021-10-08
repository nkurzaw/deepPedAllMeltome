#' Save an \code{ExpressionSet} as plain but annotated txt
#'
#' @param e_set \code{ExpressionSet}
#' @param file destination
#' @param colnames_col identifier for final column names
#' @param rownames_col identifier for final row names
#' @param additional_cols names of additional column annotations
#' @param additional_rows names of additional row annotations
#'
#' @return only file output
#' 
#' @importFrom magrittr "%>%"
#' @importFrom Biobase exprs pData fData
#' 
#' @export
#'
#' @examples
save_e_set_as_txt <- function (e_set,
                               file,
                               colnames_col = "sample_id",
                               rownames_col = "proteoform_id",
                               additional_cols = NULL,
                               additional_rows = NULL) {
  
  # get the storage directory from the file path
  directory <- gsub(pattern = paste0(basename(file), "$"),
                    replacement = "",
                    x = file)
  # ... and create it, if it does not exist
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  
  num_add_rows <- 0
  
  df_save <- e_set %>%
    exprs() %>%
    as.data.frame()
  
  colnames(df_save) <- pData(e_set)[[colnames_col]]
  
  if (!is.null(additional_rows)) {
    # check if columns exist
    if (!all(additional_rows %in% colnames(pData(eSet)))) stop("Column name not existing!")
    
    add_rows <- as.data.frame(pData(e_set)[, additional_rows])
    num_add_rows <- ncol(add_rows)
    add.rows.t <- t(add_rows)
    row.names(add_rows_t) <- additional_rows
    colnames(add_rows_t) <- colnames(df_save)
    df.save <- rbind(add_rows_t, as.matrix(df_save))
  }
  
  if (!is.null(additional_cols)) {
    # check if columns exist
    if (!all(additional_cols %in% colnames(fData(eSet)))) stop("Column name not existing!")
    
    add_cols <- as.data.frame(fData(eSet)[, additional_cols])
    colnames(add_cols) <- additional_cols
    dummy_matrix <- matrix(nrow = num_add_rows, ncol = ncol(add_cols), data = NA)
    colnames(dummy_matrix) <- additional_cols
    if (num.add.rows > 0) row.names(dummy_matrix) <- row.names(df_save[1:num_add_rows, ])
    add_cols_full <- rbind(dummy_matrix, add_cols)
    df.save <- cbind(add_cols_full, df_save)
  }
  
  row.names(df_save)[(num_add_rows + 1):length(row.names(df_save))] <- fData(e_set)[[rownames_col]]
  
  write.table(x = df_save,
              file = file,
              quote = FALSE,
              sep = "\t",
              row.names = TRUE,
              col.names = NA)
}
