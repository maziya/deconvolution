generate_simulatedbulk_countmatrix <- function(cell_types, input_dir, output_dir) {
  library(dplyr)
  library(tibble)
  library(purrr)
  library(readr)
  
  counts_list <- list()
  sample_gene_list <- list()
  
  for (cell_type in cell_types) {
    #Read per-cell-type pseudobulk matrix (genes x samples)
    pattern <- paste0("^", cell_type, ".*_count_matrix\\.csv$")
    matched_files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
    stopifnot(length(matched_files) == 1)
    file_path <- matched_files[1]
    
    count_matrix <- read.csv(file_path,row.names = 1, check.names = FALSE)
    
    #per-cell-type sample summary
    total_counts <- as.data.frame(t(colSums(count_matrix)))
    total_counts$Celltype <- cell_type
    rownames(total_counts) <- cell_type
    counts_list[[cell_type]] <- total_counts
    
    #Transpose to samples x genes and tag with sample ID
    mat_t <- as.data.frame(t(count_matrix))
    mat_t <- rownames_to_column(mat_t, var = "projid")
    sample_gene_list[[cell_type]] <- mat_t
  }
  
  #Combine all samples
  stacked_df <- bind_rows(sample_gene_list)
  sample_gene_df <- stacked_df %>%
    group_by(projid) %>%
    summarize(across(where(is.numeric), sum), .groups = "drop")
  
  gene_sample_df <- column_to_rownames(sample_gene_df, "projid") %>%
    t() %>%
    as.data.frame()
  
  write.csv(gene_sample_df, file.path(output_dir, "simulatedbulk24_countmatrix.csv"), quote = FALSE)
  return(list(
    gene_sample_df = gene_sample_df
  ))
}
