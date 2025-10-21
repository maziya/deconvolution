#cell_types used
#cell_types = c("Neu","Oligo","Astro","Micro","Endo","OPC")

generate_celltype_pseudobulk_countmatrix <- function(cell_types, output_dir) {
  library(Matrix)
  
  #Load count matrix and gene info
  gene_rows <- read.csv(file.path(output_dir, "murphy_snRNA_genesinfo.csv"), header = TRUE)
  full_count_matrix <- readMM(file.path(output_dir, "murphycount.mtx"))
  gene_counts_list <- list()
  cell_counts_list <- list()
  
  for (cell_type in cell_types) {
    #Read metadata and column indices for this cell type
    col_metadata <- read.table(file.path(output_dir, paste0(cell_type, "_meta.txt")), header = TRUE)
    col_numbers <- read.table(file.path(output_dir, paste0("col_numbers_", cell_type, ".txt")), header = TRUE)
    
    #Subset only required columns
    celltype_count_matrix <- full_count_matrix[, col_numbers$x, drop = FALSE]
    colnames(celltype_count_matrix) <- col_metadata$sample_id
    
    #number of cells per sample
    sample_cell_counts <- as.data.frame(table(colnames(celltype_count_matrix)))
    colnames(sample_cell_counts) <- c("SampleID", "NumCells")
    sample_cell_counts$Celltype <- cell_type
    cell_counts_list[[cell_type]] <- sample_cell_counts
    
    #Create pseudobulk count matrix: transpose -> rowsum -> transpose back
    pseudobulk_matrix <- t(rowsum(t(as.matrix(celltype_count_matrix)), group = colnames(celltype_count_matrix)))
    rownames(pseudobulk_matrix) <- gene_rows$ensembl_gene_id
    
    #Save pseudobulk matrix to file
    output_file <- file.path(output_dir, paste0(cell_type, "_count_matrix.csv"))
    write.csv(pseudobulk_matrix, file = output_file,quote = FALSE)
    
    #Store summary of gene counts and num of cells
    gene_counts <- as.data.frame(t(colSums(pseudobulk_matrix)))
    gene_counts$Celltype <- cell_type
    rownames(gene_counts) <- cell_type
    gene_counts_list[[cell_type]] <- gene_counts
    
    
    rm(celltype_count_matrix, pseudobulk_matrix, col_metadata, col_numbers)
    gc()
  }
  
  rm(full_count_matrix)
  gc()
  
  #Combine count summary and save
  summary_gene_counts_df <- dplyr::bind_rows(gene_counts_list, .id = "cell_type")
  summary_cell_counts_df <- dplyr::bind_rows(cell_counts_list)
  write.table(summary_gene_counts_df,
              file = file.path(output_dir, "celltype_total_genecounts_summary.csv"),
              sep = ",", quote = FALSE, row.names = FALSE)
  write.table(summary_cell_counts_df,
              file = file.path(output_dir, "celltype_num_cells_summary.csv"),
              sep = ",", quote = FALSE, row.names = FALSE)
  
  return(list(gene_counts_list = gene_counts_list,
              cell_counts_list = cell_counts_list))
}
