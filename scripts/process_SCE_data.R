library(SingleCellExperiment)
library(qs)
library(Matrix)
library(tidyverse)

process_SCE_data <- function(sce_filepath, output_dir) {
  #Load SingleCellExperiment object
  SCE <- qread(sce_filepath)

  #Extract count matrix and metadata
  count <- SCE@assays@data@listData$counts
  writeMM(count, file = file.path(output_dir, "murphycount.mtx"))
  meta <- SCE@colData
  meta.df <- as.data.frame(meta[, 1:26])
  meta.df <- meta.df[, -2]

  #Rename excitatory and inhibitory celltypes as neuronal cell types
  meta.df$cluster_celltype <- sub("^EN-.*", "Neu", meta.df$cluster_celltype)
  meta.df$cluster_celltype <- sub("^IN-.*", "Neu", meta.df$cluster_celltype)
  write.csv(meta.df, file.path(output_dir, "murphy_snRNA_metadata.csv"), quote = FALSE, row.names = FALSE)
  meta.df <- meta.df[, c(1, 2, 25)]
  
  #Extract genes info
  genes_data <- as.data.frame(SCE@rowRanges@elementMetadata@listData)
  genes_data <- genes_data[, c(1, 2, 3)]
  write.csv(genes_data, file.path(output_dir, "murphy_snRNA_genesinfo.csv"), quote = FALSE,row.names = FALSE)

  cell_types <- c("Endo", "Oligo", "Astro", "Micro", "Neu", "OPC")
  col_indices_list <- list(
    which(meta.df$cluster_celltype == "Endo"),
    which(meta.df$cluster_celltype == "Oligo"),
    which(meta.df$cluster_celltype == "Astro"),
    which(meta.df$cluster_celltype == "Micro"),
    which(meta.df$cluster_celltype == "Neu"),
    which(meta.df$cluster_celltype == "OPC")
  )

  #Save metadata and column indices for each celltype
  for (i in seq_along(cell_types)) {
    cell_type <- cell_types[i]
    cols <- col_indices_list[[i]]
    meta_subset <- meta.df[cols, ]
    write.table(meta_subset, file.path(output_dir, paste0(cell_type, "_meta.txt")), quote = FALSE, row.names = FALSE)
    write.table(cols, file.path(output_dir, paste0("col_numbers_", cell_type, ".txt")), quote = FALSE, row.names = FALSE)
  }
}


