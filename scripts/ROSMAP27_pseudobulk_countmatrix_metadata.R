
#==============================================================
#Pseudobulk 27 count matrices from the 48 sample matrix----
#==============================================================
count_files <- list.files(pattern = "_count_matrix\\.csv$")

for (file in count_files) {
  counts <- read.csv(file, row.names = 1, check.names = FALSE)
  matched_ids <- intersect(colnames(counts), ids27)
  filtered_counts <- counts[, matched_ids, drop = FALSE]
  
  celltype <- sub("_count_matrix\\.csv$", "", file)
  output_filename <- paste0(celltype,length(colnames(filtered_counts)),"_count_matrix.csv")
  
  write.csv(filtered_counts, file = output_filename)}

#=========================================================
#Metadata for pseudobulk each cell type----
#=========================================================
cell_types = c("Neu","Oligo","Astro","Micro","Endo","OPC")

metadata27 = read.csv("metadata_bulk27.csv")
pattern <- paste0("^(", paste(cell_types, collapse = "|"), ")\\d{2}_count_matrix\\.csv$")
count27_files <- list.files(pattern = pattern)

for (file in count27_files) {
  counts <- read.csv(file, row.names = 1, check.names = FALSE)
  metadata_celltype <- metadata27 %>%
    filter(projid %in% colnames(counts))  
  celltype <- sub("\\d{2}_count_matrix$", "", tools::file_path_sans_ext(file))
  output_filename <- paste0(celltype, "_metadata",length(metadata_celltype$projid),".csv")
  
  write.csv(metadata_celltype, file = output_filename, row.names = FALSE)}
