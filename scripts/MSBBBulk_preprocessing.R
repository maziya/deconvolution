library(dplyr)

#MSBB files from synapse AD Knowledge portal
msbb_filtered1cpm <- read.table("MSBB_Filtered_counts_(greater_than_1cpm).tsv",row.names = 1, header = TRUE)
msbb_rnaseq_metadata <- read.csv("MSBB_assay_rnaSeq_metadata.csv")
msbb_biospecimen <- read.csv("MSBB_biospecimen_metadata.csv")
individual_metadata <- read.csv("MSBB_individual_metadata.csv")

#protein coding genes list
protein_coding = read.csv("protein_coding_ensemble_hgnclist.csv")


#get metadata for a given tissue
get_metadata <- function(tissue_name) {
  msbb_biospecimen %>%
    filter(assay == "rnaSeq", tissue == tissue_name) %>%
    inner_join(msbb_rnaseq_metadata, by = "specimenID") %>%
    filter(exclude == FALSE, RIN >= 5.0) %>%
    select(where(~ !all(is.na(.) | . == ""))) %>%
    inner_join(individual_metadata, by = "individualID") %>%
    mutate(
      pathology = case_when(Braak >= 3 & CERAD < 3 ~ "AD", TRUE ~ "NoAD"),
      ageDeath = ifelse(ageDeath == "90+", "90", ageDeath)
    )
}

#select one sample per individual based on read mapping percentage
get_unique_samples <- function(metadata_df) {
  metadata_df %>%
    mutate(mapped_percentage = (mapped / totalReads) * 100) %>%
    group_by(individualID) %>%
    filter(mapped_percentage == max(mapped_percentage, na.rm = TRUE)) %>%
    ungroup() %>%
    as.data.frame()
}


#Get bulkcounts and metadata and save files
bulkcounts_metadata <- function(tissue_name) {
  metadata_filtered <- get_metadata(tissue_name)
  metadata_unique <- get_unique_samples(metadata_filtered)
  
  #get counts for metadata ids
  metadata_unique <- metadata_unique %>%
    filter(specimenID %in% colnames(msbb_filtered1cpm))
  bulk_counts <- msbb_filtered1cpm[, metadata_unique$specimenID]  
  bulk_counts <- bulk_counts[rownames(bulk_counts) %in% protein_coding$target_id, ]
  
  row_sums <- rowSums(bulk_counts[, 1:(ncol(bulk_counts))])
  
  #Keep only the gene row with max row sum per HGNC_symbol
  bulk_counts$rowsums <- row_sums
  bulk_counts <- as.data.frame(bulk_counts) %>%
    rownames_to_column(var = "target_id") %>%
    inner_join(protein_coding, by = "target_id")
  
  bulk_counts <- bulk_counts %>%
    group_by(HGNC_symbol) %>%
    filter(rowsums == max(rowsums)) %>%
    slice(1) %>%
    ungroup() %>%
    select(-rowsums)%>%
    column_to_rownames(var="target_id")%>%
    select(-HGNC_symbol)
  
  
  tissue_tag <- gsub(" ", "_", tolower(tissue_name))
  write.table(bulk_counts,file = paste0("bulkcounts_", tissue_tag, ".csv"),
              sep = ",", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  write.csv(metadata_unique,file = paste0("metadata_", tissue_tag, ".csv"),
            row.names = FALSE)
}


tissue_list <- c("inferior frontal gyrus", "superior temporal gyrus", 
                 "parahippocampal gyrus", "prefrontal cortex", "frontal pole")
lapply(tissue_list, bulkcounts_metadata)


#Merge bulkcounts and metadata from prefrontal cortex and frontal pole#

bulk_PFC <- read.csv("bulkcounts_prefrontal_cortex.csv", row.names = 1, check.names = FALSE)
meta_PFC <- read.csv("metadata_prefrontal_cortex.csv")

bulk_FP <- read.csv("bulkcounts_frontal_pole.csv", row.names = 1, check.names = FALSE)
meta_FP <- read.csv("metadata_frontal_pole.csv")

bulk_merged <- cbind(bulk_PFC, bulk_FP)
meta_merged <- bind_rows(meta_PFC, meta_FP)

write.csv(bulk_merged, "bulkcounts_PFC_FP.csv", quote = FALSE)
write.csv(meta_merged, "metadata_PFC_FP.csv", row.names = FALSE)



#==========================================================
#samples from common individualIDs across the 4 regions
#==========================================================
#read the metadata files from the 4 regions
metadata_files = list.files(pattern = "metadata_.*\\.csv")
metadata_list = lapply(metadata_files, read.csv)

#Merge all metadata data frames by common individualID
merged_metadata = Reduce(function(x, y) inner_join(x, y, by = "individualID"), metadata_list)
#Select specific columns from merged_metadata
metadata_common = merged_metadata %>%
  select(individualID,matches("^specimenID"),
    RIN.x,sequencingBatch.x,sequencingBatch.x.x,
    sequencingBatch.y,sequencingBatch.y.y,
    pmi.x,libraryPrep.x,
    sex.x,race.x,ageDeath.x,Braak.x,CERAD.x,pathology.x,rRNA.rate.x,ethnicity.x)
colnames(metadata_common) = c('individualID','specimenID.IFG','specimenID.PHG','specimenID.FP','specimenID.STG',
                                "RIN","sequencingBatch.IFG","sequencingBatch.FP","sequencingBatch.PHG","sequencingBatch.STG","pmi","libraryPrep","sex","race","ageDeath",
                                "Braak","CERAD","pathology","rRNA.rate","ethnicity")
write.csv(metadata_common, "metadata_common_individuals.csv", row.names = FALSE)

specimen_ids_list = list(
  PFC_FP = metadata_common$specimenID.FP,
  parahippocampal_gyrus = metadata_common$specimenID.PHG,
  inferior_frontal_gyrus = metadata_common$specimenID.IFG,
  superior_temporal_gyrus = metadata_common$specimenID.STG
)

#read the 4 region bulkcounts file
bulkcounts_files = list.files(pattern = "^bulkcounts_.*\\.csv$")
bulkcounts_list = setNames(lapply(bulkcounts_files, read.csv, row.names = 1, check.names = FALSE),
                            gsub("bulkcounts_(.*)\\.csv", "\\1", bulkcounts_files))  

#Subset common individualID bulkcounts for each region
for (region in names(specimen_ids_list)) {
  ids = specimen_ids_list[[region]]
  if (region %in% names(bulkcounts_list)) {
    counts = bulkcounts_list[[region]]
    ids_in_counts = intersect(ids, colnames(counts))
    filtered_counts = counts[, ids_in_counts, drop = FALSE]
    
    #save common individualID bulkcounts file
    output_file = paste0("bulkcounts_", region, "_common.csv")
    write.csv(filtered_counts, output_file, quote = FALSE)
  }
}
