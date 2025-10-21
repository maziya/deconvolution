library(dplyr)
library(readxl)

#metadata and count matrix obtained from synapse AD portal
rna_meta <- read.csv("ROSMAP_assay_rnaSeq_metadata.csv")
clinical <- read.csv("ROSMAP_clinical.csv")
biospecimen <- read.csv("ROSMAP_biospecimen_metadata.csv")
counts <- read.table("ROSMAP_all_counts_matrix.txt")

#===============================
#Metadata pre-processing----
#===============================

#Extract specimen IDs from counts matrix (excluding first column)
specimen_ids <- as.data.frame(t(counts[1, -1]))
colnames(specimen_ids) <- "specimenID"
rownames(specimen_ids) <- NULL

#Map specimenID to individualID, keeping only rnaSeq entries
#and those without NA in fields of pmi, RIN
indiv_info <- specimen_ids %>%
  inner_join(biospecimen, by = "specimenID") %>%
  filter(assay == "rnaSeq") %>%
  inner_join(clinical, by = "individualID")%>%
  inner_join(rna_meta, by = "specimenID") %>%
  select_if(~ !all(is.na(.) | . == "")) %>%
  filter(is.na(exclude), !is.na(pmi), !is.na(RIN))

#Classify as AD or NoAD
indiv_pathology <- indiv_info %>%
  mutate(pathology = case_when(
      braaksc >= 3 & ceradsc < 3 & cogdx >= 4 & dcfdx_lv >= 4 ~ "AD",
      TRUE ~ "NoAD"),
    age_death = ifelse(age_death == "90+", "90", age_death))

#Count number in each pathology category
category_counts <- indiv_pathology %>%
  count(pathology)

#Manually change value in library batch column as there are three values 0,6& 7
indiv_pathology[317, "sequencingBatch"] <- 6

#=================================
#Count Matrix pre-processing----
#=================================

#retain only gene IDs and sample IDs from the count matrix
gene_ids <- counts[-(1:5), 1]  
sample_ids <- as.character(counts[1, -1])
count_matrix <- counts[-(1:5), -1]

count_matrix <- apply(count_matrix, 2, as.numeric)
rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- sample_ids

#subset count matrix based on specimenID
IDs <- indiv_pathology$specimenID
count_matrix <- count_matrix[, IDs] 

#remove zero sum rows
row_sums<- rowSums(count_matrix)
count_matrix <- count_matrix[which(row_sums != 0),]

colnames(count_matrix) = indiv_pathology$projid

#extract only protein coding genes from the count matrix
protein_coding = read.csv("protein_coding_ensemble_hgnclist.csv")

count_matrix <- as.data.frame(count_matrix) %>%
  rownames_to_column(var = "target_id") %>%
  mutate(target_id = sub("\\..*", "", target_id)) %>%
  inner_join(protein_coding, by = "target_id")

count_matrix_num <- count_matrix
count_matrix_num[, 2:(ncol(count_matrix_num) - 1)] <- lapply(
  count_matrix_num[, 2:(ncol(count_matrix_num) - 1)],as.numeric)

row_sums <- rowSums(count_matrix_num[, 2:(ncol(count_matrix_num) - 1)])


#Keep only the gene row with max row sum per HGNC_symbol
count_matrix$rowsums <- row_sums
count_matrix <- count_matrix %>%
  group_by(HGNC_symbol) %>%
  filter(rowsums == max(rowsums)) %>%
  slice(1) %>%
  ungroup() %>%
  select(-rowsums)%>%
  column_to_rownames(var="target_id")%>%
  select(-HGNC_symbol)

write.csv(count_matrix, "bulk631_countmatrix.csv",row.names = TRUE,quote = FALSE)

#==================================================================
#Retrieving count matrix for samples with matched snRNA data----
#==================================================================
#mathys et.al supplementary for AD/NoAD status 
rosID_path <- readxl::read_excel("41586_2019_1195_MOESM3_ESM.xlsx", sheet = "S1")

snRNA_biospecimen <- read.csv("snRNAseqPFC_BA10_biospecimen_metadata.csv")
snRNA_idmap <- read.csv("snRNAseqPFC_BA10_id_mapping.csv")
snRNA_idmap_uniq <- snRNA_idmap %>%
  group_by(Subject,projid) %>%
  slice(1) %>%
  select(-fastq)%>%
  ungroup()
snRNA_idmap_uniq <- snRNA_idmap_uniq %>%
  inner_join(rosID_path,by="Subject") %>%
  mutate(`pathologic diagnosis of AD`= case_when(
    `pathologic diagnosis of AD` == "YES" ~ "AD",
    `pathologic diagnosis of AD` == "NO" ~ "NoAD"))

#change the status of pathology in indiv_pathology based on Mathys et.al paper
#for the 27 samples
indiv_pathology$pathology[indiv_pathology$projid %in% snRNA_idmap_uniq$projid] <- 
  snRNA_idmap_uniq$`pathologic diagnosis of AD`[match(indiv_pathology$projid[indiv_pathology$projid %in% snRNA_idmap_uniq$projid], snRNA_idmap_uniq$projid)]


#metadata df for the bulk27 and bulk631 samples
metadata27 = indiv_pathology %>%
  filter(projid %in% snRNA_biospecimen$projid)

write.csv(metadata27, "metadata_bulk27.csv",row.names = FALSE, quote = FALSE)

metadata631 = indiv_pathology
write.csv(metadata631, "metadata_bulk631.csv", row.names = FALSE,quote = FALSE)


#Count matrix for bulk27 samples protein coding
count_matrix27 = count_matrix[,as.character(metadata27$projid)]
row_sums<- rowSums(count_matrix27)
count_matrix27 <- count_matrix27[which(row_sums != 0),]

write.csv(count_matrix27, "bulk27_countmatrix.csv", row.names = TRUE, quote = FALSE)
