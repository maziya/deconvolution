library(tidyr)
library(tibble)

#================================================================
#For ROSMAP 27 samples get actual proportions of cell types----
#================================================================

celltype_samples = read.csv("celltype_num_cells_summary.csv",check.names = FALSE)
celltype_samples <- celltype_samples %>%
  pivot_wider(names_from = SampleID, values_from = NumCells, values_fill = 0)
celltype_samples = column_to_rownames(celltype_samples,var="Celltype")

ids27 = c("20149910","20249897","10536568","11399321","20179164","20275399","21189544","20977678","10222853",
            "11302830","11409232","10101327","20956867","11336574","10260309","21126823","11345331","11342432",
            "11200645","10288185","20104101","20261901","10248033","20207013","20173942","20254740","11609672")

celltypecounts_27 <- celltype_samples[, ids27, drop = FALSE]
csums_27 <- colSums(celltypecounts_27)
frequency_matrix <- celltypecounts_27 / matrix(csums_27, 
                                               nrow = nrow(celltypecounts_27), 
                                               ncol = ncol(celltypecounts_27), 
                                               byrow = TRUE)
frequency_df <- as.data.frame(t(frequency_matrix))
write.csv(frequency_df, "ROSMAPbulk27_actual_celltype_proportions.csv")
