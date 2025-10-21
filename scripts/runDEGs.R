#covariates <- c("RIN", "ageDeath", "pmi", "sex","sequencingBatch") that should be used in the design matrix
#give bulkcounts_file, metadata_file paths and cellproportions_file paths
#sampleID is "projid" for ROSMAP and "specimenID" for MSBB


runDEGs = function(bulkcounts_file, metadata_file, covariate_names,sampleID,cellproportions_file = NULL){
  library(edgeR)
  library(dplyr)
  library(tibble)
  
  #bulkcounts is a df with genes as rows and columns as samples
  bulkcounts = read.csv(bulkcounts_file, header = TRUE, row.names = 1,check.names = FALSE)
  
  #metadata_df has all metadata on the dataset in different columns with one column as sampleID
  metadata_df= read.csv(metadata_file, header = TRUE)
  protein_coding = read.csv("protein_coding_ensemble_hgnclist.csv")

  #filter bulkcounts for protein coding genes and remove zero sum rows
  bulkcounts.num = as.matrix(bulkcounts)
  bulkcounts.num = bulkcounts.num[rownames(bulkcounts.num) %in% protein_coding$target_id,]
  
  row_sums= rowSums(bulkcounts.num)
  bulkcounts.num= bulkcounts.num[which(row_sums != 0),]
  
  #reorder metadata sampleID based on colnames of bulkcounts.num
  samples = colnames(bulkcounts.num)
  metadata_df = metadata_df[match(samples,metadata_df[[sampleID]]), ]
  stopifnot(all(samples == metadata_df[[sampleID]]))
  
  group = factor(metadata_df[["pathology"]])
  group = relevel(group, ref = "NoAD")
  factor_vars = c("msex", "sex", "sequencingBatch","race")
  numeric_vars = c("RIN", "pmi", "ageDeath", "age_death", "rRNA.rate")

  for (cov in covariate_names) {
    if (cov %in% factor_vars) {
      metadata_df[[cov]] = as.factor(metadata_df[[cov]])
    } else if (cov %in% numeric_vars) {
      metadata_df[[cov]] = as.numeric(metadata_df[[cov]])
    } 
  }
  
  #If cell proportion estimates are provided for each sample
  if (!is.null(cellproportions_file)) {
    cellproportions = read.csv(cellproportions_file, header = TRUE, row.names = 1, check.names = FALSE)
    cellproportions = cellproportions[match(samples, rownames(cellproportions)), ]
    stopifnot(all(samples == rownames(cellproportions)))
    cellproportions[] = lapply(cellproportions, as.numeric)

    metadata_df = cbind(metadata_df, cellproportions)
    covariate_formula = paste(c(covariate_names, colnames(cellproportions)), collapse = "+")
  } else {
    covariate_formula = paste(covariate_names, collapse = "+")
  }
  
  design_formula = as.formula(paste("~ group +", covariate_formula))
  design = model.matrix(design_formula, data = metadata_df)
    
  #edgeR analysis
  dge = DGEList(counts = bulkcounts.num,group=group)
  y = calcNormFactors(dge, method = 'TMM')
  y = estimateDisp(y, design,robust=TRUE)
  fit = glmQLFit(y, design = design)
  qlf = glmQLFTest(fit,coef =2)
  pvals = qlf$table
  pvals$adj_pval = stats::p.adjust(pvals$PValue, method = 'BH')
  pvals$target_id = rownames(pvals)
  DEGs = left_join(pvals, protein_coding, by = "target_id")
    
  #save DEG results
  base_name = tools::file_path_sans_ext(basename(bulkcounts_file))
  suffix = if (!is.null(cellproportions_file)) "_CTC" else ""
  output_file = paste0("DEGs_", base_name, suffix, "_", format(Sys.Date(), "%Y%m%d"), ".csv")
  write.csv(DEGs, output_file, row.names = FALSE, quote = FALSE)
}
  
