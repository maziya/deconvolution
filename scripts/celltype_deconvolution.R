celltype_deconvolution_estimates = function(methods, sig_matrix, mixture_file, celltypes=NULL) {
  #sig_matrix is the MultiBrain signature matrix from Sutton et.al paper
  #sig_matrix = signature_MB.csv file
  library(tidyr)
  library(ggplot2)
  
  deconvolution_list = list()
  for (method in methods) {
    proportions_file = NULL
    if (method == "CIBERSORT") {
      source('runCIBERSORT.R')
      proportions_file = runCIBERSORT(sig_matrix, mixture_file)
    } else if (method == "dtangle") {
      source('runDtangle.R')
      proportions_file = runDtangle(sig_matrix, mixture_file)
    }
  
    proportions = read.csv(proportions_file, header = TRUE, row.names = 1,check.names = FALSE)
    if (is.null(celltypes)) {
      celltypes <- colnames(proportions)
    }
    proportions_selected = proportions[, celltypes,drop = FALSE]
    proportions_selected = as.matrix(proportions_selected)
    base_name <- tools::file_path_sans_ext(basename(proportions_file))
    output_file <- paste0(base_name, "_celltype_estimates.csv")
    write.csv(proportions_selected, file = output_file, row.names = TRUE, quote = FALSE)

    #Plot celltype proportions across samples
    proportions_long <- pivot_longer(as.data.frame(proportions_selected), cols = celltypes,
                                     names_to = "CellType", values_to = "cellproportions")
    plot<- ggplot(proportions_long, aes(x=CellType, y=cellproportions, fill=CellType)) +
      geom_boxplot()+
      labs(title = paste("Cell Type Proportions -", method),
           y = "Proportion",
           x = "Cell Type")
    plot_file <- paste0(base_name, "_celltype_proportions_plot.png")
    ggsave(plot_file, plot = plot, width = 7, height = 5)

    deconvolution_list[[method]] <- list(
      celltype_proportions = proportions_selected,
      proportions_file = proportions_file,
      plot = plot,
      plot_file = plot_file)
  }
  return(deconvolution_list)
}
