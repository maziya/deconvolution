library(dtangle)

# run Dtangle analysis with sig_matrix and mixture_file being the paths to
# the signature matrix and the expression matrix files respectively
runDtangle <- function(sig_matrix, mixture_file) {

  #Read the signature matrix and mixture file
  X <- read.table(sig_matrix, header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
  Y <- read.table(mixture_file, header = TRUE, sep = ",", row.names = 1, check.names = FALSE)

  #Filter rows where the sum of gene counts is greater than 50
  row_sums <- rowSums(Y)
  Y <- Y[row_sums > 50, ]

  Y <- data.matrix(Y)
  X <- data.matrix(X)
  X <- X[order(rownames(X)), ]
  Y <- Y[order(rownames(Y)), ]
  #Intersect genes between signature and mixture
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  Y <- Y[Ygns %in% Xgns, ]
  X <- X[Xgns %in% row.names(Y), ]

  X_new <- t(X)
  Y_new <- t(Y)

  #Log transform the matrices
  Y_log <- log2(Y_new + 0.5)
  X_log <- log2(X_new + 0.5)

  #Find markers from the expression matrix
  markers <- find_markers(Y = Y_log, references = X_log, pure_samples = NULL,
                          data_type = 'rna-seq', gamma = NULL, marker_method = "diff")

  #Compute estimates
  dt_out <- dtangle(Y_log, references = X_log, pure_samples = NULL, n_markers = 0.01,
                    data_type = 'rna-seq', gamma = NULL, markers = markers,
                    summary_fn = mean)
  base_name <- tools::file_path_sans_ext(basename(mixture_file))
  estimates_filepath <- file.path(getwd(), paste0(base_name, "_Dtangle_estimates.csv"))
  write.csv(dt_out$estimates, estimates_filepath, row.names = TRUE, quote = FALSE)
  return(estimates_filepath)
}
