#' From list to dataframe
#'
#' From list data structure to ESTIMATES data.frame and PVALUE dataframe
#'
#' @param mat1 mRNA matrix
#' @param mat2 prot matrix
#' @param mat3 miRNA matrix
#' @param corr_method Correlation method for partial correlation. The allowed methods are "pearson", "kendall", "spearman"
#' @param out_folder Name of the output folder in which results will be saved
#' @param ncores Number of cores to be used to parallelize the code. Suggested numers are 1 for Windows users, from 1 to max number of cores for Linux users.
#'
#' @return A list with two fields: 1. estimate estimate value
#' 2. p.value pvalue. The two lists have the same length
#' @export
#' @examples a <- compute_pc_values(mat1, mat2, mat3, corr_method, out_folder, ncores)
compute_pc_values <- function(mat1, mat2, mat3, corr_method, out_folder, ncores){

  genes_proteins <- rownames(mat1)
  condv <- rownames(mat3)

  print("Computing partial correlation values...................................")
  # Calculating partial correlations for ALL genes and ALL mirnas
  t_init <- Sys.time()
  corrs <- parallel::mclapply(genes_proteins, corr_func_vect, corr_method, mat1, mat2, mat3, mc.cores = ncores)

  a <- from_corrs_to_df(corrs, genes_proteins, condv)
  t_fin <- Sys.time()
  cat(sprintf("Standard correlation computed in: %s\n", t_fin - t_init))

  # Saving data
  print(paste("Saving partial correlation values in:", file.path(out_folder,
              paste0("partial_estimates_pvalue_", corr_method,".RData"))))

  write.csv(a[[1]], file.path(out_folder, paste0("partial_corr_estimates_", corr_method, ".csv")))
  write.csv(a[[2]], file.path(out_folder, paste0("partial_corr_pvalue_", corr_method, ".csv")))

  save(a, file = file.path(out_folder, paste0("partial_estimates_pvalue_", corr_method, ".RData")))
  print("Data succesfully saved.................................................")
  return(a)
}


