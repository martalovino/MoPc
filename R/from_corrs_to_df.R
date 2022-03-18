#' This function transforms corrs into ESTIMATES and PVALUE DATAFRAMES
#'
#' From list data structure to ESTIMATES data.frame and PVALUE dataframe
#'
#' @param corrs A list with two fields: 1. estimate estimate value
#' 2. p.value pvalue.
#' @param genes_proteins List of gene/proteins
#' @param condv List of miRNAs
#'
#' @return A list with two fields: 1. estimate estimate value
#' 2. p.value pvalue. The two lists have the same length
#'
from_corrs_to_df <- function(corrs, genes_proteins, condv){

  # from list data structure to ESTIMATES data.frame and PVALUE dataframe
  un <- unlist(corrs)
  n <- length(un)/2
  estimates <- un[2*(1:n)-1]
  pvalue <- un[2*(1:n)]

  est_matrix = matrix(estimates, nrow = length(genes_proteins), byrow = TRUE)
  pvalue_matrix = matrix(pvalue, nrow = length(genes_proteins), byrow = TRUE)

  est_df <- data.frame(est_matrix, row.names = genes_proteins)
  colnames(est_df) <- condv
  pval_df <- data.frame(pvalue_matrix, row.names = genes_proteins)
  colnames(pval_df) <- condv
  return(list("estimates" = est_df, "pval" = pval_df))
}
