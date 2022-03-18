#' Correlation between one gene and all miRNA
#'
#' Function to get the correlation between one gene and all mirna
#'
#' @param gene gene name (COMMON gene notation)
#' @param m corralation method ("pearson", "kendall", "spearman")
#' @param mat1 mRNA matrix
#' @param mat2 prot matrix
#' @param mat2 miRNA matrix
#'
#' @return A list with two fields: 1. estimate estimate value
#' 2. p.value pvalue.
#' @examples a <- corr_func_vect(gene, condv, m, mat1, mat2, mat3)
#'
corr_func_vect<- function(gene, m, mat1, mat2, mat3){
  # print(gene)
  corrs <- sapply(rownames(mat3), cor_fu_single, gene, m, mat1, mat2, mat3)
  return(corrs)
}
