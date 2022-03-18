#' Correlation between one gene and one miRNA
#'
#' Function to get the correlation between one gene and one mirna
#'
#' @param mi miRNA name
#' @param gene gene name (COMMON gene notation)
#' @param m correlation method ("pearson", "kendall", "spearman")
#' @param mat1 mRNA matrix
#' @param mat2 prot matrix
#' @param mat2 miRNA matrix
#'
#' @return A list with two fields: 1. estimate estimate value
#' @examples a <- cor_fu_single(mi, gene, m, mat1, mat2, mat3)
#' 2. p.value pvalue.
cor_fu_single <- function(mi, gene, m, mat1, mat2, mat3){
  # print(mi)
  a <- ppcor::pcor.test(mat1[gene,], mat2[gene,], mat3[mi,], method = m)
  return (list(a$estimate, a$p.value))
}
