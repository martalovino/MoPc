#' Function to compute standard correlation
#'
#' @param gene  Common gene name
#' @param mat1 First matrix for std corr
#' @param mat2 Second matrix for std corr
#' @param corr_method One of "pearson", "spearman", "kendall"
#'
#' @return A list with 2 fields: estimate and pvalue

normal_cor <- function(gene, mat1, mat2, corr_method){

  p <- stats::cor.test(mat1[gene,], mat2[gene,], method = corr_method)
  return(c(p$estimate, p$p.value))
}

#' Function to compute improvement
#'
#' @param gene  Common gene name
#' @param pvalue  Gene name pvalue
#' @param single_corr_df  single correlation values
#'
#' @return P value improvement

est_minus_single <- function(gene, pval, single_corr_df){
  pval[gene,] - single_corr_df[gene, "sing_cor_pval"]
}


#' Function to fix mirna names
#'
#' put miRNA name in the standard format
#'
#' @param name  miRNA name
#'
#' @return Fixed miRNA name
#' @export
fix_names <- function(name){

  # delete "hsa-", if present
  name <- gsub("hsa-", "", name)
  # check for last 3 characters 5p 3p
  if (stringr::str_sub(name, -3, -1) %in% list("-5p", "-3p")){
    new_name <- name
  } else{
    new_name <- paste0(name, "-5p")
  }
  return(new_name)
}


#' Function to get matrix with only relevant improvements
#'
#' @param diff_th  improvement threashold
#' @param mat1 mRNA matrix
#' @param mat2 prot matrix
#' @param pval pval partial correlation  matrix
#' @param corr_method "pearson", "spearman", "kendall"
#'
#' @return matrix with only relevant improvements (differences between partial correlation pvalue and
#' standard correlation pvalue)
#' @export
#' @examples mat <- get_relevant_improvement(mat1, mat2, pval, corr_method, diff_th)
get_relevant_improvement <- function(mat1, mat2, pval, corr_method, diff_th){
  print("Computing improvement of partial correlation over standard correlation")
  genes_proteins <- rownames(mat1)
  condv <- colnames(pval)
  single_corr <- sapply(genes_proteins, normal_cor, mat1, mat2, corr_method)
  single_corr_df <- data.frame(sing_cor_est=single_corr[1,],
                               sing_cor_pval=single_corr[2,],
                               row.names = genes_proteins)

  #-----------------------------------------------------------------------------
  # Differences between estimates and single correlations P-VALUES
  #-----------------------------------------------------------------------------
  est_min_sin <- sapply(genes_proteins, est_minus_single, pval=pval, single_corr_df = single_corr_df)
  est_min_sin <- t(as.data.frame(est_min_sin))
  rownames(est_min_sin) <- genes_proteins
  colnames(est_min_sin) <- condv

  #-------------------------------------------------------------------------------
  # Selecting only RELEVANT correlations: extreme_tails
  #
  # CONDITIONS: 1) partial correlation improves single correlation PVALUE!!!
  #-------------------------------------------------------------------------------
  # fdr correction on pval
  pval_fdr <- matrix(p.adjust(as.vector(as.matrix(pval)), method = "fdr"), nrow=nrow(pval))
  colnames(pval_fdr) <- colnames(pval)
  rownames(pval_fdr) <- rownames(pval)

  # hist(-log(pval_fdr))
  # hist(est_min_sin)

  # cond1 <- pval_fdr < p_val_th
  # cond2 <- estimates > est_th
  cond3 <- est_min_sin > diff_th

  # sum(cond1 & cond2 * cond3)
  # sum(cond3)

  #-----------------------------------------------------------------------------
  # Get significant corr matrix to plot
  #-----------------------------------------------------------------------------
  mat <- est_min_sin
  colnames(mat) <- colnames(pval_fdr)
  rownames(mat) <- rownames(pval_fdr)
  mat[!cond3]  <- NA

  return(mat)
}
