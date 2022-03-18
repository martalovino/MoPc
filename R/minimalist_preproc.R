#' Minimalist preprocessing to MoPc
#'
#' @param mat1 The first matrix. Usually mRNA matrix with samples on the columns and genes on the rows.
#' The gene names must be in the "common name" format (e.g., ABL1)
#' @param mat2 The second matrix. Usually prot matrix with samples on the columns and proteins on the rows.
#' The protein names must be in the "common name" format (e.g., ABL1)
#' @param mat3 The third matrix. Usually miRNA matrix with samples on the columns and miRNAs on the rows.
#' The miRNA names must be in the "hsa-miR-XXXXXX" format (e.g., hsa-miR-200c-3p)
#'
#' @return It quits the execution if there are not enough samples or genes to proceed.
#' Otherwise, it returns a list with the three matrices, properly subset to common genes and common samples.
#' @export
#'
#' @examples d <- minimalist_preproc(rna, prot, mirna)

#'
minimalist_preproc <- function(mat1, mat2, mat3){

  # Common samples
  s <- colnames(mat3)[colnames(mat3) %in% colnames(mat1)[colnames(mat1) %in% colnames(mat2)]]

  if (length(s) > 0){
    cat(sprintf("Total number of samples: %s\n", length(s)))
  } else{
    cat("Not enough samples to proceed. Please check samples' names (column names
        in mat1, mat2, mat3 matrices!\n)")
    quit(status=1)
  }

  # Common genes
  ind = which(rownames(mat2) %in% rownames(mat1))
  genes_proteins = rownames(mat2)[ind]

  if (length(genes_proteins) > 0){
    cat(sprintf("Total number of genes/proteins: %s\n", length(genes_proteins)))
  } else{
    cat("Not enough genes and proteins to proceed. Please check genes' and
    proteins' names (row names in mat1, mat2, mat3 matrices!)\n")
    quit(status=1)
  }

  # Subsetting matrices
  mat1n <- mat1[genes_proteins, s]
  mat2n <- mat2[genes_proteins, s]
  mat3n <- mat3[, s]


  return(list("mat1" = mat1n, "mat2" = mat2n, "mat3" = mat3n))
}
