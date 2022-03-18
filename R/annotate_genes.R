#' Order and annotate genes
#'
#' Function to get annotation and order of genes, returning the input data sorted
#' and the chromosome annotation.
#'
#' @param mat1 RNA o prot matrix with unordered genes
#'
#' @return A list with two fields: 1. data the input matrix with only annotated genes
#' ordered by chromosome position 2. annot list with the corresponding
#' chromosome annotation. The number of data rows and the annot list have the same length
#' @export
#'
#' @examples mat1 <- annotate_genes(d$mat1)
annotate_genes <- function(mat1){

  genes_proteins <- rownames(mat1)
  # get ensg of genes of interest
  ensg <- as.character(common_name[genes_proteins, "gene_id"])
  df_mat <- df[df$ensg %in% ensg, ]
  df_mat$common_name <- ensg2common[as.character(df_mat$ensg), "gene_name"]

  gene_order <- as.character(df_mat$common_name[df_mat$common_name %in% genes_proteins])
  gene_annot <- as.character(df_mat$chr[df_mat$common_name %in% genes_proteins])

  mat1n <- mat1[gene_order,]

  return(list("data" = mat1n, "annot" = gene_annot))
}
