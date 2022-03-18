#' Order and annotate miRNAs
#'
#' Function to get annotation and order of miRNAs, returning the input data sorted
#' and the chromosome annotation.
#'
#' @param mat3 miRNA matrix with unordered miRNAs
#'
#' @return A list with two fields: 1. data the input matrix with only annotated miRNAs
#' ordered by chromosome position 2. annot list with the corresponding
#' chromosome annotation. The number of data rows and the annot list have the same length
#' @export
#'
#' @examples mat3 <- annotate_genes(d$mat3)
annotate_mirnas <- function(mat3){

  mirnas <- rownames(mat3)
  mirna_order <- as.character(df1$mirna[df1$mirna %in% mirnas])
  mirna_annot <- as.character(df1$chr[df1$mirna %in% mirnas])

  mat3n <- mat3[mirna_order,]

  return(list("data" = mat3n, "annot" = mirna_annot))
}
