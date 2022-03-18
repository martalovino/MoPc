#' Function to plot standard correlation between two matrices in a pdf file
#'
#' @param mat1 first matrix
#' @param mat2 second matrix
#' @param mat1_name first matrix name
#' @param mat2_name second matrix name
#' @param column_annot List containing the chromosome location of each cor_matrix
#' column
#' @param row_annot List containing the chromosome location of each cor_matrix
#' rows
#' @param filename Path and name where the plot will be saved (e.g.,
#' "/results/Heatmap")
#' @return dataframe with estimates and pvalues
#' @export
#' @examples d_rna_mirna <- heatmap_corr(mat1, mat3, mat1_name, mat3_name,
#' condv_annot, gene_annot, "results/mRNA_miRNA")
#'
heatmap_corr <- function(mat1, mat2, mat1_name, mat2_name, column_annot, row_annot, filename){
  # computing standard corr
  d <- compute_standard_corr(mat1, mat2, mat1_name, mat2_name)
  print(paste("Saving values in:", paste0(filename, "_standard-correlations.RData")))
  readr::write_csv(d$table, paste0(filename, "_standard-correlations.csv"))
  save(d, file = paste0(filename, "_standard-correlations.RData"))

  print("Saving heatmap ..................................")
  plot_heatmap_corr(d$cormatrix, column_annot, row_annot, filename, mat1_name, mat2_name)
  return(d)
}


#' Function to compute standard correlation between two matrices
#'
#' @param mat1 first matrix, genes on rows, samples on columns
#' @param mat2 second matrix, genes on rows, samples on columns. The number of
#' samples must be the same
#' @param mat1_name name of the first matrix, e.g., "RNA"
#' @param mat2_name name of the second matrix, e.g., "miRNA"
#'
#' @return A list with three fields: 1. cormatrix: correlation matrix estimates
#' 2. pmatrix: pvalue matrix, without fdr correction  3. table: dataframe with
#' all previous info in a table
#' @export
#' @examples d <- compute_standard_corr(mat1, mat2, mat1_name, mat2_name)
#'
compute_standard_corr <- function(mat1, mat2, mat1_name, mat2_name){

  matr1 <- t(mat1)
  matr2 <- t(mat2)
  print(paste("Computing standard correlations ...................."))

  mat <- cbind(matr1, matr2)
  t_init <- Sys.time()
  m <- Hmisc::rcorr(mat, type=corr_method)
  t_fin <- Sys.time()
  cat(sprintf("Standard correlation computed in: %s\n", t_fin - t_init))

  cor_matrix <- m[[1]][colnames(matr1), colnames(matr2)]
  p_matrix <- m[[3]][colnames(matr1), colnames(matr2)]

  doc <- reshape2::melt(cor_matrix)
  d <- reshape2::melt(p_matrix)
  doc$pavlue <- d$value
  doc$pvaluefdr <- p.adjust(doc$pavlue, method="fdr")
  colnames(doc) <- c(mat1_name, mat2_name, "estimate", "pvalue", "pvaluefdrcorr")

  return(list("cormatrix" = cor_matrix, "pmatrix" = p_matrix, "table" = doc))
}


#' Function to plot standard correlation between two matrices in a pdf file
#'
#' @param cor_matrix correlation matrix to plot eventually with NAs
#' @param column_annot List containing the chromosome location of each cor_matrix
#' column
#' @param row_annot List containing the chromosome location of each cor_matrix
#' rows
#' @param filename Path and name where the plot will be saved (e.g.,
#' "/results/Heatmap")
#' @param row_name name to be put in the plot for rows, e.g., "mRNA"
#' @param column_name name to be put in the plot for columns, e.g., "miRNA"
#' @export
#'
plot_heatmap_corr <- function(cor_matrix, column_annot, row_annot, filename,
                              row_name, column_name){

  # identifying labels for anno_mark
  convd_pos <- sapply(ch, function(x) ceiling(mean(which(column_annot==x))))
  gene_pos <- sapply(ch, function(x) ceiling(mean(which(row_annot==x))))

  # info for the bottom part of the histogram -----------
  p_occ <- apply(cor_matrix, 2, function(x) sum(x>0, na.rm = TRUE))
  n_occ <- apply(cor_matrix, 2, function(x) sum(x<0, na.rm = TRUE)*(-1))

  ha = ComplexHeatmap::HeatmapAnnotation(positive = ComplexHeatmap::anno_barplot(p_occ, gp = grid::gpar(fill = "#CC0000", col = "#CC0000")),
                         negative = ComplexHeatmap::anno_barplot(n_occ, gp = grid::gpar(fill = "#33CC33", col = "#33CC33" )),
                         height = grid::unit(12, "cm")
  )
  names(ha) = c("Positive\ncorrelations", "Negative\ncorrelations")
  lgd = list(
    at = c(-1, 1),
    labels = c("negative", "positive"),
    title = "Correlation value",
    legend_height = grid::unit(4, "cm"),
    title_position = "lefttop-rot")

  # good picture
  pdf(file = paste0(filename, ".pdf"), height = 18, width = 15)
  htp <- ComplexHeatmap::Heatmap(cor_matrix, na_col = "white",
                 cluster_rows = FALSE, cluster_columns = FALSE,
                 col = circlize::colorRamp2(c(-1,0, 1), c("#003300", "white", "#CC0000")),
                 show_row_names = FALSE,
                 show_column_names = FALSE,
                 name = "Corr",
                 row_title = NULL,
                 column_title = paste0("Correlation ", row_name, "-", column_name),
                 column_title_gp = grid::gpar(fontsize = 30),
                 border = TRUE,
                 right_annotation = ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_mark(at = gene_pos, labels = ch)),
                 top_annotation = ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_mark(at = convd_pos, labels = ch)),
                 bottom_annotation = ha,
                 row_split = factor(row_annot, levels = ch),
                 column_split = factor(column_annot, levels = ch),
                 row_gap = grid::unit(0, "mm"),
                 column_gap = grid::unit(0, "mm"),
                 heatmap_legend_param = lgd
  )
  ComplexHeatmap::draw(htp, row_title = row_name, column_title = column_name, column_title_side = "bottom",row_title_gp = grid::gpar(fontsize = 20), column_title_gp = grid::gpar(fontsize = 20))
  dev.off()
}

