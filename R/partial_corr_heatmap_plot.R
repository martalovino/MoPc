#' Function to plot standard correlation between two matrices in a pdf file
#'
#' @param to_validate dataframe to plot eventually with NAs
#' @param genes_proteins list of all genes
#' @param condv list of all mirnas
#' @param gene_annot List containing the chromosome location of each cor_matrix
#' column
#' @param condv_annot List containing the chromosome location of each cor_matrix
#' rows
#' @param mat matrix
#' @param out_folder Path to the desired output folder
#' @export
partial_corr_heatmap_plot <- function(to_validate, genes_proteins, condv,
                                      gene_annot, condv_annot, mat, out_folder){
  print("Plot final heatmap")
  v <- to_validate[(to_validate$Validated.mirDB == "yes" |
                      to_validate$Validated.TS == "yes" |
                      to_validate$Validated.mirtarbase == "yes"), ]

  v1 <- v[, 1:2]
  v1$value <- 1

  mat_v <- as.data.frame(data.table::dcast(v1, row~col, mean))
  rownames(mat_v) <- mat_v$row
  mat_v <- mat_v[,2:ncol(mat_v)]

  # remove duplicated in gene and condv, in mat too
  ind <- !duplicated(genes_proteins)
  genes_proteins_1 <- genes_proteins[ind]
  gene_annot_1 <- gene_annot[ind]

  ind <- !duplicated(condv)
  condv_1 <- condv[ind]
  condv_annot_1 <- condv_annot[ind]

  matt <- as.data.frame(mat[genes_proteins_1, condv_1])
  p <- as.numeric(unlist(matt))
  p1 <- matrix(p, nrow = nrow(matt))
  rownames(p1) <- rownames(matt)
  colnames(p1) <- colnames(matt)
  matt <- p1

  # rescale mat_v to match with mat
  missing_genes <- genes_proteins_1[!genes_proteins_1 %in% rownames(mat_v)]
  missing_condv <- condv_1[!condv_1 %in% colnames(mat_v)]

  #adding missing rows to mat_v
  mr <- matrix(NaN, nrow = length(missing_genes), ncol = ncol(mat_v))
  rownames(mr) <- missing_genes
  colnames(mr) <- colnames(mat_v)
  mat_v <- rbind(mat_v, mr)

  # adding missing columns to mat_v
  mc <- matrix(NaN, nrow = nrow(mat_v), ncol = length(missing_condv))
  rownames(mc) <- rownames(mat_v)
  colnames(mc) <- missing_condv
  mat_v <- cbind(mat_v, mc)

  # order properly
  mat_v_o <- mat_v[genes_proteins_1, condv_1]
  matt[mat_v_o == 1] <- -1

  # identifying labels for anno_mark
  convd_pos <- sapply(ch, function(x) ceiling(mean(which(condv_annot_1==x))))
  gene_pos <- sapply(ch, function(x) ceiling(mean(which(gene_annot_1==x))))

  p_occ <- apply(matt, 2, function(x) sum(x>0, na.rm = TRUE))

  ha = ComplexHeatmap::HeatmapAnnotation(positive = ComplexHeatmap::anno_barplot(p_occ, gp = grid::gpar(fill = "gray", col = "gray")),
                                         height = grid::unit(6, "cm")
  )
  names(ha) = c("Positive\ncorrelations")

  lgd = list(
    at = c(-1, 1),
    labels = c("Validated", "Statistically\nsignificant"),
    title = "Value",
    legend_height = grid::unit(2, "cm"),
    title_position = "lefttop-rot")

  # good picture
  pdf(file = file.path(out_folder, "partial_corr_figure_validated.pdf"), height = 18, width = 15)
  htp <- ComplexHeatmap::Heatmap(matt, na_col = "white",
                                 cluster_rows = FALSE, cluster_columns = FALSE,
                                 col = circlize::colorRamp2(c(-1,0, 1), c("#003300", "white", "#CC0000")),
                                 show_row_names = FALSE,
                                 show_column_names = FALSE,
                                 name = "Corr",
                                 row_title = NULL,
                                 column_title = "Partial correlation",
                                 column_title_gp = grid::gpar(fontsize = 30),
                                 border = TRUE,
                                 right_annotation = ComplexHeatmap::rowAnnotation(foo = ComplexHeatmap::anno_mark(at = gene_pos, labels = ch)),
                                 top_annotation = ComplexHeatmap::HeatmapAnnotation(foo = ComplexHeatmap::anno_mark(at = convd_pos, labels = ch)),
                                 bottom_annotation = ha,
                                 row_split = factor(gene_annot_1, levels = ch),
                                 column_split = factor(condv_annot_1, levels = ch),
                                 row_gap = grid::unit(0, "mm"),
                                 column_gap = grid::unit(0, "mm"),
                                 heatmap_legend_param = lgd
  )
  ComplexHeatmap::draw(htp, row_title = paste0(mat1_name, "-", mat2_name), column_title = mat3_name, column_title_side = "bottom",row_title_gp = grid::gpar(fontsize = 20), column_title_gp = grid::gpar(fontsize = 20))
  dev.off()
}
