% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap_corr.R
\name{heatmap_corr}
\alias{heatmap_corr}
\title{Function to plot standard correlation between two matrices in a pdf file}
\usage{
heatmap_corr(
  mat1,
  mat2,
  mat1_name,
  mat2_name,
  column_annot,
  row_annot,
  filename
)
}
\arguments{
\item{mat1}{first matrix}

\item{mat2}{second matrix}

\item{mat1_name}{first matrix name}

\item{mat2_name}{second matrix name}

\item{column_annot}{List containing the chromosome location of each cor_matrix
column}

\item{row_annot}{List containing the chromosome location of each cor_matrix
rows}

\item{filename}{Path and name where the plot will be saved (e.g.,
"/results/Heatmap")}
}
\value{
dataframe with estimates and pvalues
}
\description{
Function to plot standard correlation between two matrices in a pdf file
}
\examples{
d_rna_mirna <- heatmap_corr(mat1, mat3, mat1_name, mat3_name,
condv_annot, gene_annot, "results/mRNA_miRNA")

}
