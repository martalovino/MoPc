
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MoPc

<!-- badges: start -->
<!-- badges: end -->

The goal of MoPc is to …

## Installation

You can install the development version of MoPc from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("martalovino/MoPc")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MoPc)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
# Parameters
corr_method <- "pearson"
mat1_name <- "rna"
mat2_name <- "prot"
mat3_name <- "mirna"
ncores = 7
fdr_th <- 0.05
improvement_th <- 0
out_folder <- "results"
dir.create(out_folder, showWarnings = FALSE)

# load breast cancer dataset
data(rna)
data(prot)
data(mirna)


# minimalist preprocecessing: selecting only common samples and common genes
# between rna and prot
d <- minimalist_preproc(rna, prot, mirna)
#> Total number of samples: 77
#> Total number of genes/proteins: 170

# annotate and sort genes/mirnas by chromosomal order
rna <- annotate_genes(d$mat1)
prot <- annotate_genes(d$mat2)
mirna <- annotate_mirnas(d$mat3)

#-------------------------------------------------------------------------------
# DEFINING ROWS and COLUMNS of partial correlation matrix
# usually genes on the rows and mirnas or CNAs on the columns
#-------------------------------------------------------------------------------
genes_proteins <- rownames(rna$data)
condv <- rownames(mirna$data)
gene_annot <- rna$annot
condv_annot <- mirna$annot

#-------------------------------------------------------------------------------
# PARTIAL CORRELATION COMPUTAION
#-------------------------------------------------------------------------------
a <- compute_pc_values(rna$data, prot$data, mirna$data, corr_method, out_folder, ncores)
#> [1] "Computing partial correlation values..................................."
#> Standard correlation computed in: 2.09219145774841
#> [1] "Saving partial correlation values in: results/partial_estimates_pvalue_pearson.RData"
#> [1] "Data succesfully saved................................................."

mat <- get_relevant_improvement(rna$data, prot$data, a$pval, corr_method, improvement_th)
#> [1] "Computing improvement of partial correlation over standard correlation"
th <- quantile(unlist(mat[!is.na(mat)]), 0.75)
mat <- get_relevant_improvement(rna$data, prot$data, a$pval, corr_method, th)
#> [1] "Computing improvement of partial correlation over standard correlation"


#-------------------------------------------------------------------------------
# GET VALIDATION info
#-------------------------------------------------------------------------------
valid <- validation_info(mat, condv, genes_proteins, 7)
#> [1] "Validating obtained miRNAs"
#> [1] "Validation process ongoing"
#> [1] "mirDB validated: 47"
#> [1] "TargetScan validated: 7"
#> [1] "mirtarbase validated: 29"

#-------------------------------------------------------------------------------
# PARTIAL CORRELATION PLOT
#-------------------------------------------------------------------------------
partial_corr_heatmap_plot(valid$to_validate, genes_proteins, condv, gene_annot, condv_annot,
                          mat, out_folder)
#> [1] "Plot final heatmap"
#> png 
#>   2

#-------------------------------------------------------------------------------
# Computing enrichment significance in mirDB, TargetScan, mirtarBase
#-------------------------------------------------------------------------------
hyper_test(valid$to_validate, valid$universe, valid$mirDB, valid$ts, valid$tarbase, file.path(out_folder, "hypergeometric_test_results.txt"))
#> [1] "Computing Fisher and hypergeometric test .............................."
#> Matrix, Fisher and hypergeometric test on mirDB
#> M1 matrix:
#>     47   1432
#>    204   6793
#> Fisher test p-value: 0.612
#> Hypergeometric test p-value: 0.263
#> 
#> 
#> Matrix, Fisher and hypergeometric test on TargetScan
#> M1 matrix:
#>      7   1472
#>     59   6938
#> Fisher test p-value: 0.191
#> Hypergeometric test p-value: 0.91
#> 
#> 
#> Matrix, Fisher and hypergeometric test on mirtarbase
#> M1 matrix:
#>     29   1450
#>    111   6886
#> Fisher test p-value: 0.312
#> Hypergeometric test p-value: 0.129


#-------------------------------------------------------------------------------
# Get plots and results on mrna-mirna and prot-mirna correlations
#-------------------------------------------------------------------------------
st_corr_mrna_mirna <- heatmap_corr(rna$data, mirna$data, mat1_name, mat3_name,
                                   condv_annot, gene_annot,
                                   file.path(out_folder, "mRNA-miRNA"))
#> [1] "Computing standard correlations ...................."
#> Standard correlation computed in: 1.15331292152405
#> [1] "Saving values in: results/mRNA-miRNA_standard-correlations.RData"
#> [1] "Saving heatmap .................................."

st_corr_prot_mirna <- heatmap_corr(prot$data, mirna$data, mat2_name, mat3_name,
                                   condv_annot, gene_annot,
                                   file.path(out_folder, "prot-miRNA"))
#> [1] "Computing standard correlations ...................."
#> Standard correlation computed in: 0.0117776393890381
#> [1] "Saving values in: results/prot-miRNA_standard-correlations.RData"
#> [1] "Saving heatmap .................................."

#-------------------------------------------------------------------------------
# Validating mrna-mirna and prot-mirna correlations
#-------------------------------------------------------------------------------
to_validate_mrna <- st_corr_mrna_mirna$table[st_corr_mrna_mirna$table$pvaluefdrcorr < fdr_th,]
valid_mrna <- validation_info(to_validate_mrna, condv, genes_proteins, 7, from_cor_mat = FALSE)
#> [1] "Validating obtained miRNAs"
#> [1] "Validation process ongoing"
#> [1] "mirDB validated: 3"
#> [1] "TargetScan validated: 0"
#> [1] "mirtarbase validated: 3"

hyper_test(valid_mrna$to_validate, valid_mrna$universe, valid_mrna$mirDB, valid_mrna$ts, valid_mrna$tarbase, file.path(out_folder, "mRNA_hypergeometric_test_results.txt"))
#> [1] "Computing Fisher and hypergeometric test .............................."
#> Matrix, Fisher and hypergeometric test on mirDB
#> M1 matrix:
#>      3    176
#>    248   8049
#> Fisher test p-value: 0.498
#> Hypergeometric test p-value: 0.782
#> 
#> 
#> Matrix, Fisher and hypergeometric test on TargetScan
#> M1 matrix:
#>      0    179
#>     66   8231
#> Fisher test p-value: 0.406
#> Hypergeometric test p-value: 0.757
#> 
#> 
#> Matrix, Fisher and hypergeometric test on mirtarbase
#> M1 matrix:
#>      3    176
#>    137   8160
#> Fisher test p-value: 0.771
#> Hypergeometric test p-value: 0.343


to_validate_prot <- st_corr_prot_mirna$table[st_corr_prot_mirna$table$pvaluefdrcorr < fdr_th,]
valid_prot <- validation_info(to_validate_prot, condv, genes_proteins, 7, from_cor_mat = FALSE)
#> [1] "Validating obtained miRNAs"
#> [1] "Validation process ongoing"
#> [1] "mirDB validated: 2"
#> [1] "TargetScan validated: 0"
#> [1] "mirtarbase validated: 1"

hyper_test(valid_prot$to_validate, valid_prot$universe, valid_prot$mirDB, valid_prot$ts, valid_prot$tarbase, file.path(out_folder, "prot_hypergeometric_test_results.txt"))
#> [1] "Computing Fisher and hypergeometric test .............................."
#> Matrix, Fisher and hypergeometric test on mirDB
#> M1 matrix:
#>      2     46
#>    249   8179
#> Fisher test p-value: 0.653
#> Hypergeometric test p-value: 0.169
#> 
#> 
#> Matrix, Fisher and hypergeometric test on TargetScan
#> M1 matrix:
#>      0     48
#>     66   8362
#> Fisher test p-value: 1
#> Hypergeometric test p-value: 0.314
#> 
#> 
#> Matrix, Fisher and hypergeometric test on mirtarbase
#> M1 matrix:
#>      1     47
#>    139   8289
#> Fisher test p-value: 0.551
#> Hypergeometric test p-value: 0.188
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
