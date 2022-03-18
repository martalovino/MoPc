#' Getting validation info
#'
#' Function to know if a gene-mirna couple is validated in mirDB, TargetScan,
#' mirtarBase
#'
#' @param mat  the correlation matrix (eventually with NAs).
#'  Genes on rows - mirnas on columns. Also, it can be a list of
#'  significant genes-mirnas with pvalue (to_validate dataframe).
#' @param condv all mirnas in correlation matrix
#' @param genes_proteins all genes/proteins in correlation matrix
#' @param ncores max number of cores to prallelize the process
#' @param from_cor_mat default to TRUE, set FALSE if you input the to_validate
#' dataframe
#'
#' @return A list with four fields: list(a=to_validate, b=universe, c =mirDB, d=ts, e=tarbase)
#' @export
#'
validation_info <- function(mat, condv, genes_proteins, ncores, from_cor_mat=TRUE){
  print("Validating obtained miRNAs")

  if (from_cor_mat){
    # to_validate <- melt(mat)
    to_validate <-data.table::data.table(
      row = rep(rownames(mat), ncol(mat)),
      col = rep(colnames(mat), each = nrow(mat)),
      value = c(mat))

  }else{
    to_validate <- mat
    colnames(to_validate)[1:2] <- c("row", "col")
  }
  if("value" %in% colnames(to_validate)){
  to_validate <- to_validate[!is.na(to_validate$value),]
  }
  to_validate$validation_name <- paste(to_validate$row,
                                       as.character(parallel::mclapply(to_validate$col, fix_names, mc.cores = ncores)),
                                       sep = "_")
  #-----------------------------------------------------------------------------
  # VALIDATION ON TARGETSCAN, mirDB and MIRTARBASE: report everything in t
  # to_validated. strings are used for the search since the search is faster
  #-----------------------------------------------------------------------------
  # define UNIVERSE of gene_miRNA
  print("Validation process ongoing")
  # universe <- c(NA, length(nrow(mat)*ncol(mat)))
  # for(i in sapply(condv, fix_names)){
  #   universe <- c(universe, paste(genes_proteins, i, sep = "_"))
  # }
  pippo <- unlist(lapply(condv, fix_names))
  universe <- c()
  for(i in pippo){
    universe <- c(universe, paste(genes_proteins, i, sep = "_"))
  }

  # define validated mirDB in UNIVERSE
  mirDB <- mirDB_s[mirDB_s %in% universe]

  # check mirDB validated in TO_VALIDATE
  to_validate$Validated.mirDB <- "no"
  to_validate$Validated.mirDB[to_validate$validation_name %in% mirDB] <- "yes"

  print(paste0("mirDB validated: ", sum(to_validate$Validated.mirDB == "yes")))

  # define validated TargetScan in UNIVERSE
  ts <- ts_s[ts_s %in% universe]

  # check TargetScan validated in TO_VALIDATE
  to_validate$Validated.TS <- "no"
  to_validate$Validated.TS[to_validate$validation_name %in% ts] <- "yes"

  print(paste0("TargetScan validated: ", sum(to_validate$Validated.TS == "yes")))

  # define validated mirTARBASE2020 in UNIVERSE
  tarbase <- tarbase_s[tarbase_s %in% universe]

  # check mirTARBASE2020 validated in TO_VALIDATE
  to_validate$Validated.mirtarbase <- "no"
  to_validate$Validated.mirtarbase[to_validate$validation_name %in% tarbase] <- "yes"

  print(paste0("mirtarbase validated: ", sum(to_validate$Validated.mirtarbase == "yes")))

  # # save to list
  # write_csv(to_validate, paste0("results/", tissue, "_validated_list.csv"))
  # print(paste0("Validation info written in: ", paste0("results/", tissue, "_validated_list.csv")))

  return(list("to_validate" = to_validate, "universe"=universe, "mirDB" =mirDB, "ts"=ts, "tarbase"=tarbase))
}
