#' Computing hypergeometric test
#'
#' Function to compute te hyper geometric test
#'
#' @param to_validate  dataframe of elements to be validated.
#' @param universe all mirnas_gene couples potentially validated.
#' @param mirDB list of mirDB validated elements
#' @param ts list of TargetScan validated elements
#' @param tarbase list of mirTarBase validated elements
#' @param outfile filename with the path of the folder which contains the output
#' @examples  hyper_test(to_validate, universe, mirDB, ts, tarbase, out_folder)
#'
#' @export
#'
hyper_test <- function(to_validate, universe, mirDB, ts, tarbase, outfile){
  #-------------------------------------------------------------------------------
  # COMPUTING HYPERGEOMETRIC and FISHER TEST for mirDB
  # strings are used for the search since it is faster
  # M1
  # x  x3
  # x2  x4
  #-------------------------------------------------------------------------------
  to_save <- c()
  print("Computing Fisher and hypergeometric test ..............................")
  x <- sum(to_validate$Validated.mirDB == "yes")
  x2 <- length(mirDB) - x
  if (x2 < 0){
    x2 <- 0
  }
  x3 <- nrow(to_validate) - x
  x4 <- length(universe) - nrow(to_validate) - x2

  # Fisher exact test
  M1<-matrix(c(x, x2, x3, x4),nrow =2)
  # writeLines("\n\nMatrix, Fisher and hypergeometric test on mirDB")
  to_save <- c(to_save, "Matrix, Fisher and hypergeometric test on mirDB\nM1 matrix:")
  to_save <- c(to_save, sprintf("%6s %6s\n%6s %6s", M1[1,1], M1[1,2], M1[2,1], M1[2,2]))

  # print(M1)
  z <- fisher.test(M1)
  to_save <- c(to_save, sprintf("Fisher test p-value: %s", round(z$p.value, 3)))

  # hypergeometric test
  k <- nrow(to_validate)
  m <- length(mirDB)
  n <- length(universe) - length(mirDB)
  # print("Dhyper pvalue:")
  # print(dhyper(x, m, n, k))
  # print("Phyper pvalue:")
  z <- phyper(x, m, n, k, lower.tail = FALSE)
  to_save <- c(to_save, sprintf("Hypergeometric test p-value: %s\n\n", round(z, 3)))

  #-------------------------------------------------------------------------------
  # COMPUTING HYPERGEOMETRIC and FISHER TEST for TargetScan
  # strings are used for the search since it is faster
  # M1
  # x1  x3
  # x2  x4
  #-------------------------------------------------------------------------------
  x <- sum(to_validate$Validated.TS == "yes")
  x2 <- length(ts) - x
  if (x2 < 0){
    x2 <- 0
  }
  x3 <- nrow(to_validate)- x
  x4 <- length(universe) - nrow(to_validate) - x2

  # Fisher exact test
  M1<-matrix(c(x, x2, x3, x4),nrow =2)
  # writeLines("\n\nMatrix, Fisher and hypergeometric test on TargetScan")
  to_save <- c(to_save, "Matrix, Fisher and hypergeometric test on TargetScan\nM1 matrix:")
  to_save <- c(to_save, sprintf("%6s %6s\n%6s %6s", M1[1,1], M1[1,2], M1[2,1], M1[2,2]))
  # print(M1)
  z <- fisher.test(M1)
  to_save <- c(to_save, sprintf("Fisher test p-value: %s", round(z$p.value, 3)))

  # hypergeometric test
  k <- nrow(to_validate)
  m <- length(ts)
  n <- length(universe) - length(ts)
  # print("Dhyper pvalue:")
  # print(dhyper(x, m, n, k))
  # print("Phyper pvalue:")
  # print(phyper(x, m, n, k, lower.tail = FALSE))
  z <- phyper(x, m, n, k, lower.tail = FALSE)
  to_save <- c(to_save, sprintf("Hypergeometric test p-value: %s\n\n", round(z, 3)))

  #-------------------------------------------------------------------------------
  # COMPUTING HYPERGEOMETRIC and FISHER TEST for mirtarbase
  # strings are used for the search since it is faster
  #-------------------------------------------------------------------------------
  x <- sum(to_validate$Validated.mirtarbase == "yes")
  x2 <- length(tarbase) - x
  if (x2 < 0){
    x2 <- 0
  }
  x3 <- nrow(to_validate) - x
  x4 <- length(universe) - nrow(to_validate) - x2

  # Fisher exact test
  M1<-matrix(c(x, x2, x3, x4),nrow =2)
  # writeLines("\n\nMatrix, Fisher and hypergeometric test on mirtarbase")
  to_save <- c(to_save, "Matrix, Fisher and hypergeometric test on mirtarbase\nM1 matrix:")
  to_save <- c(to_save, sprintf("%6s %6s\n%6s %6s", M1[1,1], M1[1,2], M1[2,1], M1[2,2]))
  # print(M1)
  z <- fisher.test(M1)
  to_save <- c(to_save, sprintf("Fisher test p-value: %s", round(z$p.value, 3)))

  # hypergeometric test
  k <- nrow(to_validate)
  m <- length(tarbase)
  n <- length(universe) - length(tarbase)
  # print("Dhyper pvalue:")
  # print(dhyper(x, m, n, k))
  # print("Phyper pvalue:")
  # print(phyper(x, m, n, k, lower.tail = FALSE))
  z <- phyper(x, m, n, k, lower.tail = FALSE)
  to_save <- c(to_save, sprintf("Hypergeometric test p-value: %s", round(z, 3)))

  writeLines(to_save)

  fileConn<-file(outfile)
  writeLines(to_save, fileConn)
  close(fileConn)

}
