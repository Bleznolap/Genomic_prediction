
#' @param geno A matrix n X m output of filter.R function 
#' @returns a matrix n x n for the genomic relationship matrix

#' @details
#' This function computes the genomic relationship matrix.

#Set G matrix

read_G_matrix <- function(g){
  p <- colMeans(g)/2
  Z <-sweep(g,2,2*p)/sqrt(sum(2*p*(1-p)))
  G <- tcrossprod(Z)
  cat("---------------------------------", "\n")
  cat(paste("Genomic relationship matrix","\n"))
  print(dim(G))
  cat("Summary for diagonal:", "\n")
  print(summary(diag(G)))
  cat("---------------------------------", "\n")
  cat("Summary for off-diagonal:", "\n")
  print(summary(G[upper.tri(G, diag = FALSE)]))
  cat("---------------------------------", "\n")
  cat(paste("Normalized Z matrix:","\n"))
  print(dim(Z))
  class(G) <- "matrix"
  return(G)
}
