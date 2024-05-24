
#' @param g A genotype matrix n X m output of filter.R function for the target population
#' @param gwa A gwas object from gwa.R for the source population
#' @param scale A condition to standardize the genotype matrix (true by default)
#' @returns a matrix n x 1

#' @details
#' This function computes the polygenic score using the genotype matrix of the source population (receiving 
#' population) and the SNP effects of the target population (s) (contributing population(s)). If scale is FALSE,
#' the standardized genotype matrix (Z matrix) can replace genotype matrix.

#Get poly

get_polys<- function(g, gwa, scale=TRUE){
  if(scale == TRUE){
    p <- colMeans(g)/2
    Z <-sweep(g,2,2*p)/sqrt(sum(2*p*(1-p)))
    int <- intersect(colnames(Z),rownames(gwa))
    poly_s <- Z[, int]%*%gwa[int, 1]
  }
  else{
    int <- intersect(colnames(g),rownames(gwa))
    poly_s <- g[, int]%*%gwa[int, 1]
  }
  return(poly_s)
}
