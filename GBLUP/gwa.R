
#' @param gb A gblup object
#' @param g A genotype matrix n X m output of filter.R function
#' @param scale A condition to standardize the genotype matrix (true by default)
#' @returns an object of the class gwas: a dataframe with two columns and length equal to number of markers
#' \itemize{ 
#'    \item {\code{ghat}} {estimated SNP effects} 
#'    \item{\code{varg}} {estimated SNP effect variances}
#'}

#' @details
#' This function estimates the SNP effects and their variances. If scale is FALSE, the standardized genotype 
#' matrix (Z matrix) can replace the genotype matrix.

#Run GWA
run_gwa<- function(gb, g, scale=TRUE){
  if(scale == TRUE){
    p <- colMeans(g)/2
    Z <-sweep(g,2,2*p)/sqrt(sum(2*p*(1-p)))
    gw <- gwas(gb, t(Z))
  }
  else {
    gw <- gwas(gb, t(g))
  }
  cat(".........................................", "\n")
  cat("Estimated SNP Effect and its variance:", "\n")
  print(gw[1:5, ])
  cat(".........................................", "\n")
  cat("Summary of GWA:", "\n")
  print(summary(gw))
  cat(".........................................", "\n")
  return(gw)
}
