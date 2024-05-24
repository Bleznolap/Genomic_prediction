
#' @param pheno A phenotype dataframe output of filter.R function for the target population
#' @param ... A list of polygenic scores computed using Get_poly.R for each source population
#' @returns a phenotype dataframe containing the polygenic score(s) 

#' @details
#' This function adds the polygenic score(s) to the target population phenotype data

#Add poly

add_poly <- function(pheno, ...) {
  poly_vars <- list(...)
  for (i in seq_along(poly_vars)) {
    column_name <- paste0("poly", i)
    pheno[[column_name]] <- poly_vars[[i]]
  }
  return(pheno)
}
