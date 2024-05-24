

#Add poly

add_poly <- function(pheno, ...) {
  poly_vars <- list(...)
  for (i in seq_along(poly_vars)) {
    column_name <- paste0("poly", i)
    pheno[[column_name]] <- poly_vars[[i]]
  }
  return(pheno)
}