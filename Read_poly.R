read_poly <-function(g, ghat, scale=TRUE){
  if(scale == TRUE){
    p <- colMeans(g)/2
    Z <-sweep(g,2,2*p)/sqrt(sum(2*p*(1-p)))
    poly_s <- Z%*%ghat
  } else{
    poly_s <- g%*%ghat
  }
  poly <- unname(poly_s)
  num <- nrow(g)[1] #Sample size in population
  cat(".........................................", "\n")
  cat("Estimated breeding value showing first 6 animals of", num, ":", "\n")
  print(head(poly))
  cat(".........................................", "\n")
  return(poly)
}