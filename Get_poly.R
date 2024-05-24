
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