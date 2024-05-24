
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