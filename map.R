
#Read Map

read_map <- function(file_map){
  map_n <- fread(file_map, data.table = FALSE)
  rownames(map_n)<- map_n[, "snp"]
  map_n<- map_n[, c("chr","pos")]
  cat("..................................","\n")
  cat("Summary of map:", "\n")
  print(head(map_n))
  print(table(map_n$chr))
  cat(".................................", "\n")
  return(map_n)
}