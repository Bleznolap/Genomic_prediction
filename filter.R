
#Filter and match genotype, phenotype, and map

filter_and_match<-function(geno,pheno, map=NULL, n_row, n_col, n_maf){
  #NA per row
  row_na <-  rowSums(is.na(geno))
  geno <- geno[row_na<=n_row, ]
  #NA per column
  col_na <-  colSums(is.na(geno))
  geno <- geno[, col_na<=n_col]
  #MAF
  Maf <- colMeans(geno)/2
  Maf <- pmin(Maf, 1-Maf)
  geno <- geno[, Maf>= n_maf]
  #match geno and pheno by rownames
  Sub_com<- intersect(rownames(geno), rownames(pheno))
  geno <- geno[Sub_com, ]
  pheno <- pheno[Sub_com, ]
  if(!is.null(map)){
    snp_idx <- match(colnames(geno), rownames(map))
    map<- map[snp_idx, ]
  }
  #number of deleted rows and columns in each step of filter and match
  cat("-------------------------------","\n")
  cat(paste("number of deleted row:",sum(row_na>n_row),"\n"))
  cat(paste("number of deleted column:",sum(col_na>n_col),"\n"))
  cat(paste("number of deleted MAF:",sum(Maf<n_maf),"\n"))
  
  cat("-------------------------------","\n")
  cat(paste("Sample Output for genotype:","\n"))
  print(dim(geno))
  print(geno[1:5,1:5])
  cat(paste("Sample Output for phenotype:","\n"))
  print(pheno[1:5,])
  cat(paste("number of subjects for genotype:",nrow(geno),"\n"))
  cat(paste("number of subjects for phenotype:",nrow(pheno),"\n"))
  cat("-------------------------------","\n")
  if(!is.null(map)){
    print("A map file was provided!")
    print(head(map))
    print(table(map$chr))
  }
  return(list(geno=geno, pheno=pheno, map=map))
}