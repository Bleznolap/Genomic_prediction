
#Read Genotype

read_genotypes<-function(file){
  geno<-fread(file,header=TRUE,data.table=FALSE)
  rownames(geno)<-geno[,1]
  geno<-as.matrix(geno[,-1])
  cat("-------------------------------","\n")
  cat(paste("Read genotype file",file,"\n"))
  cat("sample output:","\n")
  print(geno[1:5,1:5])
  cat(paste("number of subjects:",nrow(geno),"\n"))
  cat(paste("number of SNP:",ncol(geno),"\n"))
  cat("-------------------------------","\n")
  
  return(geno)
}