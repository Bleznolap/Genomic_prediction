
#' @param file_pheno A pheno file
#' @param source_pop A population title
#' @param id_name A ID for the population e.g 'animal'
#' @param fe_cov A fixed covariate variable from the phenotype file (may not be included in data)
#' @param fe_fct A fixed factor variable from the phenotype file (may not be included in data)
#' @param cg_name A random variable from the phenotype file
#' @param phenotype A response variable

#' @details
#' This function allows reading of the raw phenotype file, selecting relevant variables,
#' and makes it available for use.

#Read Phenotype

read_pheno<-function(file_pheno, source_pop, id_name,fe_cov=NULL, fe_fct=NULL,cg_name,phenotype=NULL){
  pheno<-fread(file_pheno,data.table=FALSE)
  pheno$source=source_pop
  r_p<-pheno[,id_name]
  pheno<-pheno[,c(fe_cov,fe_fct,cg_name,phenotype)]
  pheno[,c(fe_fct,cg_name)]<-apply(pheno[,c(fe_fct,cg_name), drop=FALSE],2,function(x) factor(x))
  pheno<- as.data.frame(unclass(pheno), stringsAsFactors = TRUE)
  rownames(pheno)<- r_p
  cat("-------------------------------","\n")
  cat(paste("Read phenotype file",file_pheno,"\n"))
  cat("sample output:","\n")
  print(head(pheno))
  cat(paste("number of subjects:",nrow(pheno),"\n"))
  cat("deleting rows with missing values","\n")
  cat("-------------------------------","\n")
  pheno<-na.omit(pheno)
  cat(paste("number of subjects with complete observations:",nrow(pheno),"\n"))
  cat("-------------------------------","\n")
  cat("Phenotypic summary","\n")
  print(summary(pheno))
  cat("-------------------------------","\n")
  return(pheno)
}
