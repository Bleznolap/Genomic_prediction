
#' @param Gmtr A matrix n x n output of Gmatrix.R function 
#' @param pheno A phenotype dataframe output of filter.R function containing number 
#' of individuals on rows and required phenotype variable(s) on columns 
#' @param fixed A fixed effects varaible
#' @param random A random effects variable 
#' @param y_name A response variable
#' @returns an object of class gblup

#' @details
#' This function computes genomic prediction

#Set GBLUP

model_gb<-function(Gmtr,pheno,y_name,fixed=NULL,random){
  if(!is.null(fixed)){
    model_spec<-c(formula(paste("y ~ ", paste(fixed,collapse = " + "))),
                  (formula(paste("y ~",random))))
    cat("..........................................................", "\n")
    gb<-gblup(rsp = y_name,data = pheno,
              design =  c(formula(paste(" ~ ", paste(fixed,collapse = " + "))),
                          (formula(paste("~", random)))),G = Gmtr)
    print(summary(gb, fe=TRUE))}
  else if(is.null(fixed)){
    model_spec<-c(formula(y~1), (formula(paste("y~",random))))
    cat("..........................................................", "\n")
    gb<-gblup(rsp = y_name,data = pheno,
              design =  c(formula(y~1), (formula(paste("~", random)))),G = Gmtr)
    print(summary(gb))}
  cat("..........................................................", "\n")
  return(gb)
}
