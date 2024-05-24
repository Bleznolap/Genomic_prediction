
#Set GBLUP

model_gb<-function(G_mtr,pheno,y_name,fixed=NULL,random){
  if(!is.null(fixed)){
    model_spec<-c(formula(paste("y ~ ", paste(fixed,collapse = " + "))),
                  (formula(paste("y ~",random))))
    cat("..........................................................", "\n")
    gb<-gblup(rsp = y_name,data = pheno,
              design =  c(formula(paste(" ~ ", paste(fixed,collapse = " + "))),
                          (formula(paste("~", random)))),G = G_mtr)
    print(summary(gb, fe=TRUE))}
  else if(is.null(fixed)){
    model_spec<-c(formula(y~1), (formula(paste("y~",random))))
    cat("..........................................................", "\n")
    gb<-gblup(rsp = y_name,data = pheno,
              design =  c(formula(y~1), (formula(paste("~", random)))),G = G_mtr)
    print(summary(gb))}
  cat("..........................................................", "\n")
  return(gb)
}