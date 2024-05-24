

#' @param G A matrix n X n 
#' @param pheno A phenotype dataframe 
#' @param iterate A numeric object for the number of iterations (set as 5 by default)
#' @param y_name A response variable
#' @param fixed A fixed variable
#' @param rand_var A random variable
#' @param poly_name A list object for the number of polygenic scores (e.g c("poly1", "poly2", ...)) [set as NULL by default]
#' @param seed A numeric seed object to ensure reproducibility of results (set as 1900 by default)
#'
#' @returns a matrix of n x 1,  where n equals the number of iterations
#' @output a table with a single value representing the colmeans of the number of iterations


#' @details
#' This function computes the cross-validation for the prediction accuracy as the correlation between 
#' the true phenotype value and predicted breeding value
#' 

#Cross Validation
cross_validation<- function(G,pheno,iterate=5,y_name,fixed=NULL,rand_var, poly_name=NULL, seed=1900){
  set.seed(seed)
  d_pheno <- pheno
  obs <- nrow(d_pheno)
  tr <- round(obs/10)
  if(is.null(pop)){
    nc=1
  } else {
    nc=length(table(pop))
  }
  G_cross <- matrix(nrow = iterate, ncol = nc, NA)
  for (i in 1:iterate) {
    tc <- sample(1:obs, size = tr)
    d_pheno <- pheno
    d_pheno[,y_name][tc] <- NA
    if(!is.null(fixed)){
      gb_cross<-gblup(rsp = y_name,data = d_pheno,
                      design =  c(formula(paste(" ~ ", paste(fixed,collapse = " + "))),
                                  (formula(paste("~",rand_var)))),G = G)}
    if(is.null(fixed)){
      
      gb_cross<-gblup(rsp = y_name,data = d_pheno,
                      design =  c(formula(y~1), (formula(paste("~",rand_var)))),G = G)
    }
    
    A11 <- gb_cross$model$G
    A21 <- G[rownames(d_pheno)[tc], colnames(A11)]
    Uhat<- summary(gb_cross)$uhat
    Uhat2<- A21%*%chol2inv(chol(A11))%*%Uhat
    
    y_hat_m<-0
    
    if(!is.null(poly_name)){
      fe<-gb_cross$coefm[poly_name,1]
      new_X <- pheno[, poly_name]
      rownames(new_X)<-rownames(pheno)
      y_hat_m<-as.matrix(new_X[tc,])%*%fe
    }
    
    y_hat<-y_hat_m+Uhat2
    if(is.null(pop)){
      G_cross[i, 1] <- cor(y_hat,pheno[tc,gb_cross$name], use = "complete.obs")
    } else{
      mt<-cbind(y_hat,pheno[tc,gb_cross$name])
      id<-match(rownames(mt),names(pop))
      
      cr<-as.vector(by(mt,pop[id],function(x) cor(x)[2,1]))
      G_cross[i, ] <-cr
    }
    
  }
  print(kable(colMeans(na.omit(G_cross)), digits = 3, caption = 'Cor_G'))
  return(G_cross)
}

