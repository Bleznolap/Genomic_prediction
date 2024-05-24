

#' @param file A marker file
#' @details
#' This function allows reading of the raw genotype file and converting it to 
#' a matrix format. It outputs the forst five rows and columns.

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

#' @param file A map file
#' @details
#' This function allows reading of the raw genotype file and converting it to 
#' a matrix format. It outputs the forst five rows and columns.

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

#Set G matrix
read_G_matrix <- function(g){
  p <- colMeans(g)/2
  Z <-sweep(g,2,2*p)/sqrt(sum(2*p*(1-p)))
  G <- tcrossprod(Z)
  cat("---------------------------------", "\n")
  cat(paste("Genomic relationship matrix","\n"))
  print(dim(G))
  cat("Summary for diagonal:", "\n")
  print(summary(diag(G)))
  cat("---------------------------------", "\n")
  cat("Summary for off-diagonal:", "\n")
  print(summary(G[upper.tri(G, diag = FALSE)]))
  cat("---------------------------------", "\n")
  cat(paste("Normalized Z matrix:","\n"))
  print(dim(Z))
  class(G) <- "matrix"
  return(G)
}

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

#Run GWA
run_gwa<- function(gb, g, scale=TRUE){
  if(scale){
    p <- colMeans(g)/2
    Z <-sweep(g,2,2*p)/sqrt(sum(2*p*(1-p)))
  }
  gw <- gwas(gb, t(Z))
  cat(".........................................", "\n")
  cat("Estimated SNP Effect and its variance:", "\n")
  print(gw[1:5, ])
  cat(".........................................", "\n")
  cat("Summary of GWA:", "\n")
  print(summary(gw))
  cat(".........................................", "\n")
  return(gw)
}

get_polys<- function(g, gwa, scale=TRUE){
  if(scale){
    p <- colMeans(g)/2
    Z <-sweep(g,2,2*p)/sqrt(sum(2*p*(1-p)))
  }
  int <- intersect(colnames(Z),rownames(gwa))
  poly_s <- Z[, int]%*%gwa[int, 1]
  return(poly_s)
}

add_poly <- function(pheno, ...) {
  poly_vars <- list(...)
  for (i in seq_along(poly_vars)) {
    column_name <- paste0("poly", i)
    pheno[[column_name]] <- poly_vars[[i]]
  }
  return(pheno)
}

#Cross Validation
cross_validation<- function(G,pheno,iterate=5,y_name,fixed=NULL,rand_var, poly_name=NULL, pop=NULL, seed=1900){
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

read_poly <-function(g, ghat, scale=TRUE){
  if(scale){
    p <- colMeans(g)/2
    Z <-sweep(g,2,2*p)/sqrt(sum(2*p*(1-p)))
  }

  poly_s <- Z%*%ghat
  poly <- unname(poly_s)
  num <- nrow(g)[1] #Sample size in population
  cat(".........................................", "\n")
  cat("Estimated breeding value showing first 6 animals of", num, ":", "\n")
  print(head(poly))
  cat(".........................................", "\n")
  return(poly)
}
  
