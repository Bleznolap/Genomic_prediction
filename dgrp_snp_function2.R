pkgs <- c("tidyverse", "readxl", "data.table","regress","knitr",
            "kableExtra", "magrittr", "ggplot2","GGally","ggcorrplot",
            "ggpubr", "flextable", "Metrics", "cowplot", "skimr", "caret",
            "devtools", "readr", "Matrix", "gwaR", "rrBLUP", "flextable",
            "Rcpp","qqman","Hmisc","Matrix", "caret")

for (i in seq_along(pkgs)){
    pkg <- as.character(pkgs[[i]])
    tryCatch(if(pkg%in%rownames(installed.packages())){
        #cat("Available:", pkg, "\n")
        suppressPackageStartupMessages(library(pkg, character.only = TRUE))
        #cat("Package loaded:", pkg,"\n")
    }
    else{
        #cat("Not available:",pkg,">> now installing...", "\n")
        install.packages(pkg, repos =  "https://cloud.r-project.org")
        cat("Now available:", pkg, "\n")
        suppressPackageStartupMessages(library(pkg, character.only = TRUE))
        #cat("Package loaded:", pkg,"\n")
    }
    )
}

version

## Function to compute GRM
GRM <- function(x){x1 <- as.matrix(x); p0 <- colMeans(x1)/2; maf <- pmin(p0,1-p0) >0.05; gen <- x1[,maf]; p1 <- colMeans(gen)/2; z0 <- sweep(gen, 2, 2*p1)/sqrt(sum(4*p1*(1-p1))); G0 <- tcrossprod(z0); return(list(GRM=G0, freq=p1, zmat=z0, freg0=p0))}

## Function to compute marker effects and their variances
Ghat <- function(pheno, gp, zm){
  Vinve <- gp$W%*%(pheno-gp$fitted);
  itx <- intersect(rownames(zm), rownames(Vinve));
  ghat <- gp$sigma[[1]]*crossprod(zm[itx,], Vinve[itx,])
  
  # Calculate Vinv Q V t(Q) Vinv
  VQ <- gp$W %*% gp$Q %*% gp$Sigma %*%t(gp$Q) %*% gp$W;
  colnames(VQ) <- rownames(VQ)
    
  # apply the function cros prod
  tZVQ <- diagprod(t(zm[itx,]), VQ[itx, itx]);
    
  # Calculate the Variance of the G hat
  var_ghat <- tZVQ * (gp[["sigma"]][[1]])^2;
  ghat_var <- data.frame(ghat = ghat, var_ghat = var_ghat);
return(ghat_var)
}



## PCs grid search
PCs <- function(y, gp, sm,snp,t){
GG_pcs <- list()
GG_var <- list()
for(i in seq_along(snp)){
Zm <- sm[,snp[[i]], drop=F]
pc <- prcomp(Zm)
px <- pc$x
pc_cumsum <- cumsum(pc$sdev^2)/sum(pc$sdev^2)
sel <- pc_cumsum[pc_cumsum > t][1]
cat("Sel_pc_",i,":",sel,sep="", "\n\n") #Variance explained by the PCs that explained over t of the variation
n_pcs <- which(pc_cumsum==sel)
cat("No_PCs_",i,":",n_pcs,sep="","\n\n") #No of PCs that explained over t of the variation
GG_var[[i]] <- round(sel,3)
GG_pcs[[i]] <- n_pcs
}
return(list(pcs=GG_pcs, var=GG_var))
}

snp_map <- function(gene_df, snp_df, chr_lst){
  Go_fd <- list()
  Go_fd_0 <- list()
  for(j in seq_along(chr_lst)){
    chR <- gene_df %>% filter(chromosome_name==chr_lst[[j]])
    snP <- snp_df %>% filter(chr==chr_lst[[j]])
    GG_chr <- list()
    for(i in 1:nrow(chR)){
      GG_chr[[i]] <- snP[snP$pos>=chR$start_position[i] & 
                           snP$pos<=chR$end_position[i],]
    }
    Go_fd[[j]] <- GG_chr
    Go_fd_0[[j]] <- base::Filter(function(x) nrow(x)>0, GG_chr)
    
  }
  tot <- lapply(Go_fd, function(x) sapply(x, function(y) y$snp))
  tot_0 <- lapply(Go_fd_0, function(x) sapply(x, function(y) y$snp))
  #funct_snp <- unlist(c(tot[[1]], tot[[2]],tot[[3]], tot[[4]], tot[[5]]))
  return(list(tot_chr=tot,tot_chr_0=tot_0, Go_chr=Go_fd, Go_chr_0=Go_fd_0))
}


##########################################################
## Mapping gene to SNPs position ##
##########################################################

## output: an unlisted gene ID per SNP per chromosome
gene_map <- function(gene_df,snp_lst, chr_lst){
  GG_snpsCHR <- list()
  GG_snpsCHR_0 <- list()
  for(j in seq_along(chr_lst)){
    chR <- gene_df |> 
      filter(chromosome_name==chr_lst[[j]])
    snP <- snp_lst |> 
      filter(chr==chr_lst[[j]])
    GG_chr2 <- list()
    for(i in 1:nrow(snP)){
      GG_chr2[[i]] <- chR[chR$start_position<=snP$pos[i] & chR$end_position>=snP$pos[i],]
    }
    GG_snpsCHR[[j]] <- GG_chr2
  }
  n_gene <- lapply(GG_snpsCHR, function(x) sapply(x, nrow)) 
 
  return(list(n_gene_chr=n_gene, GG_snp_chr = GG_snpsCHR))
}



 
