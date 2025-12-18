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

GRM_VR <- function(x){x1 <- as.matrix(x); p0 <- colMeans(x1)/2; maf <- pmin(p0,1-p0) >0.05; gen <- x1[,maf]; p1 <- colMeans(gen)/2; z0 <- sweep(gen, 2, 2*p1)/sqrt(sum(2*p1*(1-p1))); G0 <- tcrossprod(z0); return(list(GRM=G0, freq=p1, zmat=z0, freg0=p0))}
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


cppFunction(depends = "RcppArmadillo", '
arma::mat f_LS_beta( arma::mat X,
                     arma::mat y
) {
  arma::mat beta = X.t() * y * X;
  return beta;
}
')

cppFunction(depends = "RcppArmadillo", '
arma::vec Z_sq( arma::vec X,
                     arma::mat y
) {
  arma::vec Z = X.t() * y * X;
  return Z;
}
')

## Function for SNP correlation pruning for single SNP
snp_corr <- function(zm, t){
cor_matrix <- rcorr(zm)
r_abs <- round(abs(cor_matrix$r),2)
no_v <- names(which(r_abs[,1]>t))[-1]
r_df <- r_abs[!rownames(r_abs)%in%no_v,!rownames(r_abs)%in%no_v]
n=1
while(n>0){
lst <- sapply(1:ncol(r_df), function(x) length(which(r_df[,x]>t)))
n=which(unlist(lst)>1)

if(length(n)<1){
n <- 0
next
}
n <- head(n)[1]
no_v <- names(which(r_df[,n]>t))[-1]
r_df <- r_df[!rownames(r_df)%in%no_v,!rownames(r_df)%in%no_v]
}
return(r_final <- r_df)
}
 
## Function for GWAS based on correlation pruning
vg_gwas2 <- function(corr, grm, y){
col_snp <- colnames(corr)
grm <- gsm$GRM
zm <- gsm$zmat[,col_snp]
gp <- regress(y~1, ~grm)
Vin_e <- gp$W%*%(y-gp$fitted)
itx <- intersect(rownames(zm), rownames(Vin_e));
ghat <- gp$sigma[[1]]*crossprod(zm[itx,], Vin_e[itx,])
VQ <- gp$W %*% gp$Q %*% gp$Sigma %*%t(gp$Q) %*% gp$W
colnames(VQ) <- rownames(VQ)
tvq <- f_LS_beta(zm, VQ)
Vghat <- tvq*(gp[["sigma"]][[1]])^2
return(list(Vg=Vghat, gh=ghat))
}

## Necessary function for GWAS test
vg_gwas <- function(grm=NULL, y, gp, zm){
if(is.null(grm)){
gp <- gp
}else{
gp<- regress(y~1, ~grm)
}
Vin_e <- gp$W%*%(y-gp$fitted)
#itx <- intersect(rownames(zm), rownames(Vin_e));
#ghat <- gp$sigma[[1]]*crossprod(zm[itx,], Vin_e[itx,])
ghat <- gp$sigma[[1]]*crossprod(zm, Vin_e)
VQ <- gp$W %*% gp$Q %*% gp$Sigma %*%t(gp$Q) %*% gp$W
colnames(VQ) <- rownames(VQ)
tvq <- f_LS_beta(zm, VQ)
Vghat <- tvq*(gp[["sigma"]][[1]])^2
sol <- solve(Vghat)
C_f1 <- Z_sq(ghat, sol)
pval <- 1-pchisq(C_f1,nrow(sol))
return(list(Vg=Vghat, gh=ghat, sol=sol, pval=pval, chi=C_f1))
}

## SNP pruning by correlation thresholding (t-0.85) for multiple SNPs
part_corr <- function(zm,f,chk){
sz <- chk
tc <- ncol(zm)
st <- seq(1, tc, by=sz)
ed <- seq(sz, tc, by=sz)

if(ed[length(ed)]<tc){
ed <- c(ed, tc)
}##
corr_ls <- list()
for(i in 1:length(st)){
st_col <-  st[i]
ed_col <- ed[i]
zm_subt <- zm[, st_col:ed_col, drop=FALSE]
corr_ls[[i]] <- snp_corr(zm_subt, t=0.85)
}
assign(paste("cor_",f,sep=""), corr_ls)
save(list=paste("cor_",f,sep=""), file=paste("id_corr_ls_", f,".RData", sep=""))
}



## GWAS test
## y - phenotype; gp - gblup outputs; sm - standardize genotype; 
## snp - unlist SNP ID (SNP ID per gene); th - % cumulative variance explained threshold
GGwas2 <- function(y, gp, sm,snp, th){
GG_tst <- list()
GG_pcs <- list()
GG_var <- list()
for(i in seq_along(snp)){
Zm <- sm[,snp[[i]], drop=F]
pc <- prcomp(Zm)
px <- pc$x
pc_cumsum <- cumsum(pc$sdev^2)/sum(pc$sdev^2)
sel <- pc_cumsum[pc_cumsum > th][1]
#cat("Sel_pc_",i,":",sel,sep="", "\n\n") #Variance explained by the PCs that explained over 98.9 of the variation
n_pcs <- which(pc_cumsum==sel)
#cat("No_PCs_",i,":",n_pcs,sep="","\n\n") #No of PCs that explained over 98.9 of the variation
GG_var[[i]] <- round(sel,3)
GG_pcs[[i]] <- n_pcs[1]
GG_tst[[i]] <- vg_gwas(y=y,gp=gp,zm=px[,1:n_pcs, drop=F])
}
return(list(tst=GG_tst, pcs=GG_pcs, var=GG_var))
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

##########################################################
## Mapping SNPs to gene start and end positions ##
##########################################################

## output: an unlisted SNP ID per gene per chromosome

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
                                              
  n_gene <- list()
    for(i in seq_along(tot)){n_gene[[i]] = tot[[i]]}
  snpid <- unique(unlist(n_gene))
              
  n_gene_0 <- list()
    for(i in seq_along(tot_0)){n_gene_0[[i]] = tot_0[[i]]}
  snpid_0 <- unique(unlist(n_gene_0))
                                           
  return(list(tot_chr=tot,tot_chr_0=tot_0, Go_chr=Go_fd, Go_chr_0=Go_fd_0, snp_id = snpid, snp_id_0 = snpid_0))
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
                   
  cbn <- lapply(GG_snpsCHR, function(x) sapply(x, function(y) y$ensembl_gene_id))
  n_snp <- list()
    for(i in seq_along(cbn)){n_snp[[i]] <- cbn[[i]]}
  geneid <- unique(unlist(n_snp))
                                               
  return(list(n_gene_chr=n_gene, GG_snp_chr = GG_snpsCHR,gene_id=geneid))
}




 



