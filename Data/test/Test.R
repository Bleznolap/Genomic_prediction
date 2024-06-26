
library(data.table)
library(gwaR)
library(regress)
library(Matrix)
library(knitr)
library(magrittr)
source("Toolkit_function.R")

source_pop = 'marc'
pheno_name = pheno
geno_name = geno
cg_name2 = 'cg'
id_name2 = 'animal'
fe_cov2 = c('age')
fe_fct2 = c('sex')
phenotype = 'imf'
na_geno = NA
na_pheno = NA
na_per_row = 0.05
na_per_col = 0.05
MAF = 0.05

genot <-read_genotypes(geno_name)
phenot <-read_pheno(file_pheno = pheno_name, source_pop, id_name, fe_cov=fe_cov, 
                    cg_name=cg_name, phenotype = phenotype)
flt <- filter_and_match(genot, phenot, n_maf = MAF,
                        n_row=na_per_row, n_col=na_per_col)
mat <- read_G_matrix(flt$geno)
gb <- model_gb(G= mat, pheno=flt$pheno, y_name=phenotype,fixed=c(fe_cov),random=cg_name)
Gw<-run_gwa(gb, flt$geno)