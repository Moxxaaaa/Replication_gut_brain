#if(!require("remotes"))
#  install.packages("remotes")
#install.packages("dplyr")
#library(remotes)
#install_github("chr1swallace/coloc",build_vignettes=TRUE)

library(data.table)
library("coloc")
library(dplyr)
gwas <- fread("C:\\Users\\86132\\Desktop\\coloc\\coloc\\GWAS.txt")
eqtl <- fread("C:\\Users\\86132\\Desktop\\coloc\\coloc\\eQTL.txt")
input <- merge(eqtl, gwas, by="rs_id", all=FALSE, suffixes=c("_eqtl","_gwas"))
result <- coloc.abf(dataset1=list(snp=input$rs_id,pvalues=input$pval_nominal_gwas, type="cc", s=0.33, N=50000), 
                    dataset2=list(snp=input$rs_id,pvalues=input$pval_nominal_eqtl, type="quant", N=10000), 
                    MAF=input$maf)
need_result=result$results %>% filter(SNP.PP.H4 > 0.95)
a<-result$results