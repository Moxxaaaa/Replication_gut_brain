library(data.table)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(readxl)
###创建文件夹
if(!dir.exists(paste0(path_prefix,file_name))){
  dir.create(paste0(path_prefix,file_name))
} 
#####文件夹创建PLACO子文件夹
if(!dir.exists(paste0(path_prefix,file_name,"\\PLACO"))){
  dir.create(paste0(path_prefix,file_name,"\\PLACO"))
} 

####循环#####
for (m in sumstats1_for) {
  for (n in sumstats2_for) {
    load(PLACO_path)#####载入PLACO
    sumstats1=m
    sumstats2=n
    if(!dir.exists(paste0(path_prefix,file_name,"/PLACO/",sumstats1,"-",sumstats2))){
      dir.create(paste0(path_prefix,file_name,"/PLACO/",sumstats1,"-",sumstats2))
    }####创建子文件夹
    setwd(paste0(path_prefix,file_name,"/PLACO/",sumstats1,"-",sumstats2))
    
    f1<-code[code$表型==m,"文件位置"]####获取文件路径
    f2<-code[code$表型==n,"文件位置"]####获取文件路径
    f1<-fread(f1)
    f2<-fread(f2)
    f12<-inner_join(f1,f2,by=c("SNP","CHR","POS","A1","A2"),suffix=c(".f1", ".f2"))
    rm(f1)
    rm(f2)
    gc()
    
    ###增加Z_col,并剔除Z^2>80
    f12 %<>%
      mutate(Z.f1=BETA.f1/SE.f1) %>% 
      mutate(Z.f2=BETA.f2/SE.f2) %>%   
      mutate(Z2.f1=Z.f1^2) %>% 
      mutate(Z2.f2=Z.f2^2) %>%  
      filter(Z2.f1<80) %>% 
      filter(Z2.f2<80)
    
    ###剔除MHC(chr6 25000000-35000000)
    f12 %<>% 
      filter(!(CHR==6&POS>25000000&POS<35000000)) %>% 
      select(SNP,CHR,POS,A1,A2,
             BETA.f1,Z.f1,P.f1,N.f1,EAF.f1,
             BETA.f2,Z.f2,P.f2,N.f2,EAF.f1)
    fwrite(f12,
           file = paste0(sumstats1, "-",sumstats2,".txt.gz"),
           compress = "gzip",
           sep = "\t",
           quote = F,
           col.names = T,
           row.names = F)
    
    ####PLACO#####
    k <- 2 
    p <-nrow(f12)
    Z.matrix<-as.matrix(f12[,c("Z.f1","Z.f2")])
    P.matrix<-as.matrix(f12[,c("P.f1","P.f2")])
    colnames(Z.matrix) <- paste("Z",1:k,sep="")
    colnames(P.matrix) <- paste("P",1:k,sep="")
    
    
    #### decorrelating the Z-scores
    R <- cor.pearson(Z.matrix, P.matrix, p.threshold=1e-4)###计算Z的相关度
    "%^%" <- function(x, pow)
      with(eigen(x), vectors %*% (values^pow * t(vectors)))####创建公式
    Z.matrix.decor <- Z.matrix %*% (R %^% (-0.5)) ###decor
    colnames(Z.matrix.decor) <- paste("Z",1:k,sep="")
    
    VarZ <- var.placo(Z.matrix.decor, P.matrix, p.threshold=1e-4)
    
    ####删除不相关的
    rm(f12)
    rm(Z.matrix)
    rm(P.matrix)
    rm(R)
    rm(k)
    rm("%^%")
    rm(cor.pearson)
    rm(var.placo)
    gc()
    save.image(file =paste0(sumstats1,"-",sumstats2,".RData"))
    print(paste(m,n,sep = "-"))
    
  }
}