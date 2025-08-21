library(data.table)
library(dplyr)
library(parallel)
for (j in sumstats1_for) {
  for (k in sumstats2_for) {
    sumstats1=j
    sumstats2=k
    setwd(paste0(path_prefix,file_name,"/PLACO/",sumstats1,"-",sumstats2))
    load(paste0(sumstats1,"-",sumstats2,".RData"))
    ###设置多核并行跑
    cl <- makeCluster(ncore)
    clusterExport(cl,c("j","k","sumstats1","sumstats2"))
    clusterEvalQ(cl, {load(paste0(sumstats1,"-",sumstats2,".RData"))})
    out <- parSapply(cl,1:p, function(i) placo(Z=Z.matrix.decor[i,], VarZ=VarZ))
    stopCluster(cl)
    ####整合结果到原始文件
    f12<-fread(paste0(sumstats1, "-",sumstats2,".txt.gz"))
    out %>% 
      t() %>% 
      as.data.table() %>% 
      mutate(p.placo=as.numeric(p.placo)) %>% 
      mutate(T.placo=as.numeric(T.placo)) %>% 
      cbind(f12)->f12
    f12 %>%
      filter(P.f1<10^-4) %>% 
      filter(P.f2<10^-4) %>%   
      filter(p.placo<5*10^-8)->f12_sig###输出显著的点
    fwrite(f12,
           file = paste0(sumstats1, "-",sumstats2,".txt.gz"),
           compress = "gzip",
           sep = "\t",
           quote = F,
           col.names = T,
           row.names = F)
    write.table(f12_sig,
           file = paste0(sumstats1, "-",sumstats2,"_sig.txt"),
           sep = "\t",
           quote = F,
           col.names = T,
           row.names = F)

    rm(f12)
    rm(f12_sig)
    rm(out)
    gc()
    print(paste0(j,"-",k))    
    
  }
  
}