rm(list = ls()) 
gc()
setwd("D:/Rworkplace/Example/COLOC")
library(coloc)
library(data.table)
library(tidyverse)
library(magrittr)
library(susieR)
IBS_MD<-fread(".\\IBS_MD.txt")
IBS_MD %>% 
  filter(EAf_IBS>=0.5) %>% 
  mutate(MAF=1-EAf_IBS) %>% 
  mutate(nbeta_IBS=beta_IBS/-1) %>% 
  mutate(nbeta_MD=beta_MD/-1) %>%
  mutate(nZ_IBS=Z_IBS/-1) %>% 
  mutate(nZ_MD=Z_MD/-1)->f1

IBS_MD %>% 
  filter(EAf_IBS<0.5) %>% 
  mutate(MAF=EAf_IBS) %>% 
  mutate(nbeta_IBS=beta_IBS) %>% 
  mutate(nbeta_MD=beta_MD) %>%
  mutate(nZ_IBS=Z_IBS) %>% 
  mutate(nZ_MD=Z_MD)->f2

IBS_MD <-rbind(f1,f2)
rm(f1);rm(f2);gc()
IBS_MD %<>%select(-c("EAf_IBS","EAf_MD","beta_MD","Z_MD","beta_IBS","Z_IBS"))

IBS_MD %>% 
  filter(chr==1) %>% 
  filter(pos>=175902660&pos<=176406835) %>%
  mutate(varbeta_IBS=se_IBS^2) %>% 
  mutate(varbeta_MD=se_MD^2) %>% 
  rename(A1=EA,A2=OA) %>% 
  select(rsid,chr,pos,A1,A2,
         nbeta_IBS,varbeta_IBS,
         nbeta_MD,varbeta_MD,
         N_IBS,N_MD,
         P_IBS,P_MD,varbeta_IBS,varbeta_MD,
         MAF)->region1

dataset1=list(snp=region1$rsid,
              beta=region1$nbeta_IBS, 
              varbeta=region1$varbeta_IBS,
              type="cc", 
              #s=0.1097408, 
              N=max(region1$N_IBS))



dataset2=list(snp=region1$rsid,
              beta=region1$nbeta_MD, 
              varbeta=region1$varbeta_MD,
              type="cc",  
              #s=0.3413761,
              N=max(region1$N_MD))

check_dataset(dataset1)
check_dataset(dataset2)
result <- coloc.abf(dataset1, 
                    dataset2, 
                    MAF=region1$MAF_MD)
result
a<-result$results


#------------susieR---------------#

#plink --bfile ref_panel --extract snplist.txt --r --ld-window-r2 0 --out ld_output

snplist<-region1[,1]
write.table(snplist,
            file = "snplist.txt",
            sep = "",
            col.names =F,
            row.names = F,
            quote = F)



ref_panel="./Ref_Panel/g1000_eur/g1000_eur"
snplist_dir="./Example/COLOC/snplist.txt"
ld_output="./Example/COLOC/example"
command<-paste0(
  "plink --bfile ",ref_panel,
  " --extract ",snplist_dir,
  " --r --ld-window-kb 3000 --ld-window 99999 --ld-window-r2 0  --out ",ld_output)
#--ld-window-kb 3000 --ld-window 99999 --ld-window-r2 0
# --r square  --out 
system(command)


f1="D:/Rworkplace/Example/COLOC/example.ld"
f1<-fread(f1)

f1 %>% distinct(SNP_A) %>% dplyr::rename(rsid=SNP_A)->c1
f1 %>% distinct(SNP_B) %>% dplyr::rename(rsid=SNP_B)->c2
f1<-as.data.frame(f1)
c<-bind_rows(c1,c2)
c %<>% distinct(rsid)
region1 %<>% inner_join(c,by="rsid") 

snplist<-region1$rsid
snp<-region1[,"rsid"]
snp %<>% rename(SNP_B=rsid)
snp$SNP_A=NA
result<-list()
for (i in 1:length(snplist)) {
  snp_name=snplist[i]
  snp %<>% mutate(SNP_A=snp_name)
  result[[i]]=snp
  print(i)
}
result<-bind_rows(result)
result %<>% mutate(id=paste0(SNP_A,":",SNP_B))
f1 %>% 
  rename(V1=SNP_A,V2=SNP_B) %>% 
  rename(SNP_A=V2,SNP_B=V1) %>% 
  select("CHR_A","BP_A","SNP_A","CHR_B","BP_B","SNP_B","R") ->f2

f<-rbind(f1,f2)
f %<>% mutate(id=paste0(SNP_A,":",SNP_B)) %>% select(id,R)
f3<-data.frame(snplist,snplist)
f3 %<>%
  mutate(id=paste0(snplist,":",snplist.1)) %>% 
  mutate(R=1) %>% 
  select(id,R)
f<-rbind(f,f3)

result_all<-left_join(result,f,by="id")
R_c=result_all$R



Nsnps=nrow(region1)

LD_matrix <- matrix(R_c,nrow = Nsnps, ncol = Nsnps,byrow = TRUE)


rownames(LD_matrix)<-snplist;colnames(LD_matrix)<-snplist

dataset.1=list(snp=region1$rsid,
              beta=region1$nbeta_IBS, 
              varbeta=region1$varbeta_IBS,
              type="cc", 
              #s=0.1097408, 
              N=max(region1$N_IBS),
              LD=LD_matrix,
              pos=region1$pos)



dataset.2=list(snp=region1$rsid,
              beta=region1$nbeta_MD, 
              varbeta=region1$varbeta_MD,
              type="cc",  
              #s=0.3413761,
              N=max(region1$N_MD),
              LD=LD_matrix,
              pos=region1$pos)

check_dataset(dataset.1)
check_dataset(dataset.2)

result1<-runsusie(dataset.1)
result2<-runsusie(dataset.2)
result<-coloc.susie(result1,result2)
summary(result1)
summary(result2)
print(result$summary)

if(requireNamespace("susieR",quietly=TRUE)) {
  susie.res=coloc.susie(dataset1,dataset2)
  print(susie.res$summary)
}

if(requireNamespace("susieR",quietly=TRUE)) {
  sensitivity(susie.res,"H4 > 0.9",row=1,dataset1=dataset1,dataset2=dataset2)
  sensitivity(susie.res,"H4 > 0.9",row=2,dataset1=dataset1,dataset2=dataset2)
}
