#----------------------------------1.设置表型以及生成Rdata---------------------------------
rm(list = ls()) 
gc()
####-----------------设置表型---------###
library(readxl)
code<-read_xlsx("D:/研究生/202309/疾病编号【230729】.xlsx",sheet = "EUR")
code<-as.data.frame(code)
#匹配的是code的表型那一列！！！！
sumstats1_for=c("PGC102_no23andMe")
sumstats2_for=c("IBS")
#--------设置参数-------#
path_prefix="D:/研究生/202309/0918验证PLACO/" ###这路径要存在
file_name="IBS_MDD"###会创造这个文件夹，PLACO结果放在这个文件夹的PLACO子文件夹内
PLACO_path="D:/研究生/202309/0918验证PLACO/PLACO.RData"###PLACO包
#------运行代码-----#
source("D:/研究生/202309/0918验证PLACO/1.PLACO生成Rdata.R")

###---------------------PLACO循环--------------------------
ncore=4#并行几核
source("D:/研究生/202309/0918验证PLACO/2.PLACO计算循环.R")
