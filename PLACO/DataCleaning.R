library(magrittr)
library(data.table)
library(R.utils)
library(dplyr)
f1 <- fread(file = "./IBS_sumstats_cleaned.txt.gz", header = T)
m1 <- fread(file = "./MDD/mddoriginal.gz", header = T)
m1 <- m1[!is.na(m1$SNP),]
attach(m1)
m1$EAF <- m1$FRQ_U_507679
m1 <- m1[which(EAF > 0.01 & EAF < 0.99),]
m1 <- m1[A1 == "A" | A1 == "T" | A1 == "C" | A1 == "G" ,]
m1 <- m1[A2 == "A" | A2 == "T" | A2 == "C" | A2 == "G" ,]
m1 <- m1[INFO > 0.8 ,]
detach(m1)
m1 <- transform(m1, BETA = log(OR), N = Nca + Nco)
m1 %<>%
  mutate(Z = BETA/SE)
names(m1)[names(m1) == "BP"] <- "POS"
m1 %<>%
  mutate(Case_Percentage = Nca/N)
m1 %<>%
  mutate(MLOG10P = -log10(P))
# Keep columns
mv <- c("SNP", "CHR", "POS", "A1", "A2", "EAF", "BETA", "SE", "P", "MLOG10P", "N", "Case_Percentage")
m1 <- m1[, mv, with = F]
# Export as txt.gz or txt format, quote=F, sep="\t", no row names, but include column names
write.table(m1, gzfile("./MDD/mdd_sumstats_cleaned.txt.gz"),
            quote = F, sep = "\t", row.names = F, col.names = T)