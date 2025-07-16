source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/4.module_pathway/statistic")

level <- "pathway" # "pathway","module"
meta <- "high" # high, v1, v2, normal
scFEA_raw<-read_in_FC(level, meta)
scFEA <- ct_adj(scFEA_raw)

ct <- levels(as.factor(scFEA[,4]))
pw <- levels(as.factor(scFEA[,3]))
sam <- levels(as.factor(scFEA[,2]))
Nct <- as.numeric(table(scFEA$celltype)/length(pw))

# 在平均水平上/下数量统计
sumU <- matrix(0,nrow = length(ct),ncol = length(pw),dimnames = list(ct,pw))
sumD <- matrix(0,nrow = length(ct),ncol = length(pw),dimnames = list(ct,pw))
for (i in ct) {
  for (j in pw) {
    a <- scFEA[scFEA[,1] %in% j & scFEA$celltype %in% i,]
    sumU[i,j] <- length(which(a$log2FC>0))
    sumD[i,j] <- length(which(a$log2FC<0))
  }
}
sumU <- sumU/Nct*100
sumD <- sumD/Nct*100
# t(sumU)、t(sumD)，看上下调的通路/模块
write.table(sumU,file=paste0(level,"_mean_up_",meta,".txt"),sep="\t")
write.table(sumD,file=paste0(level,"_mean_down_",meta,".txt"),sep="\t") # log2FC < 0 in more than 50% of samples

re <- matrix(NA,nrow = length(ct),ncol = 7,dimnames = list(ct,c("# Total samples","80%","50%","Upregulated >= 80% samples",
             "Upregulated >= 50% samples","Downregulated >= 80% samples","Downregulated >= 50% samples")))
re[,1] <- Nct
re[,2] <- 0.5*Nct
re[,3] <- 0.8*Nct
for (r in 1:nrow(sumU)) {
  a <- length(which(sumU[r,]>50))
  ba <- round(as.numeric(a)/as.numeric(length(pw))*100,2)
  b <- length(which(sumU[r,]>80))
  bb <- round(as.numeric(b)/as.numeric(length(pw))*100,2)
  re[r,5] <- paste0(a,"(",ba,"%)")
  re[r,4] <- paste0(b,"(",bb,"%)")
}
for (r in 1:nrow(sumD)) {
  a <- length(which(sumD[r,]>50))
  ba <- round(as.numeric(a)/as.numeric(length(pw))*100,2)
  b <- length(which(sumD[r,]>80))
  bb <- round(as.numeric(b)/as.numeric(length(pw))*100,2)
  re[r,7] <- paste0(a,"(",ba,"%)")
  re[r,6] <- paste0(b,"(",bb,"%)")
}
write.table(re,paste0(level,"_up_down_statistic.txt"),sep = "\t",quote = F)


######
# 仅在某一细胞类型特异性在平均水平上/下
level <- "module" # "pathway","module"
meta <- "high" # high, v1, v2, normal
scFEA_raw<-read_in_FC(level, meta)
scFEA <- ct_adj(scFEA_raw)

ct <- levels(as.factor(scFEA[,4]))
pw <- levels(as.factor(scFEA[,3]))
sam <- levels(as.factor(scFEA[,2]))
Nct <- as.numeric(table(scFEA$celltype)/length(pw))

sumSu <- matrix(0,nrow = length(ct),ncol = length(pw),dimnames = list(ct,pw))
sumSd <- matrix(0,nrow = length(ct),ncol = length(pw),dimnames = list(ct,pw))
for (s in sam) {
  for (j in pw) {
    a <- scFEA[scFEA$module_1 %in% j & scFEA$sample %in% s,]
    if(length(which(a$log2FC>0))==1){
      i <- a$celltype[which(a$log2FC>0)]
      sumSu[i,j] <- sumSu[i,j]+1
    }
    if(length(which(a$log2FC<0)) == 1){
      i <- a$celltype[which(a$log2FC<0)]
      sumSd[i,j] <- sumSd[i,j]+1
    }
  }
}
sumSu <- sumSu/Nct*100
sumSd <- sumSd/Nct*100
for (r in 1:length(ct)) {
  print(paste0(ct[r],":",colnames(sumSu)[which(sumSu[ct[r],]>40)])) # log2FC > 0 exclusively in malignant cells in more than 50% of samples
  print(paste0(ct[r],":",colnames(sumSd)[which(sumSd[ct[r],]>40)]))
}
write.table(sumSu,file="exclusive_up.txt",sep="\t")
write.table(sumSd,file="exclusive_down.txt",sep="\t")
