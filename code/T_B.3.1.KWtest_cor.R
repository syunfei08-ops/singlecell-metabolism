######
# 一、表示各细胞类型趋势一致的模块
# 读入数据
source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/2-2.tumor/KWtest")

library(reshape2)
library(Hmisc)

level <- "module" # "pathway","module"
meta <- "high" # high, v1, v2, normal
scFEA_raw<-read_in_FC(level, meta)
scFEA <- ct_adj(scFEA_raw)
kw <- as.matrix(read.csv(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/4.module_pathway/KWtest/KWtest_flux.csv"))) 
row.names(kw) <- kw[,1]
kw <- kw[,-1,drop=F]
kw[is.na(kw)] <- 1

m <- "M_20" # M_20, M_77
for(m in levels(as.factor(scFEA$module_1))){

k <- as.numeric(kw[,m])
names(k) <- rownames(kw)

scFEA_M <- scFEA[which(scFEA$module_1==m & scFEA$celltype=="Macrophages"),c("sample","log2FC"),drop=F]
rownames(scFEA_M) <- scFEA_M[,1]
scFEA_M <- scFEA_M[,-1,drop=F]
scFEA_MC <- scFEA[which(scFEA$module_1==m & scFEA$celltype=="Malignant cells"),c("sample","log2FC"),drop=F]
rownames(scFEA_MC) <- scFEA_MC[,1]
scFEA_MC <- scFEA_MC[,-1,drop=F]

matr <- data.frame("conserve"=k)
matr <- merge(matr,scFEA_M,by="row.names",all=F)
rownames(matr) <- matr[,1]
matr <- matr[,-1]
matr <- merge(matr,scFEA_MC,by="row.names",all=F)
rownames(matr) <- matr[,1]
matr <- matr[,-1]
colnames(matr) <- c("conserve","Macrophages","Malignant cells")

cor_ma <- as.matrix(matr)
cor_ma <- apply(cor_ma,2,as.numeric)
reP <- rcorr(cor_ma, type = "pearson")
rP <- reP$r
pP <- reP$P
reS <- rcorr(cor_ma, type = "spearman")
rS <- reS$r
pS <- reS$P
if(abs(rS[2,1])>0.5 || abs(rS[3,1])>0.5){
  print(paste0(m,": "))
  print(rS)
}
}

# 比较两组中malignant cells的代谢强弱
source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/2-2.tumor/KWtest")

g <- 6
cho_ct <- "Malignant cells" # "Macrophages", "Malignant cells"
meta_data <- build_meta(cho_ct) # 1:sample, 2:celltype, 3:tumor, 4:primary_site, 5:tissue, 6:group1, 7:group2
group <- levels(as.factor(meta_data[,as.numeric(g)]))

level <- "module" # "pathway","module"
meta <- "high" # high, v1, v2, normal
scFEA_raw<-read_in_flux(level, meta)
scFEA <- ct_adj(scFEA_raw)

sam1 <- meta_data[which(meta_data[,as.numeric(g)]==group[1]),"sample",drop=T]
scFEA1 <- scFEA[scFEA$sample %in% sam1,]
scFEA1 <- scFEA1[which(scFEA1$celltype==cho_ct),c("sample","flux"),drop=F]
m1 <- aggregate(flux~sample,scFEA1,mean)
m1$group <- "group1"

sam2 <- meta_data[which(meta_data[,as.numeric(g)]==group[2]),"sample",drop=T]
scFEA2 <- scFEA[scFEA$sample %in% sam2,]
scFEA2 <- scFEA2[which(scFEA2$celltype==cho_ct),c("sample","flux"),drop=F]
m2 <- aggregate(flux~sample,scFEA2,mean)
m2$group <- "group2"

sam3 <- meta_data[which(meta_data[,as.numeric(g)]==group[3]),"sample",drop=T]
scFEA3 <- scFEA[scFEA$sample %in% sam3,]
scFEA3 <- scFEA3[which(scFEA3$celltype==cho_ct),c("sample","flux"),drop=F]
m3 <- aggregate(flux~sample,scFEA3,mean)
m3$group <- "group3"

matr <- rbind(m1,m2)
t.test(flux~group, matr, paired = FALSE, alternative = 'two.sided')
