######
# 代谢能力较为相似的细胞类型：主要考虑Macrophages, Fibroblasts, Dendritic cells和Malignant cells的相关性
source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/3.metabolic_landscape/overall")

library(Hmisc) # rcorr函数
library(corrplot)
library(RColorBrewer) # 配色包
library(paletteer) # 配色包
# 读入数据
level <- "module" # "pathway","module"
meta <- "high" # high, v1, v2, normal
scFEA_raw<-read_in_FC(level, meta)
scFEA <- ct_adj(scFEA_raw)
ct <- levels(as.factor(scFEA[,4]))
pw <- levels(as.factor(scFEA[,3]))
dim <- list(c(ct),c(pw))
re_ct <- matrix(nrow = length(ct),ncol = length(pw),dimnames = dim)
for (i in ct) {
  for (j in pw) {
    re_ct[i,j] <- mean(scFEA[which(scFEA[,4]==i & scFEA[,3]==j),5])
  }
}

# 建立用于计算相关性的矩阵
cor_ma  <- t(re_ct)

res_p <- rcorr(cor_ma, type = "pearson")
r1 <- res_p$r
r1[which(is.na(r1))] <- 1
r1[which(is.infinite(r1))] <- 0
r1 <- round(r1,2)
p1 <- res_p$P
p1[which(is.na(p1))] <- 0
res_s <- rcorr(cor_ma, type = "spearman")
r2 <- res_s$r
r2[which(is.na(r2))] <- 1
r2[which(is.infinite(r2))] <- 0
r2 <- round(r2,2)
p2 <- res_s$P
p2[which(is.na(p2))] <- 0
save(r1,p1,r2,p2,file="cor_metabolism.Rdata") 

# load("cor_metabolism.Rdata")
png(file=paste("corplot_celltype_module.png",sep = ""),width=850,height=570,res=96)
corrplot(r1, method = "color",
         order = "hclust", addrect = 6, 
         hclust.method = "mcquitty",
         tl.cex=0.7, tl.col="black",
         type = "full",
         diag = TRUE,
         cl.cex = 0.7,
         col = colorRampPalette(colors = c("#0da9ce","white","#f49537"))(100),
         p.mat = p1, insig = "label_sig", sig.level = c(0.01), 
         pch.cex = .7, pch.col = "white")
dev.off()
pdf("corplot_celltype_module.pdf",width=9,height=9)
corrplot(r1, method = "color",
         order = "hclust", addrect = 6, 
         hclust.method = "mcquitty",
         tl.cex=0.7, tl.col="black",
         type = "full",
         diag = TRUE,
         cl.cex = 0.7,
         col = colorRampPalette(colors = c("#0da9ce","white","#f49537"))(100),
         p.mat = p1, insig = "label_sig", sig.level = c(0.01), 
         pch.cex = .7, pch.col = "white")
dev.off()
# hclust.method = c("complete", "ward", "ward.D", "ward.D2", "single", "average",
#                   "mcquitty", "median", "centroid"),
# rev(paletteer_c("ggthemes::Red-Blue-White Diverging", 30))
# rev(paletteer_c("ggthemes::Classic Red-Blue", 30))

