# 检查转录组和代谢组的一致性
source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
# memory.limit(10240)
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/3.metabolic_landscape/transcriptome")

library(Hmisc) # rcorr函数
library(corrplot)
library(RColorBrewer) # 配色包
library(paletteer) # 配色包

level <- "module" # "pathway","sub_pathway","module"
meta <- "high" # high, v1, v2, normal
A<-read_in_FC(level, meta)
ct <- sort(unique(A$celltype)) # 细胞类型统计
# 细胞类型的样本数目
nSamCt <- c()
for (i in ct) {
  ct_sam <- A[which(A$celltype==i),]
  nSamCt <- c(nSamCt, length(unique(ct_sam$sample)))
}
names(nSamCt) <- ct
re <- ct[which(as.numeric(nSamCt) < 5)]

a <- read.table(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/fluxandcelltype_",meta,".list"))[,1,drop=T]
n <- length(a)

len <- read.table("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/gene_length.txt")
rownames(len) <- len[,1]
len <- metabolics_gene(len)[, 2, drop=F]
colnames(len) <- "length"

for (i in 1:n){
  sam <- unlist(strsplit(a[i], ".", fixed = TRUE))[1]
  
  counts <- as.matrix(read.table(paste("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/count/",sam,".counts.csv",sep=""), header = T, sep = ","))
  colnames(counts) <- gsub('.',"-",colnames(counts),fixed = T)
  colnames(counts) <- gsub(':',"-",colnames(counts),fixed = T)
  colnames(counts) <- paste(sam,colnames(counts),sep="-")
  counts <- counts[rownames(counts) %in% rownames(len),]
  # 转换TPM
  C <- merge(counts,len,by="row.names",all=F)
  kb <- C$length / 1000
  rpk <- C[,2:(ncol(C)-1)]
  rpk <- rpk / kb
  tpm <- t(t(rpk)/colSums(rpk)*1000000)
  rownames(tpm) <- C$Row.names
  counts <- tpm
  
  metadata <- read.table(paste("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/celltype/",sam,".celltype",sep=""), header = F, sep = "\t")
  metadata[,1] <- gsub('.',"-",metadata[,1],fixed = T)
  metadata[,1] <- gsub(':',"-",metadata[,1],fixed = T)
  metadata[,1] <- paste(sam,metadata[,1],sep="-")
  rownames(metadata) <- metadata[,1]
  metadata <- metadata[,-1,drop=F]
  colnames(metadata) <- "celltype"
  metadata$sample <- sam
  
  C <- merge(t(counts),metadata,by="row.names",all=F)
  C <- C[which(!C$celltype %in% c(re,"Unassigned")),]
  cName <- C[,1,drop=T]
  counts <- t(C[,2:(ncol(C)-2)])
  # counts <- metabolics_gene(counts)
  colnames(counts) <- cName
  counts <- as.data.frame(counts)
  counts$gene_symbol <- rownames(counts)
  metadata <- metadata[which(!metadata$celltype %in% c(re,"Unassigned","celltype")),]
  
  if(i==1){
    allSample <- counts
    allMeta <- metadata
  }else{
    allSample <- merge(allSample,counts,by="gene_symbol",all=TRUE) # all参数控制是否求交集
    # allSample <- merge(allSample,counts,by="gene_symbol",all=FALSE)
    allMeta <- rbind(allMeta,metadata)
  }
  print(paste0(i,": ",nrow(metadata)))
}
save(allSample,allMeta,file="allData_RNA.Rdata")
# r<-rownames(allMeta[which(allMeta$sample=="NB-033-15-1A"),])
# r1<-intersect(r,colnames(allSample))
# r2<-setdiff(r,colnames(allSample))
# r2<-gsub("^([A-Z]+)-1$","\\1-1.x",r2)
# r2<-gsub("^([A-Z]+)-11$","\\1-1.y",r2)
# r2<-gsub("^([A-Z]+)-12$","\\1-1",r2)
# r2<-gsub("^([A-Z]+)$","\\1.x",r2,perl=T)
# r2<-gsub("^([A-Z]+)1$","\\1.y",r2,perl=T)
# r2<-gsub("^([A-Z]+)2$","\\1",r2,perl=T)
# r<-c(r1,r2)

# load("allData_RNA.Rdata")
######
# 总体样本
ct <- levels(as.factor(allMeta$celltype))
gs <- allSample$gene_symbol
dim <- list(c(ct),c(gs))
re <- matrix(NA,nrow = length(ct),ncol = nrow(allSample),dimnames = dim)
allSample$gene_symbol <- NULL
for (i in ct) {
  tem <- allSample[,rownames(allMeta[which(allMeta$celltype %in% i),])]
  tem <- apply(tem,2,as.numeric)
  re[i,] <- t(rowMeans(tem))
}

cor_ma  <- t(re)

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
save(r1,p1,r2,p2,file="cor_transcriptome.Rdata") 

png(file=paste("corplot_celltype_gene.png",sep = ""),width=850,height=570,res=96)
# c("complete", "ward", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid")
corrplot(r1, method = "color",
         order = "hclust", addrect = 6, # 当order为hclust时，可以为添加相关系数图添加矩形框
         hclust.method = "ward.D2",
         tl.cex=0.7, tl.col="black",  # 字体
         type = "full",  # 指定展示的方式
         diag = TRUE,  # 展示对角线结果
         cl.cex = 0.7,  # 图列
         col = colorRampPalette(colors = c("#0da9ce","white","#f49537"))(100),
         p.mat = p1, insig = "label_sig", sig.level = c(0.01), 
         pch.cex = .7, pch.col = "white")
dev.off()
pdf("corplot_celltype_gene.pdf",width=9,height=9)
corrplot(r1, method = "color",
         order = "hclust", addrect = 6, # 当order为hclust时，可以为添加相关系数图添加矩形框
         hclust.method = "ward.D2",
         tl.cex=0.7, tl.col="black",  # 字体
         type = "full",  # 指定展示的方式
         diag = TRUE,  # 展示对角线结果
         cl.cex = 0.7,  # 图列
         col = colorRampPalette(colors = c("#0da9ce","white","#f49537"))(100),
         p.mat = p1, insig = "label_sig", sig.level = c(0.01), 
         pch.cex = .7, pch.col = "white")
dev.off()
