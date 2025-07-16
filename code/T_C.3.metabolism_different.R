# /xtdisk/xiaojf_group/qianqh/tool/miniconda3/envs/rbase/bin/R

# 肿瘤类型的代谢描述
source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/2-2.tumor/DESeq2")

level <- "module" # "pathway","module"
meta <- "high" # high, v1, v2, normal
scFEA_raw<-read_in_flux(level, meta)
scFEA <- ct_adj(scFEA_raw) # 小于5个样本的细胞类型去掉

scFEA$tumor <- unlist(strsplit(scFEA$sample,"-"))[seq(1,length(scFEA$sample)*4,by=4)]
sam <- levels(as.factor(scFEA$sample))
tum <- levels(as.factor(scFEA$tumor))
pw <- levels(as.factor(scFEA[,3]))

cho_ct <- "Macrophages" # "Macrophages" "Malignant cells"
scFEA_ct <- scFEA[scFEA$celltype %in% cho_ct,]

# 构建代谢的样本bulk表达矩阵
for (s in sam){
  s_scFEA <- scFEA_ct[scFEA_ct$sample %in% s,]
  # a <- aggregate(log2FC ~ module_1, data=scFEA_ct, FUN=mean)
  a <- s_scFEA[,"flux",drop=F]
  rownames(a) <- s_scFEA[,1]
  colnames(a) <- c(s)
  if(s==sam[1]){
    expr <- a
  }else{
    expr <- merge(expr,a,by="row.names",all=T)
	rownames(expr) <- expr[,1]
	expr <- expr[,-1]
  }
}
rN <- rownames(expr)
expr <- apply(as.matrix(expr),2,as.numeric)
rownames(expr) <- rN

# 构建样本信息
args=commandArgs(T)
g <- args[1]
meta_data <- build_meta(cho_ct) # 1:sample, 2:celltype, 3:tumor, 4:primary_site, 5:tissue, 6:group1, 7:group2
group <- levels(as.factor(meta_data[,as.numeric(g)]))

sam1 <- meta_data[which(meta_data[,as.numeric(g)]==group[1]),"sample",drop=T]
sam2 <- meta_data[which(meta_data[,as.numeric(g)]==group[2]),"sample",drop=T]

# 计算差异表达基因的方法DESeq2等行不通，不能有负数，不能有非整数
# library(DESeq2)
# # 构建dds对象
# ## 制作dds对象，构建差异基因分析所需的数据格式
# dds <- DESeqDataSetFromMatrix(countData = expr, colData = meta, design = ~ group2)

# # countData = B，readscount矩阵
# # colData = coldata,分组信息，根据这个才能在两组之间比较
# # design = ~ condition，公式，表示按照condition进行分析

# ## 5.差异分析结果
# dds <- DESeq(dds)	#正式进行差异分析

# 用差异倍数(fold change)来计算
expr1 <- expr[,sam1]
if(any(is.na(expr1))){
expr1 <- expr1[,-which(is.na(colSums(expr1)))]
}
expr2 <- expr[,sam2]
if(any(is.na(expr2))){
expr2 <- expr2[,-which(is.na(colSums(expr2)))]
}

re <- log2(rowMeans(expr1)/rowMeans(expr2))

print(re[which(re>1)])
print(re[which(re< -1)])
