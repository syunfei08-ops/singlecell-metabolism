# /xtdisk/xiaojf_group/qianqh/tool/miniconda3/envs/rbase/bin/R

# 肿瘤类型的代谢描述
source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/2-3.normal/DESeq2")
options(digits) # 控制小数点

# 构建样本信息
meta_data <- normal_meta() 

level <- "pathway" # "pathway","module"
meta <- "normal" # high, v1, v2, normal
scFEA_raw<-read_in_flux(level, meta)
scFEA <- ct_normal(scFEA_raw,meta_data) # 小于5个样本的细胞类型去掉

ct <- levels(as.factor(scFEA$celltype))
type <- levels(as.factor(meta_data$class))
ds <- levels(as.factor(meta_data$dataset))
tum <- levels(as.factor(meta_data$tumor))
sam <- levels(as.factor(meta_data$sample))

t <- tum[4]
samT <- meta_data$sample[meta_data$tumor %in% t]
scFEAT <- scFEA[scFEA$sample %in% samT,]
metaT <- meta_data[meta_data$sample %in% samT,]

res <- c()
for (i in ct){
  # 构建代谢的样本bulk表达矩阵
  scFEA_ct <- scFEAT[scFEAT$celltype %in% i,]
  temSam <- levels(as.factor(scFEA_ct$sample))
  for (s in temSam){
    # s_scFEA <- scFEA_ct[scFEA_ct$sample %in% s,]
	A <- scFEA[scFEA$sample %in% s,]
	if(min(A$flux)<0){
      A$flux <- A$flux-min(A$flux)
    }
	s_scFEA <- A[A$celltype %in% i,]
    # a <- aggregate(log2FC ~ module_1, data=scFEA_ct, FUN=mean)
    a <- s_scFEA[,"flux",drop=F]
    rownames(a) <- s_scFEA[,1]
    colnames(a) <- c(s)
    if(s==temSam[1]){
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
  sam1 <- intersect(colnames(expr),metaT$sample[which(metaT$class=="Normal")])
  sam2 <- intersect(colnames(expr),metaT$sample[which(metaT$class=="Tumour")])
  expr1 <- expr[,sam1,drop=F]
  if(any(is.na(expr1))){
    expr1 <- expr1[,-which(is.na(colSums(expr1))),drop=F]
  }
  expr2 <- expr[,sam2,drop=F]
  if(any(is.na(expr2))){
    expr2 <- expr2[,-which(is.na(colSums(expr2))),drop=F]
  }
  expr1 <- cbind(expr1,expr1) # 正常样本
  expr2 <- cbind(expr2,expr2) # 肿瘤样本
  re <- log2(rowMeans(expr2)/rowMeans(expr1)) # 肿瘤相对于正常
  # print(paste0(t,": ", i))
  # print(paste0("re > 1"))
  # print(sort(re[which(re>1)],decreasing=T))
  # print(paste0("re < -1"))
  # print(sort(re[which(re< -1)],decreasing=T))
  if(i==ct[1]){
    res <- re
  }else{
    res <- cbind(res,re)
  }
}
colnames(res) <- ct
rownames(res) <- rN
write.table(res, file=paste0("normal-",level,"_",t,"_differential.txt"),sep="\t")

# 画图
library(pheatmap)
library(paletteer)
library(corrplot)
library(RColorBrewer)

t <- "BCC" # "BCC", "BRCA", "ccRCC", "PDAC"
level <- "module" # "pathway","module"
# res <- read.table(paste0("normal-",level,"_",t,"_differential.txt"),sep="\t", check.names=FALSE)
res <- read.table(paste0("normal-",level,"_differential.txt"),sep="\t", check.names=FALSE)
res <- res[ , colSums(is.na(res)) == 0]
res[sapply(res, is.infinite)] <- NA

# 数量统计的函数
f<-function(x) sum(abs(x)>1)
# apply(res,2,f)
# apply(res,1,f)

annotation_col = data.frame(
  nDiffM = apply(res,1,f),
  MetChangeM = rowSums(res)
)
# rownames(annotation_col) = paste("Test", 1:10, sep = "")

annotation_row = data.frame(
  nDiffCT = apply(res,2,f),
  MetChangeCT = colSums(res)
)
# rownames(annotation_row) = paste("Gene", 1:20, sep = "")

# 自定注释信息的颜色列表
ann_colors = list(
  nDiffCT = colorRampPalette(colors = c("#008277","white","#9F7200"))(100),
  nDiffM = colorRampPalette(colors = c("#008277","white","#9F7200"))(100),
  MetChangeM = colorRampPalette(colors = c("white","#00A1A4"))(100),
  MetChangeCT = colorRampPalette(colors = c("white","#5F54A2"))(100)
)

p <- pheatmap(t(res),
              color = colorRampPalette(colors = c("#005691","white","#CF191A"))(100),
              # color = colorRampPalette(colors = c(rep("royalblue",3),"white",rep("firebrick3",3))(100),
              border_color = "white",
              # scale = "column",
              fontsize_row = 12,fontsize_col = 12,
              # cellheight = 10,cellwidth = 7,
              cluster_row = T,cluster_cols = T,
              clustering_distance_rows = "euclidean",
              clustering_method="single",
              annotation_row = annotation_row, annotation_col = annotation_col,
              annotation_legend = T, annotation_colors = ann_colors,
			  na_col = "#DDDDDD"
)
if(level=="module"){
  wid=36
}else{
  wid=24
}

# pdf(paste0("heatmap-normal-",level,"_",t,"_differential.pdf"),width=wid,height=9)
pdf(paste0("heatmap-normal-",level,"_differential.pdf"),width=wid,height=9)
print(p)
dev.off()
