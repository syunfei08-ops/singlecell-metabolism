######
# 检验样本批次
source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
args=commandArgs(T)
disease <- args[1] # "LUAD"

library(Seurat)
library(ggplot2)
library(paletteer) # 配色包

# sam <- read.csv("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/fluxandcelltype_high.list",header = F)
# sam$tumType <- sapply(strsplit(as.character(sam[,1]), "-"), "[[", 1)
# sam$dataset <- sapply(strsplit(as.character(sam[,1]), "-"), "[[", 2)
# sam$Sample <- sapply(strsplit(as.character(sam[,1]), ".", fixed = T), "[[", 1)

# tt <- unique(sam$tumType)

# s=sam[which(sam$tumType==disease),]
# scRNAlist <- list()
# for (i in 1:nrow(s)){
	# scRNAlist[[i]] <- readRDS(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/seurat/",s$Sample[i],"_seurat.rds"))
# }


# count.combined <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]], scRNAlist[[4]], scRNAlist[[5]], scRNAlist[[6]], scRNAlist[[7]], scRNAlist[[8]], scRNAlist[[9]], 
						# scRNAlist[[10]], scRNAlist[[11]], scRNAlist[[12]], scRNAlist[[13]], scRNAlist[[14]], scRNAlist[[15]], scRNAlist[[16]], scRNAlist[[17]], scRNAlist[[18]],
						# scRNAlist[[19]], scRNAlist[[20]], scRNAlist[[21]], scRNAlist[[22]], scRNAlist[[23]], scRNAlist[[24]], scRNAlist[[25]], scRNAlist[[26]], scRNAlist[[27]],
						# scRNAlist[[28]], scRNAlist[[29]], scRNAlist[[30]], scRNAlist[[31]], scRNAlist[[32]], scRNAlist[[33]], scRNAlist[[34]], scRNAlist[[35]], scRNAlist[[36]],
						# scRNAlist[[37]], scRNAlist[[38]], scRNAlist[[39]], scRNAlist[[40]], scRNAlist[[41]], scRNAlist[[42]], scRNAlist[[43]], scRNAlist[[44]], scRNAlist[[45]],
						# scRNAlist[[46]], scRNAlist[[47]], scRNAlist[[48]], scRNAlist[[49]], scRNAlist[[50]], scRNAlist[[51]], scRNAlist[[52]], scRNAlist[[53]], scRNAlist[[54]],
						# scRNAlist[[55]], scRNAlist[[56]], scRNAlist[[57]], scRNAlist[[58]], scRNAlist[[59]], scRNAlist[[60]], scRNAlist[[61]], scRNAlist[[62]], scRNAlist[[63]], 
						# scRNAlist[[64]]))

# save(count.combined,disease,file=paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/2.tumor_heterogeneity/",disease,".Rdata"))

load(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/2.tumor_heterogeneity/",disease,".Rdata"))
setwd(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/2.tumor_heterogeneity"))
# 全部细胞
Idents(count.combined) <- "celltype"
count.combined <- subset(x=count.combined,
                           idents = c("Unassigned"),invert=TRUE)
count.combined <- NormalizeData(count.combined, verbose = FALSE)
count.combined <- FindVariableFeatures(count.combined, selection.method = "vst", verbose = FALSE)
count.combined <- ScaleData(count.combined, verbose = FALSE)
count.combined <- ScaleData(count.combined, verbose = FALSE)
count.combined <- RunPCA(count.combined, npcs = 30, verbose = FALSE)
count.combined <- RunUMAP(count.combined, reduction = "pca", dims = 1:30)

ct.col <- color_ct_val() # 细胞类型的配色
ct.col <- ct.col[which(names(ct.col) %in% unique(count.combined$celltype))] # 特别注意大的集合在前面
# sam.col <- paletteer_d("khroma::soil")
sam.col <- c(paletteer_d("ggsci::category20_d3"),paletteer_d("ggsci::category20b_d3"),paletteer_d("ggsci::default_ucscgb")) # 样本的配色 66个
n_ct <- length(ct.col)
n_sam <- length(unique(count.combined$orig.ident))
  
p1 <- DimPlot(count.combined, reduction = "umap",
              cols = sam.col,
              group.by = "orig.ident") + 
  labs(title = "")
p2 <- DimPlot(count.combined, reduction = "umap", 
              cols = ct.col,
              group.by = "celltype", label = TRUE, repel = TRUE) + 
  labs(title = "")
if(n_sam > 20){
  png(filename = paste0(disease,"_sample.png"), 
      width = 1070, height = 749, res = 108)
}else{
  png(filename = paste0(disease,"_sample.png"), 
      width = 856, height = 749, res = 108)
}
plot(p1)
dev.off()
if(n_sam > 20){
  pdf(file = paste0(disease,"_sample.pdf"),
      width = 14, height = 12)
}else{
  pdf(file = paste0(disease,"_sample.pdf"),
      width = 12, height = 12)
}
plot(p1)
dev.off()
  
if(n_ct > 20){
  png(filename = paste0(disease,"_celltype.png"), 
      width = 1177, height = 749, res = 108)
}else{
  png(filename = paste0(disease,"_celltype.png"), 
      width = 963, height = 749, res = 108)
}
plot(p2)
dev.off()
if(n_ct > 20){
  pdf(file = paste0(disease,"_celltype.pdf"),
      width = 14, height = 12)
}else{
  pdf(file = paste0(disease,"_celltype.pdf"),
      width = 12, height = 12)
}
plot(p2)
dev.off()

# 恶性细胞
count.Malignant <- subset(x=count.combined,
                           idents = "Malignant cells")
count.Malignant <- NormalizeData(count.Malignant, verbose = FALSE)
count.Malignant <- FindVariableFeatures(count.Malignant, selection.method = "vst", verbose = FALSE)
count.Malignant <- ScaleData(count.Malignant, verbose = FALSE)
count.Malignant <- RunPCA(count.Malignant, npcs = 30, verbose = FALSE)
count.Malignant <- RunUMAP(count.Malignant, reduction = "pca", dims = 1:30)

sam.col <- c(paletteer_d("ggsci::category20_d3"),paletteer_d("ggsci::category20b_d3"),paletteer_d("ggsci::default_ucscgb")) # 样本的配色 66个
n_sam <- length(unique(count.Malignant$orig.ident))

p <- DimPlot(count.Malignant, reduction = "umap",
             cols = sam.col,
             group.by = "orig.ident") + 
  labs(title = "")
if(n_sam > 20){
  png(filename = paste0(disease,"_Malignant_sample.png"), 
      width = 1070, height = 749, res = 108)
}else{
  png(filename = paste0(disease,"_Malignant_sample.png"), 
      width = 856, height = 749, res = 108)
}
plot(p)
dev.off()
if(n_sam > 20){
  pdf(file = paste0(disease,"_Malignant_sample.pdf"),
      width = 14, height = 12)
}else{
  pdf(file = paste0(disease,"_Malignant_sample.pdf"),
      width = 12, height = 12)
}
plot(p)
dev.off()

