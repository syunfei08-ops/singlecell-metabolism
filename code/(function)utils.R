# 保存函数
# FC数据读取
read_in_FC <- function(lev,meta){
  if(lev=="pathway"){
    dat <- read.csv(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/input_df/samples.34pathways_",meta,".log2FC"),
             sep = "\t", quote = "")
  }
  if(lev=="module"){
    dat <- read.csv(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/input_df/samples.168modules_",meta,".log2FC"),
             sep = "\t", quote = "")
  }
  return(dat)
}

# 标准化后的flux值读取 scFEA+min
read_in_flux <- function(lev,meta){
  if(lev=="pathway"){
    dat <- read.csv(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/input_df/samples.34pathways_",meta,".flux"),
                    sep = "\t", quote = "")
  }
  if(lev=="module"){
    dat <- read.csv(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/input_df/samples.168modules_",meta,".flux"),
                    sep = "\t", quote = "")
  }
  return(dat)
}

metabolics_gene <- function(se){
  # 获取代谢基因
  library(fgsea)
  pathway_file <- "/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/KEGG_metabolism.gmt"
  pathways <- gmtPathways(pathway_file)
  metabolics <- unique(as.vector(unname(unlist(pathways))))
  
  # 只保留代谢基因
  # counts <- GetAssayData(object = se, slot = "counts")
  # filtered_counts <- counts[rownames(se) %in% metabolics, ]
  # se <- CreateSeuratObject(filtered_counts, meta.data = se@meta.data)
  se <- se[rownames(se) %in% metabolics, ]
  
  return(se)
}

# 定义细胞类型的颜色
color_ct_val <- function(){
  values=c("Malignant cells" = "#4F8AC1", "CD4+ naive T cells" = "#00FFFF", "CD4+ effector T cells" = "#87CEFA", "CD4+ Tcm" = "#000080", "CD4+ Tem" = "#4169E1", "Tregs" = "#7FFF00", "CD8+ naive T cells" = "#DAA520", "CD8+ effector T cells" = "#7BCD7B", "CD8+ Tcm" = "#218A21", "CD8+ Tem" = "#2F4F4F", "NK cells" = "#8B0000", "Naive B cells" = "#87481F", "Activated B cells" = "#FFB6C1", "Memory B cells" = "#7F7F7F", "Plasma cells" = "#FFD700", "Macrophages" = "#FAB17A", "Monocytes" = "#FF6246", "Dendritic cells" = "#FF1493", "Mast cells" = "#FF0000", "Neutrophils" = "#FFDEAD", "Endothelial cells" = "#0000FF", "Epithelial cells" = "#708090", "Fibroblasts" = "#440053", "Erythrocytes" = "#483D8B", "Oligodendrocytes" = "#D2691E", "Microglial cells" = "#000000", "Neural progenitor cells" = "#D8BFD8", "Neurons" = "#D2B48C", "Astrocytes" = "#7E5094", "Acinar cells" = "#800080", "Progenitor cells" = "#9370DB", "HSCs" = "#FF00FF", "GMPs" = "#E6E6FA", "Unassigned" = "#BEBEBE", "Megakaryocytes" = "#BC8F8F", "Basophils" = "#6A5ACD", "Basal cells" = "#F6D7BE", "Smooth muscle cells" = "#F660C6", "Melanocytes" = "#FFFACD", "Myelocytes" = "#F0FFFF", "Neuroblasts" = "#C1CDCD", "Oligodendrocyte precursor cells" = "#8B668B", "Pericytes" = "#CDB5CD", "Plasmacytoid dendritic cells" = "#008B8B", "CLPs" = "#96CDCD", "CMPs" = "#FAF0E6", "Chondrocytes" = "#53868B", "Endocrine cells" = "#00F5FF", "Enterocytes" = "#FFDAB9", "Erythroid progenitor cells" = "#8B3626", "Hepatocytes" = "#8B8378", "Keratinocytes" = "#E0EEE0", "MEPs" = "#8B1C62", "Mesenchymal cells" = "#836FFF", "Muscle cells" = "#27408B", "Pancreatic ductal cells" = "#6CA6CD", "Platelets" = "#CAE1FF", "Retinal ganglion cells" = "#B4CDCD", "Astrocyte precursor cells" = "#8B7B8B", "Cajal-retzius cells" = "#1C1C1C", "Ciliated cells" = "#FFE1FF", "Club cells" = "#5D478B", "Eosinophils" = "#4A708B", "Glial cells" = "#B9D3EE", "Intermediated cells" = "#4F4F4F", "Interneurons" = "#6E7B8B", "Leukocytes" = "#7A8B8B", "Loop of Henle cells" = "#838B8B", "Mesothelial cells" = "#FFE4E1", "Mucosa cells" = "#8B8386", "Pinealocytes" = "#CD8C95", "Pneumocytes" = "#8B3A62", "Purkinje cells" = "#EEE9BF", "Pyramidal cells" = "#CDC5BF", "Renal principal cells" = "#FFA500", "Stem cells" = "#CD8500", "Adipocytes" = "#778899", "Plasmablasts" = "#ADD8E6", "MPPs" = "#B0C4DE", "Thyroid follicular cells" = "#9ACD32")
  return(values)
}

color_tumor_val <- function(){
  values=c("ATC"="#1F78B4FF","BCC"="#A6CEE3FF","BRCA"="#2FE4C0","ccRCC"="#BD26E3","CRC"="#CAB2D6FF",
         "DCIS"="#738BED","ESCC"="#F45B5BFF","HNSCC"="#FB9A99FF","LUAD"="#FF7F00FF",
         "LUSC"="#19BE4A","MCC"="#B2DF8AFF","MIUBC"="#18BC9CFF","NB"="#8CBE19","NBL"="#BE9D19",
         "NSCLC"="#CCBE93FF","PDAC"="#BE4A19","PRAD"="#F4E5A9","pRCC"="#E6C43D","STAD"="#93A2CC",
         "TNBC"="#93BECC","UCEC"="#9C8749","other"="#FFFFFF")
  return(values)
}

# 去除细胞类型（小于5个样本的细胞类型去掉）
ct_adj <- function(A){
  ct <- sort(unique(A$celltype)) # 细胞类型统计
  # 细胞类型的样本数目
  nSamCt <- c()
  for (i in ct) {
    ct_sam <- A[which(A$celltype==i),]
    nSamCt <- c(nSamCt, length(unique(ct_sam$sample)))
  }
  names(nSamCt) <- ct
  re <- ct[which(as.numeric(nSamCt) < 5)]
  B<-A[which(!A$celltype %in% c(re,"Unassigned")),]
  return(B)
}


# fpkm转化tpm
fpkmToTpm <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

# meta信息生成
build_meta <- function(cho_ct){
  # 样本根据malignant cells和macrophages的代谢能力分组
  groupMC <- c(1,1,6,4,1,2,6,24,1,6,3)
  names(groupMC) <- c("ATC", "DCIS", "LUSC", "MCC", "MIUBC", "NBL", "NSCLC", "PDAC", "pRCC", "TNBC", "UCEC")
  groupM <- c(4,21,19,16,7,3)
  names(groupM) <- c("BCC", "BRCA", "ccRCC", "CRC", "HNSCC", "NB")
  groupMid <- c(56,64,32)
  names(groupMid) <- c("ESCC", "LUAD", "PRAD")
  groupNO <- c(1)  # 不包含macrophages的STAD样本(1)被排除
  names(groupNO) <- c("STAD")
  
  info <- read.csv("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/data_meta_new.csv",header=F)
  colnames(info) <- info[1,]
  info <- info[-1,]
  meta <- data.frame("sample"=info[,c(3)])
  meta$celltype <- cho_ct
  meta$tumor <- unlist(strsplit(meta$sample,"-"))[seq(1,length(meta$sample)*4,by=4)]
  meta$primary_site <- info[,c(17)]
  meta$tissue <- info[,c(18)]
  
  meta[which(meta$tumor %in% names(groupMC)),"group1"] <- "High_Malignant cells_tumor"
  meta[which(meta$tumor %in% names(groupM)),"group1"] <- "High_Macrophages_tumor"
  meta[which(meta$tumor %in% names(groupMid)),"group1"] <- "Middle_tumor"
  meta[which(meta$tumor %in% names(groupNO)),"group1"] <- "no_Macrophages_tumor"
  
  level <- "module" # "pathway","module"
  m <- "high" # high, v1, v2, normal
  scFEA <-read_in_FC(level, m)
  scFEA_MC <- scFEA[scFEA$celltype %in% "Malignant cells",]
  scFEA_M <- scFEA[scFEA$celltype %in% "Macrophages",]
  temMC <- aggregate(log2FC~sample,scFEA_MC,mean)
  temM <- aggregate(log2FC~sample,scFEA_M,mean)
  tem <- merge(temMC,temM,by="sample")
  groupMC <- tem$sample[which(tem[,2]>tem[,3])]
  groupM <- tem$sample[which(tem[,2]<tem[,3])]
  meta$group2 <- "no_Macrophages_sample"
  meta[which(meta$sample %in% groupM),"group2"] <- "High_Macrophages_sample"
  meta[which(meta$sample %in% groupMC),"group2"] <- "High_Malignant cells_sample"
  
  return(meta)
}

# 正常样本的meta信息
normal_meta <- function(){
  info <- read.csv("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/data_meta_normal_new.csv",header=F)
  colnames(info) <- info[1,]
  info <- info[-1,]
  meta <- data.frame("sample"=info[,c(3)])
  meta$dataset <- info[,c(2)]
  meta$class <- info[,c(4)]
  meta$tumor <- info[,c(16)]
  meta$primary_site <- info[,c(17)]
  meta$tissue <- info[,c(18)]
  return(meta)
}

# 正常样本的细胞类型
ct_normal <- function(scFEA_raw,meta_data){
  samN <- meta_data$sample[meta_data$class %in% "Normal"]
  scFEAN <- scFEA_raw[scFEA_raw$sample %in% samN,]
  scFEAN <- ct_adj(scFEAN) # 小于5个样本的细胞类型去掉
  ctN <- levels(as.factor(scFEAN$celltype))

  samT <- meta_data$sample[meta_data$class %in% "Tumour"]
  scFEAT <- scFEA_raw[scFEA_raw$sample %in% samT,]
  scFEAT <- ct_adj(scFEAT) # 小于5个样本的细胞类型去掉
  ctT <- levels(as.factor(scFEAT$celltype))

  ct <- intersect(ctN,ctT)
  scFEA<-scFEA_raw[which(scFEA_raw$celltype %in% c(ct)),]
  
  return(scFEA)
}