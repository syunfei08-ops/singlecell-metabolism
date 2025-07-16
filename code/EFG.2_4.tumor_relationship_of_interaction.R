# 转录和代谢相关性与互作的关系
# 各肿瘤类型分别计算
library(Hmisc) # rcorr函数
library(corrplot)
library(RColorBrewer) # 配色包
library(paletteer) # 配色包
library(reshape2)
library(ggplot2)
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/3.metabolic_landscape/interaction/tumor")
source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")

# 代谢和转录数据读入
load("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/3.metabolic_landscape/transcriptome/allData_RNA.Rdata") 
gs <- allSample$gene_symbol
sam <- levels(as.factor(allMeta$sample))
tum <- unique(unlist(strsplit(sam,"-"))[seq(1,length(sam)*4,by=4)])

level <- "module" # "pathway","sub_pathway","module"
meta <- "high" # high, v1, v2, normal
scFEA_raw<-read_in_FC(level, meta)
ct <- sort(unique(scFEA_raw$celltype))
# 去除样本占比少的细胞类型
nSamCt <- c()
for (i in ct) {
  ct_sam <- scFEA_raw[which(scFEA_raw$celltype==i),]
  nSamCt <- c(nSamCt, length(unique(ct_sam$sample)))
}
names(nSamCt) <- ct
remct <- ct[which(as.numeric(nSamCt) < 5)]

scFEA <- ct_adj(scFEA_raw)
pw <- levels(as.factor(scFEA[,3]))
ct <- sort(unique(scFEA$celltype)) # 细胞类型统计

for(t in tum){
  print(t)
  # 转录本数据
  a <- which(unlist(strsplit(allMeta$sample,"-"))[seq(1,length(allMeta$sample)*4,by=4)] %in% t)
  metaTum <- allMeta[a,]
  matrTum <- allSample[,rownames(metaTum)]
  ct <- levels(as.factor(metaTum$celltype))
  dim <- list(c(ct),c(gs))
  re <- matrix(NA,nrow = length(ct),ncol = nrow(matrTum),dimnames = dim)
  for (i in ct) {
    if(length(which(metaTum$celltype %in% i)) > 0){
      tem <- matrTum[,which(metaTum$celltype %in% i)]
      tem <- apply(as.data.frame(tem),2,as.numeric)
      re[i,] <- t(rowMeans(tem))
    }
  }
  cor_ma  <- t(re)
  
  res_p <- rcorr(cor_ma, type = "pearson")
  r1 <- res_p$r
  r1[which(is.na(r1))] <- 1
  r1[which(is.infinite(r1))] <- 0
  # r1 <- round(r1,2)
  rR <- r1
  pR <- res_p$P
  pR[which(is.na(pR))] <- 0
  
  png(file=paste("corplot_",t,"_gene.png",sep = ""),width=850,height=570,res=96)
  # c("complete", "ward", "ward.D", "ward.D2", "single", "average", "mcquitty", "median", "centroid")
  corrplot(r1, method = "color",
           order = "hclust", addrect = 6, # 当order为hclust时，可以为添加相关系数图添加矩形框
           hclust.method = "median",
           tl.cex=0.7, tl.col="black",  # 字体
           type = "full",  # 指定展示的方式
           diag = TRUE,  # 展示对角线结果
           cl.cex = 0.7,  # 图列
           col = colorRampPalette(colors = c("#0da9ce","white","#f49537"),bias = 0.25)(100),
           p.mat = pR, insig = "label_sig", sig.level = c(0.01), 
           pch.cex = .7, pch.col = "white")
  dev.off()
  print(1)
  # 代谢数据
  a <- which(unlist(strsplit(scFEA$sample,"-"))[seq(1,length(scFEA$sample)*4,by=4)] %in% t)
  scFEATum <- scFEA[a,]
  ct <- levels(as.factor(scFEATum[,4]))
  dim <- list(c(ct),c(pw))
  re_ct <- matrix(nrow = length(ct),ncol = length(pw),dimnames = dim)
  for (i in ct) {
    for (j in pw) {
      re_ct[i,j] <- mean(scFEA[which(scFEA[,4]==i & scFEA[,3]==j),5])
    }
  }
  cor_ma  <- t(re_ct)
  
  res_p <- rcorr(cor_ma, type = "pearson")
  r1 <- res_p$r
  r1[which(is.na(r1))] <- 1
  r1[which(is.infinite(r1))] <- 0
  # r1 <- round(r1,2)
  rM <- r1
  pM <- res_p$P
  pM[which(is.na(pM))] <- 0
  
  png(file=paste("corplot_",t,"_module.png",sep = ""),width=850,height=570,res=96)
  corrplot(r1, method = "color",
           order = "hclust", addrect = 6, 
           hclust.method = "mcquitty",
           tl.cex=0.7, tl.col="black",
           type = "full",
           diag = TRUE,
           cl.cex = 0.7,
           col = colorRampPalette(colors = c("#0da9ce","white","#f49537"))(100),
           p.mat = pM, insig = "label_sig", sig.level = c(0.01), 
           pch.cex = .7, pch.col = "white")
  dev.off()
  print(2)
  # 互作数据
  tem_sam=levels(as.factor(scFEATum[,2]))
  n = length(tem_sam)  # 读取dir长度，也就是文件夹下的文件个数
  matrNum <- matrix(0,nrow = length(ct),ncol = length(ct),dimnames = list(c(ct),c(ct)))
  matrSum <- matrix(0,nrow = length(ct),ncol = length(ct),dimnames = list(c(ct),c(ct)))
  for (i in 1:n){
    dir = paste("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/cellphoneDB/",tem_sam[i],"_count_network.txt",sep="")
    intCount <- read.table(dir, header = T, sep = "\t")
    intCount <- intCount[which(!(intCount[,1] %in% c(remct,"Unassigned") | intCount[,2] %in% c(remct,"Unassigned"))),]
    for (j in 1:nrow(intCount)) {
      matrNum[intCount[j,1],intCount[j,2]] <- matrNum[intCount[j,1],intCount[j,2]]+1
      matrSum[intCount[j,1],intCount[j,2]] <- matrSum[intCount[j,1],intCount[j,2]]+intCount[j,3]
    }
  }
  matr <- as.data.frame(matrSum/matrNum)
  matr$ctX <- row.names(matr)
  inter <- melt(matr,id.vars=c("ctX"),variable.name="ctY",value.name="interaction")
  rM <- as.data.frame(rM)
  rM$ctX <- row.names(rM)
  rMet <- melt(rM,id.vars=c("ctX"),variable.name="ctY",value.name="cor")
  rR <- as.data.frame(rR)
  rR$ctX <- row.names(rR)
  rRNA <- melt(rR,id.vars=c("ctX"),variable.name="ctY",value.name="cor")
  
  matr <- data.frame("inter"=inter$interaction,"rMet"=rMet$cor,"rRNA"=rRNA$cor)
  # matr <- data.frame("inter"=inter$value,"rMet"=rMet$value,"rRNA"=rRNA$value)
  cor_ma <- as.matrix(matr)
  cor_ma <- apply(cor_ma,2,as.numeric)
  reP <- rcorr(cor_ma, type = "pearson")
  rP <- reP$r
  pP <- reP$P
  reS <- rcorr(cor_ma, type = "spearman")
  rS <- reS$r
  pS <- reS$P
  
  # dat <- as.data.frame(scale(matr))
  # dat$x <- rownames(dat)
  dat <- as.data.frame(matr)
  dat <- dat[which(!is.nan(dat$inter)),]
  
  col1<-"#2166AC"
  col2<-"#7BCD7B" # 4393C3,7BCD7B
  colS<-"#696773"
  p2<-ggplot(dat,aes(x=rRNA,y=inter))+
    geom_point(shape=20,colour=col2,size=3,stroke = 1)+
    geom_smooth(data=dat,method = "lm",colour=colS)+
    labs(x="Correlation of Transcriptome",y="Interaction")+
    theme_classic() + 
    theme(legend.position="none",
          axis.text.x=element_text(colour="black", size = 10,hjust=1,vjust=1),
          axis.text.y=element_text(colour="black", size = 10),
          axis.line=element_line(size=0.2,color="black"),
          axis.ticks = element_line(colour = "black",size=0.2),
          panel.border = element_blank(), panel.background = element_blank(),
          axis.ticks.length= unit(.5, "mm"))
  p <- p2+annotate("text", x = (max(dat$rRNA)+min(dat$rRNA))/2, y = max(dat$inter[which(!is.nan(dat$inter))]), 
                   label = paste0("(Pearson)R = ",round(rP["inter","rRNA"],2),",p = ",round(pP["inter","rRNA"],2)),
                   color="#000000",size = 3.5, angle=0, fontface="bold")+
    annotate("text", x = (max(dat$rRNA)+min(dat$rRNA))/2, y = max(dat$inter[which(!is.nan(dat$inter))])-7, 
             label = paste0("(Spearman)R = ",round(rS["inter","rRNA"],2),",p = ",round(pS["inter","rRNA"],2)),
             color="#000000",size = 3.5, angle=0, fontface="bold" )
  png(file=paste0("dotplot_",t,"_RNAInter.png"),width=850,height=570,res=96)
  print(p)
  dev.off()
  print(3)
  p1<-ggplot(dat,aes(x=rMet,y=inter))+
    geom_point(shape=20,colour=col1,size=3,stroke = 1)+
    geom_smooth(data=dat,method = "lm",colour=colS)+
    labs(x="Correlation of Metabolism",y="Interaction")+
    theme_classic() + 
    theme(legend.position="none",
          axis.text.x=element_text(colour="black", size = 10,hjust=1,vjust=1),
          axis.text.y=element_text(colour="black", size = 10),
          axis.line=element_line(size=0.2,color="black"),
          axis.ticks = element_line(colour = "black",size=0.2),
          panel.border = element_blank(), panel.background = element_blank(),
          axis.ticks.length= unit(.5, "mm"))
  p <- p1+annotate("text", x = mean(dat$rMet)/2, y = max(dat$inter[which(!is.nan(dat$inter))]), 
                   label = paste0("(Pearson)R = ",round(rP["inter","rMet"],2),",p = ",round(pP["inter","rMet"],2)),
                   color="#000000",size = 3.5, angle=0, fontface="bold")+
    annotate("text", x = mean(dat$rMet)/2, y = max(dat$inter[which(!is.nan(dat$inter))])-7, 
             label = paste0("(Spearman)R = ",round(rS["inter","rMet"],2),",p = ",round(pS["inter","rMet"],2)),
             color="#000000",size = 3.5, angle=0, fontface="bold")
  png(file=paste0("dotplot_",t,"_MetInter.png"),width=850,height=570,res=96)
  print(p)
  dev.off()
  print(4)
}

