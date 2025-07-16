# 肿瘤类型的代谢特征
source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/2-2.tumor/comparation")

# 构建样本信息
args=commandArgs(T)
g <- args[1]
cho_ct <- "Malignant cells" # "Macrophages"
meta_data <- build_meta(cho_ct) # 1:sample, 2:celltype, 3:tumor, 4:primary_site, 5:tissue, 6:group
group <- levels(as.factor(meta_data[,as.numeric(g)]))

library(dplyr)
library(ggplot2)
library(Hmisc)
library(paletteer)

level <- "pathway" # "pathway","module"
meta <- "high" # high, v1, v2, normal
scFEA_raw<-read_in_FC(level, meta)
scFEA <- ct_adj(scFEA_raw) # 小于5个样本的细胞类型去掉
scFEA <- scFEA[-which(is.infinite(scFEA$log2FC)),] # 去掉-Inf，影响下游计算

ct <- levels(as.factor(scFEA$celltype))
pw <- levels(as.factor(scFEA$pathway_1))
pw_cho <- c("hsa00010","hsa00020","hsa00230","hsa00240")
sam <- levels(as.factor(scFEA$sample))
scFEA$tumor <- sapply(strsplit(as.character(scFEA$sample), "-"), "[[", 1)
tum <- levels(as.factor(scFEA$tumor))
nSamTumor <- c()
for (i in tum){
  tem <- scFEA$sample[scFEA$tumor %in% i]
  nSamTumor <- c(nSamTumor, length(unique(tem)))
}
names(nSamTumor) <- tum

######
# 记录相关系数
library(dplyr)
library(ggplot2)
library(Hmisc)
library(paletteer)
library(corrplot)
# library(ggcorrplot)
library(ggthemes)

# for (t in tum) {
for (t in group) {
  corMatr <- matrix(NA,nrow = length(pw),ncol = length(pw),dimnames = list(pw,pw))
  pMatr <- matrix(NA,nrow = length(pw),ncol = length(pw),dimnames = list(pw,pw))
  for (cho1 in pw) {
    for (cho2 in pw) {
	  sCho <- meta_data[which(meta_data[,as.numeric(g)] %in% t),1]
      # tem <- scFEA[grep(t,scFEA$sample),]
	  tem <- scFEA[which(scFEA$sample %in% sCho),]
	  
      col1 <- tem[tem$pathway_1 %in% cho1,c(4,5)] %>% group_by(celltype) %>% summarize_each(funs(mean))
      col2 <- tem[tem$pathway_1 %in% cho2,c(4,5)] %>% group_by(celltype) %>% summarize_each(funs(mean))
      A <- merge(col1,col2,"celltype")
      A$label <- paste0(t,"-",A$celltype)
      A$tumor <- t

      colnames(A)<-c("celltype","flux1","flux2","label","tumor")
    
      # 将cho_ct的细胞类型标注实心
      cho_ct <- c("Malignant cells")
      A$group <- "other" 
      A$group[which(A$celltype %in% cho_ct)] <- A$tumor[which(A$celltype %in% cho_ct)]
    
      # 添加相关系数
      cor <- rcorr(A$flux1,A$flux2,type ="pearson")
      r<- cor$r[2]
      pi <- cor$P[2]
      corMatr[cho1,cho2] <- r
      pMatr[cho1,cho2] <- pi
    }
  }
  # cho <- c("hsa00010","hsa00020","hsa00230","hsa00240")
  # data <- corMatr[colnames(corMatr) %in% cho,]
  # pdf(paste0("corrplot_four.pdf"),width=4,height=10)
  data <- corMatr
  pdf(paste0(g,"-",t,"_corrplot_all.pdf"),width=10,height=10)
  corrplot(t(data), method = "color",outline = "white",
           hclust.method = "mcquitty",addCoef.col = "white",
           tl.col = "black",tl.cex = 0.7,
           number.digits = 2,number.cex = 0.6,number.font = NULL,
           cl.pos = "n")
  dev.off()
}

# ggcorrplot(corMatr,method = "circle",type = "upper", ggtheme = ggplot2::theme_bw(), title = "", 
#            show.legend = TRUE, legend.title = "Corr", show.diag = T, 
#            colors = c("#839EDB", "white", "#FF8D8D"), outline.color = "white", 
#            hc.order = T, hc.method = "complete", lab = T, lab_col = "black", 
#            lab_size = 2, p.mat = pMatr, sig.level = 0.01, insig = "blank", pch = 4, pch.col = "black", pch.cex = 5, tl.cex = 12, 
           # tl.col = "black", tl.srt = 45, digits = 2)

######
library(cowplot)

level <- "pathway" # "pathway","sub_pathway","module"
scFEA_raw<-read_in_flux(level,meta)
scFEA <- ct_adj(scFEA_raw)

ct <- levels(as.factor(scFEA$celltype))
pw <- levels(as.factor(scFEA$pathway_1))
pw_cho <- c("hsa00010","hsa00020","hsa00230","hsa00240")
sam <- levels(as.factor(scFEA$sample))
scFEA$tumor <- sapply(strsplit(as.character(scFEA$sample), "-"), "[[", 1)
tum <- levels(as.factor(scFEA$tumor))

for (t in group) {
  n=0
  sCho <- meta_data[which(meta_data[,as.numeric(g)] %in% t),1]
  t_scFEA <- scFEA[which(scFEA$sample %in% sCho),]
  for (cho1 in pw_cho) {
    for (cho2 in pw_cho) {
      n=n+1
      if(cho1==cho2)
        next
      A <- c()
      for (i in ct) {
        tem1 <- t_scFEA[which(t_scFEA$celltype %in% i & t_scFEA[,1] %in% cho1),]
        tem2 <- t_scFEA[which(t_scFEA$celltype %in% i & t_scFEA[,1] %in% cho2),]
        B1 <- c(i,cho1,summary(tem1$flux)[2:5])
        B2 <- c(i,cho2,summary(tem2$flux)[2:5])
        A <- rbind(A,B1,B2)
      }
      colnames(A)<-c("celltype","pathway","quantile_25","median","mean","quantile_75")
    
      # 画图
      dat <- as.data.frame(A)
      dat[,3:6] <- apply(dat[,3:6],2,as.numeric)
      p <- dat%>%
      mutate(pathway = factor(pathway,levels=c(cho1,cho2)))%>%
      ggplot(aes(x=celltype,color=pathway)) + 
        geom_errorbar(aes(ymin = quantile_25, ymax=quantile_75,group=pathway), #The error bar represents a 95% confidence interval
                      width=0.2, #The width of the short horizontal line at the end of the error bar
                      position=position_dodge(0.9), 
                      alpha = 1,
                      size=1) + 
        theme_bw()+
        #Plot the median as a dot plot
        geom_point(aes(x=celltype, y=median),pch=19,position=position_dodge(0.9),size=2) +
        scale_color_manual(values = c('#4393C3','#D6604D') )+ 
        #Fill color : c("#56B4E9", "#E69F00") c("#1F78B4","#CB181D") c('#4393C3','#D6604D') 
        theme(
          panel.grid.major = element_blank(),   #Grid lines are not displayed
          panel.grid.minor = element_blank(),   #Grid lines are not displayed
          axis.text.x=element_text(angle=300,hjust=0))+  
        xlab("Celltype")+ylab("log2FC")+ #title of X&Y,TRPM8 for figures3A, TRPA1 for figures3B, TRPM3 for figures3C
        geom_vline(xintercept=c(seq(1.5,(length(ct)-0.5),by=1)), linetype="dotted") #Adds a dotted line at the specified position
      assign(paste0("p",n),p)
      # png(file=paste0("errorbar_",cho1,cho2,"_allTumor.png"),width=1300,height=700,res=144)
      # print(p)
      # dev.off()
      # pdf(paste0("errorbar_",cho1,cho2,"_allTumor.pdf"),width=9,height=6)
      # print(p)
      # dev.off()
  	print(paste0(n,": ",cho1," and ",cho2))
    }
  }
  ######
  # 组合errorbar图
  p <- plot_grid(p2,p3,p4,p7,p8,p12, labels = c("A","B","C","D","E","F"), nrow = 3, align = "v")
  pdf(paste0(g,"-",t,"_groupSixErrorbar.pdf"),width=13,height=11)
  print(p)
  dev.off()
}
