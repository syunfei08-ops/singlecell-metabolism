# 各细胞类型的总体代谢能力
source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/3.metabolic_landscape/overall")

level <- "module" # "pathway","sub_pathway","module"
meta <- "high" # high, v1, v2, normal
scFEA_raw<-read_in_FC(level, meta)
scFEA <- ct_adj(scFEA_raw) # 小于5个样本的细胞类型去掉

library(ggplot2)
ct <- levels(as.factor(scFEA$celltype))
re <- c()
for (i in ct) {
  scFEA_tem <- scFEA[scFEA$celltype %in% i,] 
  tem <- aggregate(log2FC~sample,scFEA_tem,mean)
  tem$celltype <- i
  re <- rbind(re,tem)
}
# mean(scFEA[scFEA$celltype %in% "Malignant cells",]$log2FC)
# mean_Malignant=0.06510131, mean_Macrophage=0.0827283, mean_DC=0.06001711
m<-c()
for (i in ct) {
  t <- mean(scFEA[scFEA$celltype %in% i,]$log2FC)
  m <- c(m,t)
}
write.csv(m,file="overall_celltype_mean.csv")

# 细胞类型排序
# til <- ct[c(1,18,20,25,2:9,28,23)]
lev <- c("Activated B cells","Memory B cells","Naive B cells","Plasma cells","CD4+ effector T cells","CD4+ naive T cells","CD4+ Tcm","CD4+ Tem","CD8+ effector T cells","CD8+ naive T cells","CD8+ Tcm","CD8+ Tem","Tregs","NK cells","Neutrophils","Mast cells","Plasmacytoid dendritic cells","Dendritic cells","Monocytes","Macrophages","Malignant cells","Fibroblasts","Endothelial cells","Epithelial cells","Smooth muscle cells","Neurons","Oligodendrocytes","Erythrocytes")

# 画图
sc <- as.data.frame(re)
sc$celltype <- factor(sc$celltype, levels = lev)
color <- color_ct_val()
# yin <- mean(scFEA$log2FC)
p<-ggplot(sc,aes(x=celltype,y=log2FC,fill=celltype)) +
  scale_y_continuous(limits=c(-1,1),breaks=-1:1,labels=-1:1)+
  geom_violin(trim=F,size=0.2,show.legend = F,width=1.0) + labs(y=NULL,x=NULL) + 
  stat_summary(fun.y = mean,geom="point",size=1,color="#1B1919")+
  geom_hline(aes(yintercept=0),color="#999d9c",linetype="dashed",size=0.2)+
  scale_fill_manual(name="Cell type", values=color)+
  scale_color_manual(name="Cell type", values=color)+
  labs(y="Log2(FoldChange)")+
  theme_classic() + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 9,angle=60,hjust=1,vjust=1),
        axis.text.y=element_text(colour="black", size = 9),
        axis.line=element_line(size=0.2,color="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"))

pdf("Overall_metabolic_landscape.pdf",width=11,height=6)
print(p)
dev.off()
jpeg(filename="Overall_metabolic_landscape.png",width=1000,height=600,res=120)
print(p)
dev.off()


######
# 癌种的各细胞类型的比较
library(pheatmap)
library(paletteer)
sam <- levels(as.factor(scFEA$sample))
tum <- unique(unlist(strsplit(sam,"-"))[seq(1,length(sam)*4,by=4)])
ct <- lev
dat <- matrix(NA,nrow = length(tum),ncol = length(ct),dimnames = list(tum,ct))
for (j in ct) {
  tem <- scFEA[scFEA$celltype %in% j,]
  for (i in tum) {
    dat[i,j] <- mean(tem$log2FC[grepl(i,tem$sample)])
  }
}
write.csv(dat,file="overall_tumor_mean.csv")

p <- pheatmap(t(dat),
              # color = paletteer_c("grDevices::Oslo", 30),
              # color = rev(paletteer_c("ggthemes::Red-Blue-White Diverging", 30)),
              color = rev(paletteer_c("ggthemes::Red-Blue-White Diverging", 30)),
              border_color = "white",
              # scale = "column",
              fontsize_row = 10,fontsize_col = 8,
              # cellheight = 10,cellwidth = 7,
              cluster_row = F,cluster_cols = F,
              # clustering_distance_rows = "euclidean",
              # clustering_method="single",
      )
pdf("overall_tumor.pdf",width=6,height=11)
print(p)
dev.off()
jpeg(filename="overall_tumor.png",width=1000,height=600,res=120)
print(p)
dev.off()

# which(dat[,"Macrophages"]>dat[,"Malignant cells"])
# "BCC","BRCA","ccRCC","CRC","ESCC","HNSCC","LUAD","PRAD"

######
# 随机扰动，取1/3样本，看总体趋势是否和所有样本一样
# color <- color_ct_val()
# for (i in 1:100) {
  # tem <- sort(sample(unique(scFEA$sample),100))
  # reTem <- re[re$sample %in% tem,]
  # sc <- as.data.frame(reTem)
  # sc$celltype <- factor(sc$celltype, levels = lev)
  # p<-ggplot(sc,aes(x=celltype,y=log2FC,fill=celltype)) +
    # scale_y_continuous(limits=c(-1,1),breaks=-1:1,labels=-1:1)+
    # geom_violin(trim=F,size=0.2,show.legend = F,width=1.0) + labs(y=NULL,x=NULL) + 
    # stat_summary(fun.y = mean,geom="point",size=1,color="#1B1919")+
    # geom_hline(aes(yintercept=0),color="#999d9c",linetype="dashed",size=0.2)+
    # scale_fill_manual(name="Cell type", values=color)+
    # scale_color_manual(name="Cell type", values=color)+
    # labs(y="Log2(FoldChange)")+
    # theme_classic() + 
    # theme(legend.position="none",
          # axis.text.x=element_text(colour="black", size = 9,angle=60,hjust=1,vjust=1),
          # axis.text.y=element_text(colour="black", size = 9),
          # axis.line=element_line(size=0.2,color="black"),
          # axis.ticks = element_line(colour = "black",size=0.2),
          # panel.border = element_blank(), panel.background = element_blank(),
          # axis.ticks.length= unit(.5, "mm"))
  # jpeg(filename=paste0("random",i,"_overall_landscape.png"),width=1000,height=600,res=120)
  # print(p)
  # dev.off()
# }
# 随机扰动，取1/3样本，利用平均值表示每次的代谢结果，做100次的分布
color <- color_ct_val()
re100 <- c()
for (i in 1:100) {
  tem <- sort(sample(unique(scFEA$sample),100))
  reTem <- re[re$sample %in% tem,]
  sc <- as.data.frame(reTem)
  sc$celltype <- factor(sc$celltype, levels = lev)
  tem <- aggregate(log2FC~celltype,sc,mean)
  tem$a <- i
  if(i==1){
    re100 <- tem
  }else{
    re100 <- rbind(re100,tem)
  }
}  
p<-ggplot(re100,aes(x=celltype,y=log2FC,fill=celltype)) +
  # scale_y_continuous(limits=c(-1,1),breaks=-1:1,labels=-1:1)+
  geom_violin(trim=F,size=0.2,show.legend = F,width=1.0) + labs(y=NULL,x=NULL) + 
  stat_summary(fun.y = mean,geom="point",size=1,color="#1B1919")+
  geom_hline(aes(yintercept=0),color="#999d9c",linetype="dashed",size=0.2)+
  scale_fill_manual(name="Cell type", values=color)+
  scale_color_manual(name="Cell type", values=color)+
  labs(y="Log2(FoldChange)")+
  coord_cartesian(ylim = c(-0.3,0.1))+
  theme_classic() + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 9,angle=60,hjust=1,vjust=1),
        axis.text.y=element_text(colour="black", size = 9),
        axis.line=element_line(size=0.2,color="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"))
pdf("100sample_overall.pdf",width=8,height=6)
print(p)
dev.off()
jpeg(filename="100sample_overall.png",width=1000,height=600,res=120)
print(p)
dev.off()
