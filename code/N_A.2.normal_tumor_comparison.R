source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")

library(ggplot2)
library(ggunchained) # geom_split_violin
library(ggpubr) # stat_compare_means
library(gghalves)
library(ggsci)
library(scales)

level <- "module" # "pathway","sub_pathway","module"
meta <- "normal" # high, v1, v2, normal
meta_data <- normal_meta() 
scFEA_raw<-read_in_FC(level, meta)
scFEA <- ct_normal(scFEA_raw,meta_data) # 小于5个样本的细胞类型去掉
sample <- read.table("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/fluxandcelltype_normal.list")
sample$sam <- sapply(strsplit(as.character(sample$V1), ".", fixed = T), "[[", 1)
sample$type <- "tumor"
sample[1:33,"type"] <- "normal"
nSam <- sample[which(sample$type=="normal"),"sam"]
scFEAN <- scFEA[scFEA$sample %in% nSam,]

level <- "module" # "pathway","sub_pathway","module"
meta <- "high" # high, v1, v2, normal
scFEAT <- read_in_FC(level, meta)

ct <- intersect(unique(scFEAN$celltype),unique(scFEAT$celltype))
lev <- c("Activated B cells","Memory B cells","Naive B cells","Plasma cells","CD4+ effector T cells","CD4+ naive T cells","CD4+ Tcm","CD4+ Tem","CD8+ effector T cells","CD8+ naive T cells","CD8+ Tcm","CD8+ Tem","Tregs","NK cells","Neutrophils","Mast cells","Plasmacytoid dendritic cells","Dendritic cells","Monocytes","Macrophages","Malignant cells","Fibroblasts","Endothelial cells","Epithelial cells","Smooth muscle cells","Neurons","Oligodendrocytes","Erythrocytes")

reN <- c()
m <- c()
for (i in ct) {
  scFEA_tem <- scFEAN[scFEAN$celltype %in% i,] 
  tem <- aggregate(log2FC~sample,scFEA_tem,mean)
  tem$celltype <- i
  reN <- rbind(reN,tem)
}
reN$type <- "normal"

reT <- c()
m <- c()
for (i in ct) {
  scFEA_tem <- scFEAT[scFEAT$celltype %in% i,] 
  tem <- aggregate(log2FC~sample,scFEA_tem,mean)
  tem$celltype <- i
  reT <- rbind(reT,tem)
}
reT$type <- "tumor"

# scFEAN$type <- "normal"
# scFEAT$type <- "tumor"

dat <- rbind(reT,reN)
# dat <- rbind(scFEAN,scFEAT)
dat$celltype <- factor(dat$celltype, levels = lev)

p <- ggplot(dat,aes(x=celltype,y=log2FC,fill=type))+
  geom_split_violin(trim = T,colour="white", scale = 'width')+
  geom_point(stat = 'summary',fun=median,size=0.5,position = position_dodge(width = 0.5))+
  scale_fill_manual(values = c("#5773CC","#FFB900"))+
  stat_summary(fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = 'errorbar',
               color='black',width=0.18,size=0.2,
               position=position_dodge(width = 0.5))+
  ylab("scFEA predicted flux")+xlab("")+
  scale_y_continuous(limits=c(-0.5,0.5),breaks=-0.5:0.5,labels=-0.5:0.5)+
  theme_classic()+
  theme(panel.grid.major=element_line(color='gray',linetype='dotted',size=0.1),
        axis.text.x=element_text(angle=300,hjust=0))
  # stat_compare_means(label="p.signif",method="wilcox.test",size = 3,label.y=10)
  # stat_compare_means(data = dat, aes(x = celltype,y = log2FC),
  #                    symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
  #                                     symbols = c("***", "**", "*", "-")),
  #                    label = "p.signif",
  #                    label.y = 0.5,
  #                    hide.ns = F)

pdf("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/2-3.normal/landscape_split.pdf",width=30,height=12)
print(p)
dev.off()