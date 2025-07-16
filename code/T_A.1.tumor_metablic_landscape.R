# 各细胞类型的总体代谢能力
source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/2-2.tumor/landscape")

g <- 7
cho_ct <- "Macrophages" # "Macrophages", "Malignant cells"
meta_data <- build_meta(cho_ct) # 1:sample, 2:celltype, 3:tumor, 4:primary_site, 5:tissue, 6:group1, 7:group2
group <- levels(as.factor(meta_data[,as.numeric(g)]))

level <- "module" # "pathway","sub_pathway","module"
meta <- "high" # high, v1, v2, normal
scFEA_raw<-read_in_FC(level, meta)
scFEA <- ct_adj(scFEA_raw) # 小于5个样本的细胞类型去掉

library(ggplot2)

scFEA$tumor <- unlist(strsplit(scFEA$sample,"-"))[seq(1,length(scFEA$sample)*4,by=4)]
tum <- levels(as.factor(scFEA$tumor))

# 细胞类型排序
# til <- ct[c(1,18,20,25,2:9,28,23)]
lev <- c("Activated B cells","Memory B cells","Naive B cells","Plasma cells","CD4+ effector T cells","CD4+ naive T cells","CD4+ Tcm","CD4+ Tem","CD8+ effector T cells","CD8+ naive T cells","CD8+ Tcm","CD8+ Tem","Tregs","NK cells","Neutrophils","Mast cells","Plasmacytoid dendritic cells","Dendritic cells","Monocytes","Macrophages","Malignant cells","Fibroblasts","Endothelial cells","Epithelial cells","Smooth muscle cells","Neurons","Oligodendrocytes","Erythrocytes")



for (t in tum) {
	t_scFEA <- scFEA[scFEA$tumor %in% t,]
	ct <- levels(as.factor(t_scFEA$celltype))
	re <- c()
	for (i in ct) {
	  scFEA_tem <- t_scFEA[t_scFEA$celltype %in% i,] 
	  tem <- aggregate(log2FC~sample,scFEA_tem,mean)
	  tem$celltype <- i
	  re <- rbind(re,tem)
	}
	# mean(scFEA[scFEA$celltype %in% "Malignant cells",]$log2FC)
	# mean_Malignant=0.06510131, mean_Macrophage=0.0827283, mean_DC=0.06001711
	m<-c()
    for (i in ct) {
      mt <- mean(t_scFEA[t_scFEA$celltype %in% i,]$log2FC)
      m <- c(m,mt)
    }
    write.csv(m,file=paste0(t,"_celltype_mean.csv"))

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

	pdf(paste0(t,"_metabolic_landscape.pdf"),width=11,height=6)
	print(p)
	dev.off()
	jpeg(filename=paste0(t,"_metabolic_landscape.png"),width=1000,height=600,res=120)
	print(p)
	dev.off()

}
