source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/4.module_pathway/statistic")

######
# 柱状图，显示已选通路的样本数
library(ggplot2)
library(magrittr) # %>%
library(dplyr) # mutate
library(ggbreak)

level <- "pathway" # "pathway","module"
meta <- "high" # high, v1, v2, normal

sumU <- read.table(paste0(level,"_mean_up_",meta,".txt"),sep="\t")
sumD <- read.table(paste0(level,"_mean_down_",meta,".txt"),sep="\t") # log2FC < 0 in more than 50% of samples

# "Dendritic cells","Fibroblasts","Macrophages","Malignant cells","Monocytes"
cho <- "Malignant cells" 
a <- as.numeric(sumU[cho,which(sumU[cho,]>50)])
names(a) <- colnames(sumU)[which(sumU[cho,]>50)]
b <- as.numeric(sumD[cho,which(sumD[cho,]>50)])
names(b) <- colnames(sumD)[which(sumD[cho,]>50)]
order <- c(names(sort(b)),names(sort(a)))
datA <- data.frame(pw=names(a),value=a,group="a")
datB <- data.frame(pw=names(b),value=b,group="b")
dat <- rbind(datA,datB)

p <- dat %>%
  mutate(pw = factor(pw,levels=order))%>%
  ggplot()+
  geom_bar(aes(x=value,y=pw,fill=group),stat = 'identity')+
  # geom_bar(aes(x=value,y=pw),stat = 'identity',fill="#00467DFF")+
  labs(x="Pathway",y="Number of Samples")+
  scale_fill_manual(name="Metabolic capacity", values=c("#FFB900","#5773CC"))+
  theme_classic() + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 10,hjust=0.5,vjust=1),
        axis.text.y=element_text(colour="black", size = 10),
        axis.line=element_line(size=0.2,color="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"))+
  scale_x_break(c(5,50),scales = 20)
pdf(paste0("numofsamples_pw_",cho,".pdf"),width=12,height=4)
print(p)
dev.off()
png(file=paste0("numofsamples_pw_",cho,".png"),width=850,height=570,res=96)
print(p)
dev.off()

library(pheatmap)
library(paletteer)
p <- pheatmap(sumD,
         color = rev(paletteer_c("ggthemes::Red-Blue-White Diverging", 30)),
         border_color = "white",
         cluster_row = T,cluster_cols = F)
pdf(paste0("numofsamples_pw_heatmap.pdf"),width=12,height=6)
print(p)
dev.off()
