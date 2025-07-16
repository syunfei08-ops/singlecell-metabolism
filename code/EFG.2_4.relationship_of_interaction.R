source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/3.metabolic_landscape/interaction")

level <- "module" # "pathway","sub_pathway","module"
meta <- "high" # high, v1, v2, normal
A<-read_in_FC(level, meta)
ct <- sort(unique(A$celltype)) # 细胞类型统计

a <- read.table(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/fluxandcelltype_",meta,".list"))[,1,drop=T]
n <- length(a)

# 去除样本占比少的细胞类型
nSamCt <- c()
for (i in ct) {
  ct_sam <- A[which(A$celltype==i),]
  nSamCt <- c(nSamCt, length(unique(ct_sam$sample)))
}
names(nSamCt) <- ct
re <- ct[which(as.numeric(nSamCt) < 5)]
scFEA<-A[which(!A$celltype %in% c(re,"Unassigned")),]
ct <- levels(as.factor(scFEA$celltype))

matrNum <- matrix(0,nrow = length(ct),ncol = length(ct),dimnames = list(c(ct),c(ct)))
matrSum <- matrix(0,nrow = length(ct),ncol = length(ct),dimnames = list(c(ct),c(ct)))
for (i in 1:n){
  sam <- unlist(strsplit(a[i], ".", fixed = TRUE))[1]
  dir = paste("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/cellphoneDB/",sam,"_count_network.txt",sep="")
  intCount <- read.table(dir, header = T, sep = "\t")
  intCount <- intCount[which(!(intCount[,1] %in% c(re,"Unassigned") | intCount[,2] %in% c(re,"Unassigned"))),]
  for (j in 1:nrow(intCount)) {
    matrNum[intCount[j,1],intCount[j,2]] <- matrNum[intCount[j,1],intCount[j,2]]+1
    matrSum[intCount[j,1],intCount[j,2]] <- matrSum[intCount[j,1],intCount[j,2]]+intCount[j,3]
  }
}
matr <- as.data.frame(matrSum/matrNum)

library(reshape2)
library(Hmisc)
library(ggplot2)
# r1:EFG.2_2.通量分布规律的相关性的结果
# r1:EFG.2_3.转录组和代谢组的异同的结果
matr$ctX <- row.names(matr)
inter <- melt(matr,id.vars=c("ctX"),variable.name="ctY",value.name="interaction")
load("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/3.metabolic_landscape/overall/cor_metabolism.Rdata")
rM <- as.data.frame(r1)
rM$ctX <- row.names(rM)
rMet <- melt(rM,id.vars=c("ctX"),variable.name="ctY",value.name="cor")
load("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/3.metabolic_landscape/transcriptome/cor_transcriptome.Rdata")
rR <- as.data.frame(r1)
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
save(matr,rP,pP,rS,pS,file="interaction_met_rna.RData")

# dat <- as.data.frame(scale(matr))
# dat$x <- rownames(dat)
dat <- as.data.frame(matr)
dat <- dat[which(!is.nan(dat$inter)),]

col1<-"#2166AC"
col2<-"#7BCD7B" # 4393C3,7BCD7B
colS<-"#696773"
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
png(file=paste0("dotplot_MetInter.png"),width=850,height=570,res=96)
print(p)
dev.off()
pdf(paste0("dotplot_MetInter.pdf"),width=8,height=6)
print(p)
dev.off()
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
png(file=paste0("dotplot_RNAInter.png"),width=850,height=570,res=96)
print(p)
dev.off()
pdf(paste0("dotplot_RNAInter.pdf"),width=8,height=6)
print(p)
dev.off()
