# 使用相关性表示模块在样本中是否具有保守性
library(Hmisc) # rcorr函数
library(corrplot)
library(PerformanceAnalytics)
# library(reshape)
library(reshape2)
library(ggplot2)
library(paletteer)
library(ggrepel)
library(pheatmap)
library(magrittr) # %>%
library(dplyr) # mutate

source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/2-3.normal/KWtest")

# 定义函数：提取矩阵相关及其P值
CorMatrix <- function(cor,p) {
  ut <- upper.tri(cor)
  data.frame(row = rownames(cor)[row(cor)[ut]],
             column = rownames(cor)[col(cor)[ut]], 
             cor =(cor)[ut],
             p = p[ut] )
}

meta_data <- normal_meta()

level <- "module" # "pathway","sub_pathway","module"
meta <- "normal" # high, v1, v2, normal
scFEA_raw<-read_in_FC(level, meta)
scFEA <- ct_normal(scFEA_raw, meta_data) # 小于5个样本的细胞类型去掉

samN <- meta_data$sample[meta_data$class %in% "Normal"]
scFEA <- scFEA[scFEA$sample %in% samN,]

l <- levels(as.factor(scFEA$module_1))
# mod <- c("M_20","M_59","M_83","M_96","M_97","M_128")
# mod_2nd <- c("M_49","M_77","M_144","M_146","M_147")
# mod_ext <- setdiff(l,c(mod,mod_2nd))
sam <- levels(as.factor(scFEA$sample))
pval <- 0.01
rval <- 0.9
for (i in l) {
  rm(list = paste0("cor_",i))
}
percentage <- matrix(NA, nrow = length(l),ncol = 2, dimnames = list(l,c("pos","neg")))
for (i in l) {
  for (j in sam) {
    if(exists(paste0("cor_",i))){
      tem <- scFEA[which(scFEA$module_1 %in% i & scFEA$sample %in% j),c(4,5)]
      colnames(tem) <- c("celltype",j)
      assign(paste0("cor_",i),merge(get(paste0("cor_",i)),tem,by = "celltype", all=T))
    }else{
      tem <- scFEA[which(scFEA$module_1 %in% i & scFEA$sample %in% j),c(4,5)]
      colnames(tem) <- c("celltype",j)
      assign(paste0("cor_",i),tem)
    }
  }
  cor_ma <- get(paste0("cor_",i))[,-1]
  row.names(cor_ma) <- get(paste0("cor_",i))[,1]
  # 相关系数的显著性水平
  res_p <- rcorr(as.matrix(cor_ma), type = "pearson")
  re_pea <- CorMatrix(res_p$r, res_p$P)
  n_NaN <- nrow(re_pea[is.nan(re_pea$cor) | is.nan(re_pea$p),]) # NaN的原因是如果某个样本完全无变化
  n_all <- nrow(re_pea)
  re_pea <- re_pea[!(is.na(re_pea$cor) | is.na(re_pea$p)),] # na的原因是两个样本细胞类型无法重合
  n_pos <- length(which(re_pea$p<pval & re_pea$cor>rval))
  n_neg <- length(which(re_pea$p<pval & re_pea$cor<(-rval)))
  percentage[i,"pos"] <- (n_pos+n_NaN)/n_all*100
  percentage[i,"neg"] <- n_neg/n_all*100

  # mod_3rd <- setdiff(names(which(percentage[,1]>50)),c(mod,mod_2nd))
  
  # 画相关性点图和热图
  # re_pea$p <- -log10(re_pea$p)
  # re_pea$group <- 0
  # re_pea$group <- ifelse(re_pea$p>-log10(0.05),1,0)
  # p <- ggplot(data = re_pea,aes(x=cor,y=p))+
  #   geom_point(aes(fill=group),shape=21,colour="black")+
  #   geom_hline(aes(yintercept=-log10(0.05)),color="#999d9c",linetype="dashed",size=0.2)+
  #   geom_vline(aes(xintercept=0),color="#999d9c",linetype="dashed",size=0.2)+
  #   scale_fill_continuous(low = "#5773CC", high = "#FFB900")+
  #   theme_classic() +
  #   theme(legend.position="none",
  #         axis.text.x=element_text(colour="black", size = 6,angle=45,hjust=1,vjust=1),
  #         axis.text.y=element_text(colour="black", size = 6),
  #         axis.line=element_line(size=0.2,color="black"),
  #         axis.ticks = element_line(colour = "black",size=0.2),
  #         panel.border = element_blank(), panel.background = element_blank(),
  #         axis.ticks.length= unit(.5, "mm"))
  # png(file=paste(i,"_pearson_dot.png",sep = ""),width=850,height=570,res=96)
  # print(p)
  # # corrplot(apply(r1,2,as.numeric), order = "AOE", method = "color",tl.cex=0.5)
  # dev.off()
}

######
# 计算样本的sd
rm(sd_df)
for (i in l) {
  sd <- c()
  mb <- c()
  for (n in 1:length(sam)) {
    a <- sd(get(paste0("cor_",i))[!is.na(get(paste0("cor_",i))[,n+1]),n+1])
    # sd[n] <- a
    b <- abs(mean(get(paste0("cor_",i))[!is.na(get(paste0("cor_",i))[,n+1]),n+1]))
    if(b==0){
      sd[n] <- 0
      mb[n] <- b
    }else{
      sd[n] <- a/b
      mb[n] <- b
    }
  }
  if(exists("sd_df")){
    tem <- data.frame(model=i,sam=sam,sd=sd)
    sd_df <- rbind(sd_df,tem)
  }else{
    sd_df <- data.frame(model=i,sam=sam,sd=sd)
  }
}
# sd_df[is.na(sd_df$sd),3] <- 0
sd_mat <- dcast(sd_df,sam~model)
row.names(sd_mat) <- sd_mat[,1]
sd_mat <- sd_mat[,-1]

# sd <- c()
# for (n in 1:7) {
#   sd[n] <- mean(sd_df[which(sd_df$model %in% mod[n]),1])
# }

######
# 画柱状图
# 统计各模块sd为0的样本数
a <- matrix(ncol = 6,nrow = ncol(sd_mat))
len <- matrix(ncol = 2, nrow = ncol(sd_mat))
le <- c()
for (i in 1:ncol(sd_mat)) {
  a[i,] <- summary(sd_mat[,i])
  len[i,2] <- length(which(sd_mat[,i] == 0))
  len[i,1] <- colnames(sd_mat)[i]
  le[i] <- length(which(sd_mat[,i] == 0))
}
colnames(a) <- names(summary(dat$M_20))
rownames(a) <- colnames(sd_mat)
colnames(len) <- c("model","length")
names(le) <- colnames(sd_mat)

dat1 <- as.data.frame(len)
dat1$length <- as.numeric(dat1$length)
order <- names(sort(le))[(length(le)-19):length(le)]
# order <- colnames(dat)
mod <- names(sort(le[which(le>0.5*length(sam))]))
mod_2nd <- setdiff(order,c(mod))
dat1$group <- "unconservative"
dat1$group[which(dat1$model %in% mod)] <- "conservative"
dat1$percentage <- dat1$length/length(sam)*100

# 模块取sd小的样本数的top20
p <- dat1[which(dat1$model %in% order),] %>%
  mutate(model = factor(model,levels=order))%>%
  ggplot()+
  geom_bar(aes(x=model,y=percentage,fill=group),stat = 'identity')+
  scale_fill_manual(name="Conservative", values=c("#FFB900","#5773CC"))+
  labs(x="model",y="Percentage of conservative samples")+
  theme_classic() + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 10,hjust=0.5,vjust=1),
        axis.text.y=element_text(colour="black", size = 10),
        axis.line=element_line(size=0.2,color="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"))

pdf("similarity_bar.pdf",width=8.5,height=4)
print(p)
dev.off()
png("similarity_bar.png",width=8.5,height=4,unit="in",res=300)
print(p)
dev.off()
    
# write.table(sd_df,file = "sd_df.txt",sep = "\t",quote = F,row.names = F)
# ggplot(data = sd_df,aes(x=model,y=sd,fill=model))+
#   geom_violin(trim=F,size=0.2,show.legend = F,scale = "width", adjust =.1,outlier.shape = NA,)+
#   scale_color_paletteer_d("ggsci::default_locuszoom")+
#   # geom_jitter()+
#   coord_cartesian(ylim = c(0,0.05)*1.05)+
#   theme_classic() + 
#   theme(legend.position="none",
#         axis.text.x=element_text(colour="black", size = 6,angle=45,hjust=1,vjust=1),
#         axis.text.y=element_text(colour="black", size = 6),
#         axis.line=element_line(size=0.2,color="black"),
#         axis.ticks = element_line(colour = "black",size=0.2),
#         panel.border = element_blank(), panel.background = element_blank(),
#         axis.ticks.length= unit(.5, "mm"))

######
dat <- sd_mat[,c(rev(mod),rev(mod_2nd))]
# dat <- sd_mat

bk <- c(seq(0,1.99,by=0.01),seq(2,10,by=0.1))
# bk <- c(seq(-9,-0.1,by=0.01),seq(0,9,by=0.01))
p <- pheatmap(t(dat),
              # color = paletteer_c("grDevices::Blues 3", 30),
              # color = rev(paletteer_c("ggthemes::Red-Blue-White Diverging", 30)),
              # color=c(colorRampPalette(colors = c("#2E5A87FF","white","#A90C38FF"))(100)),
              color = c(colorRampPalette(colors = c("#2E5A87FF","white"))(length(bk)/2),colorRampPalette(colors = c("white","#A90C38FF"))(length(bk)/2)),
              legend_breaks=seq(0,10,2),
              breaks=bk,
              border_color = "white",
              # scale = "column",
              fontsize_row = 10,fontsize_col = 8,
              cellheight = 10,cellwidth = 7,
              cluster_row = F,cluster_cols = F,
              # clustering_distance_rows = "euclidean",
              # clustering_method="single",
              gaps_row = c(7))

pdf("similarity_pheatmap.pdf",width=8,height=5)
print(p)
dev.off()
png("similarity_pheatmap.png",width=8,height=5,unit="in",res=300)
print(p)
dev.off()

######
# 画正相关和负相关的趋势图
data <- as.data.frame(percentage)
data$group <- "a"
data[mod,"group"] <- "b"
data[mod_2nd,"group"] <- "c"
data["M_77","group"] <- "c"
data$label <- NA
data[mod,"label"] <- mod
data[mod_2nd,"label"] <- mod_2nd
p<-ggplot(data,aes(x=neg,y=pos))+
  geom_point(aes(fill=group),shape=21,colour="black",size=3,stroke = 1)+
  geom_text_repel(label=data$label,check_overlap = F,
            hjust = 0, nudge_x = 0.05)+
  scale_fill_manual(values=c("a" = "#E5E5E5", "b" = "#D82632", "c" = "#BB00BB"))+
  theme_classic() + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 10,hjust=1,vjust=1),
        axis.text.y=element_text(colour="black", size = 10),
        axis.line=element_line(size=0.2,color="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"))
pdf("similarity_dot.pdf",width=6,height=4)
print(p)
dev.off()
png("similarity_dot.png",width=6,height=4,unit="in",res=300)
print(p)
dev.off()
