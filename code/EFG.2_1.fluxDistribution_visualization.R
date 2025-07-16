######
# 特别长的大小提琴图
# 通路和子通路是否保守，相对于代码：课题_细胞类型统计和样本代谢趋势表示
source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")

library(ggplot2)
options(stringsAsFactors=FALSE)
# 对通路的保守性进行评估，使用KW检验
# 读入数据
level <- "module" # "pathway","module"
meta <- "high" # high, v1, v2, normal
scFEA_raw<-read_in_FC(level, meta)
scFEA <- ct_adj(scFEA_raw)

######
# 统计数据
ct <- levels(as.factor(scFEA[,4]))
pw <- levels(as.factor(scFEA[,3]))

# 各细胞类型的总体平均值
dim <- list(c(ct),c("median","mean",pw))
re_ct <- matrix(nrow = length(ct),ncol = length(pw)+2,dimnames = dim)
for (i in ct) {
  re_ct[i,"median"] <- median(scFEA[which(scFEA[,4]==i),5])
  re_ct[i,"mean"] <- mean(scFEA[which(scFEA[,4]==i),5])
  for (j in pw) {
    re_ct[i,j] <- mean(scFEA[which(scFEA[,4]==i & scFEA[,3]==j),5])
  }
}
# 各通路的代谢水平
dim <- list(c(pw),c("median","mean",ct))
re_pw <- matrix(nrow = length(pw),ncol = length(ct)+2,dimnames = dim)
for (i in pw) {
  re_pw[i,"median"] <- median(scFEA[which(scFEA[,3]==i),5])
  re_pw[i,"mean"] <- mean(scFEA[which(scFEA[,3]==i),5])
  for (j in ct) {
    re_pw[i,j] <- mean(scFEA[which(scFEA[,3]==i & scFEA[,4]==j),5])
  }
}
write.table(re_pw,file = "/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/3.metabolic_landscape/overall/celltype_pathway_FC.txt",
            sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
######
# 画图
color <- color_ct_val()

# 通路名称更改
# 字符串切割
# swr = function(string, nwrap=20) {
  # paste(strwrap(string, width=nwrap), collapse="\n")
# }
# swr = Vectorize(swr)

# 字符串排序
# scFEA$pathway_2 = swr(scFEA$pathway_2)
# scFEA[,3] <- factor(scFEA[,3], levels = unique(scFEA[,3]))

filename <- paste0(level,"_distribution.jpg")
jpeg(filename=paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/3.metabolic_landscape/overall/",filename),width=300+(200*length(levels(as.factor(scFEA[,3])))-1),height=1000,res=160)
p <- ggplot(scFEA, aes(x=scFEA[,5], y=scFEA[,4], fill=scFEA[,4]))+
  geom_violin(aes(color=scFEA[,4],fill=scFEA[,4]))+
  facet_grid(.~scFEA[,3]) + # 划分图像
  labs(y=NULL,x="Log2(FoldChange)") +
  coord_cartesian(xlim = c(-5,5))+
  stat_summary(fun.y = median,geom="point",size=1,color="black")+
  geom_vline(aes(xintercept=0),color="#999d9c",linetype="dashed")+
  theme(legend.position="none",
          panel.border = element_blank(),
          panel.background = element_rect(fill = "white" ),
          panel.grid.major = element_line(color = '#EAEBEB'),
          panel.grid.minor = element_line(color = '#EAEBEB'),
          axis.text.x=element_text(colour="black",angle=0,hjust=0.5,size = 5),
          axis.text.y=element_text(colour="#222222",angle=0,hjust=1,size = 11),
          strip.text.x = element_text(size = 13, color = "#1E1E1E"),
          strip.background = element_rect(
            color="#CFD1D0", fill="white", size=1.5, linetype="solid")
          )+
  scale_fill_manual(name="Cell type", values=color)+
  scale_color_manual(name="Cell type", values=color)
print(p) # dev.off无法在循环中保存图片，加上print，解决问题
dev.off()


