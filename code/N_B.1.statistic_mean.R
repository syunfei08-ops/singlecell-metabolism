source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/2-3.normal/statistic")

# 构建样本信息
meta_data <- normal_meta() 

level <- "pathway" # "pathway","sub_pathway","module"
meta <- "normal" # high, v1, v2, normal
scFEA_raw<-read_in_FC(level, meta)
scFEA <- ct_normal(scFEA_raw,meta_data) # 小于5个样本的细胞类型去掉

ct <- levels(as.factor(scFEA$celltype))
type <- levels(as.factor(meta_data$class))
ds <- levels(as.factor(meta_data$dataset))
tum <- levels(as.factor(meta_data$tumor))

pw <- levels(as.factor(scFEA[,3]))

sam <- meta_data$sample[which(meta_data$class=="Normal")]
n_scFEA <- scFEA[scFEA$sample %in% sam,]
meta <- meta_data[meta_data$class %in% "Normal",]

# 上下调的通路/模块
for (t in tum) {
  sCho <- meta[which(meta[,"tumor"] %in% t),1]
  t_scFEA <- n_scFEA[which(n_scFEA$sample %in% sCho),]

  ct <- levels(as.factor(t_scFEA[,4]))
  sam <- levels(as.factor(t_scFEA[,2]))
  Nct <- as.numeric(table(t_scFEA$celltype)/length(pw))
  
  sumU <- matrix(0,nrow = length(ct),ncol = length(pw),dimnames = list(ct,pw))
  sumD <- matrix(0,nrow = length(ct),ncol = length(pw),dimnames = list(ct,pw))
  for (i in ct) {
    for (j in pw) {
      a <- t_scFEA[t_scFEA[,1] %in% j & t_scFEA$celltype %in% i,]
      sumU[i,j] <- length(which(a$log2FC>0))
      sumD[i,j] <- length(which(a$log2FC<0))
    }
  }
  sumU <- sumU/Nct*100
  sumD <- sumD/Nct*100
  # t(sumU)、t(sumD)，看上下调的通路/模块
  write.table(sumU,file=paste0("normal-",t,"_",level,"_mean_up.txt"),sep="\t")
  write.table(sumD,file=paste0("normal-",t,"_",level,"_mean_down.txt"),sep="\t") # log2FC < 0 in more than 50% of samples
  
  re <- matrix(NA,nrow = length(ct),ncol = 7,dimnames = list(ct,c("# Total samples","50%","80%","Upregulated >= 80% samples",
             "Upregulated >= 50% samples","Downregulated >= 80% samples","Downregulated >= 50% samples")))
  re[,1] <- Nct
  re[,2] <- 0.5*Nct
  re[,3] <- 0.8*Nct
  for (r in 1:nrow(sumU)) {
    a <- length(which(sumU[r,]>50))
    ba <- round(as.numeric(a)/as.numeric(length(pw))*100,2)
    b <- length(which(sumU[r,]>80))
    bb <- round(as.numeric(b)/as.numeric(length(pw))*100,2)
    re[r,5] <- paste0(a,"(",ba,"%)")
    re[r,4] <- paste0(b,"(",bb,"%)")
  }
  for (r in 1:nrow(sumD)) {
    a <- length(which(sumD[r,]>50))
    ba <- round(as.numeric(a)/as.numeric(length(pw))*100,2)
    b <- length(which(sumD[r,]>80))
    bb <- round(as.numeric(b)/as.numeric(length(pw))*100,2)
    re[r,7] <- paste0(a,"(",ba,"%)")
    re[r,6] <- paste0(b,"(",bb,"%)")
  }
  write.table(re,paste0("normal-",t,"_",level,"_up_down_statistic.txt"),sep = "\t",quote = F)
}

# ######
# # 柱状图，显示已选通路的样本数
# library(ggplot2)
# library(magrittr) # %>%
# library(dplyr) # mutate
# library(ggbreak)

# # "Dendritic cells","Fibroblasts","Macrophages","Malignant cells","Monocytes"
# cho <- "Malignant cells" 
# for (t in tum) {
  # sumU <- read.table(paste0(t,"_mean_up.txt"),sep="\t")
  # sumD <- read.table(paste0(t,"_mean_down.txt"),sep="\t")
  
  # a <- sumU[cho,which(sumU[cho,]>50)]
  # names(a) <- colnames(sumU)[which(sumU[cho,]>50)]
  # b <- sumD[cho,which(sumD[cho,]>50)]
  # names(b) <- colnames(sumD)[which(sumD[cho,]>50)]
  # order <- c(names(sort(b)),names(sort(a)))
  # datA <- data.frame(pw=names(a),value=a,group="a")
  # datB <- data.frame(pw=names(b),value=b,group="b")
  # dat <- rbind(datA,datB)
 
  # p <- dat %>%
    # mutate(pw = factor(pw,levels=order))%>%
    # ggplot()+
    # geom_bar(aes(x=value,y=pw,fill=group),stat = 'identity')+
    # # geom_bar(aes(x=value,y=pw),stat = 'identity',fill="#00467DFF")+
    # labs(x="Pathway",y="Number of Samples")+
    # scale_fill_manual(name="Metabolic capacity", values=c("#FFB900","#5773CC"))+
    # theme_classic() + 
    # theme(legend.position="none",
          # axis.text.x=element_text(colour="black", size = 10,hjust=0.5,vjust=1),
          # axis.text.y=element_text(colour="black", size = 10),
          # axis.line=element_line(size=0.2,color="black"),
          # axis.ticks = element_line(colour = "black",size=0.2),
          # panel.border = element_blank(), panel.background = element_blank(),
          # axis.ticks.length= unit(.5, "mm"))+
    # scale_x_break(c(5,50),scales = 20)
  # pdf(paste0(t,"_numofsamples_pw_",cho,".pdf"),width=12,height=4)
  # print(p)
  # dev.off()
  # png(file=paste0(t,"_numofsamples_pw_",cho,".png"),width=850,height=570,res=96)
  # print(p)
  # dev.off()
# }

######
# 仅在某一细胞类型特异性在平均水平上/下
for (t in tum) {
  sCho <- meta[which(meta[,"tumor"] %in% t),1]
  t_scFEA <- n_scFEA[which(n_scFEA$sample %in% sCho),]
  
  ct <- levels(as.factor(t_scFEA[,4]))
  sam <- levels(as.factor(t_scFEA[,2]))
  Nct <- as.numeric(table(t_scFEA$celltype)/length(pw))

  sumSu <- matrix(0,nrow = length(ct),ncol = length(pw),dimnames = list(ct,pw))
  sumSd <- matrix(0,nrow = length(ct),ncol = length(pw),dimnames = list(ct,pw))
  for (s in sam) {
    for (j in pw) {
      a <- scFEA[scFEA$module_1 %in% j & scFEA$sample %in% s,]
      if(length(which(a$log2FC>0))==1){
        i <- a$celltype[which(a$log2FC>0)]
        sumSu[i,j] <- sumSu[i,j]+1
      }
      if(length(which(a$log2FC<0)) == 1){
        i <- a$celltype[which(a$log2FC<0)]
        sumSd[i,j] <- sumSd[i,j]+1
      }
    }
  }
  sumSu <- sumSu/Nct*100
  sumSd <- sumSd/Nct*100
  # sumSu <- read.table(paste0(t,"_exclusive_up.txt"),sep="\t")
  # sumSd <- read.table(paste0(t,"_exclusive_down.txt"),sep="\t")
  for (r in 1:length(ct)) {
    print(paste0(ct[r],":",colnames(sumSu)[which(sumSu[ct[r],]>40)])) # log2FC > 0 exclusively in malignant cells in more than 50% of samples
    print(paste0(ct[r],":",colnames(sumSd)[which(sumSd[ct[r],]>40)]))
  }
  write.table(sumSu,file=paste0("normal-",t,"_",level,"_exclusive_up.txt"),sep="\t")
  write.table(sumSd,file=paste0("normal-",t,"_",level,"_exclusive_down.txt"),sep="\t")
}
