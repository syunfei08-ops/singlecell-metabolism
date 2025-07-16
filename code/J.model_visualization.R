source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/5.model/output")

level <- "module" # "pathway","module"
for(n in 1:8){

load(paste0("round",n,"_allOutput",level,".Rdata"))

mod = list.files(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/5.model/round",n))
dir = paste("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/5.model/round",n,"/",mod,sep="")
n1 <- length(mod)

# method <- c('treebag','bagEarth','bagEarthGCV','C5.0','rpart','rpart1SE','rpart2',
            # 'rpartScore','fda','avNNet','earth','gcvEarth','pcaNNet','pda',
            # 'multinom','ranger','C5.0Rules','C5.0Tree','evtree','AdaBag','parRF',
            # 'rFerns','ordinalRF','rf','rfRules','RRF','RRFglobal','wsrf','naive_bayes',
            # 'LogitBoost','bstTree','sparseLDA','sdwd','kernelpls','pls','glmnet',
            # 'svmLinear3','nnet','RFlda','gpls','mlp','mlpWeightDecay',
            # 'snn','gaussprRadial','svmRadialWeights')
method <- c("parRF","treebag","ranger","rf","gaussprRadial","RRFglobal")

# 找到有问题的模型
for(i in 1:n1){
  for(k in 1:3)
    if(re[[mod[i]]][[k]] %in% "not enough"){
      print(paste0(i,"-",j,"-",k))
    }
}

# 统计100次结果的平均值
# reM_rf <- list()
# reM_tb <- list()
# reM <- list()
# mod <- list.files(dir[i])
reM_rf <- matrix(NA,nrow = length(mod),ncol = 8,
                          dimnames = list(mod,c("method","Accuracy_train","Kappa_train","Accuracy","Kappa","Precision","Recall","F1 score")))
reM_tb <- matrix(NA,nrow = length(mod),ncol = 8,
                          dimnames = list(mod,c("method","Accuracy_train","Kappa_train","Accuracy","Kappa","Precision","Recall","F1 score")))
reM <- matrix(NA,nrow = length(mod),ncol = 8,
                          dimnames = list(mod,c("method","Accuracy_train","Kappa_train","Accuracy","Kappa","Precision","Recall","F1 score")))
for (j in 1:n1) {
  for (k in 1:3) {
    re[[mod[j]]][[k]][re[[mod[j]]][[k]]<0] <- 0
    re[[mod[j]]][[k]][is.na(re[[mod[j]]][[k]])] <- 0
  }
  reM_rf[j,2:8] <- colMeans(apply(re[[mod[j]]][[1]][,2:8],2,as.numeric))
  reM_rf[j,1] <- mod[j]
  reM_tb[j,2:8] <- colMeans(apply(re[[mod[j]]][[2]][,2:8],2,as.numeric))
  reM_tb[j,1] <- mod[j]
  reM[j,2:8] <- colMeans(apply(re[[mod[j]]][[3]][,2:8],2,as.numeric))
  reM[j,1] <- mod[j]
}
save(reM_rf, reM_tb, reM, file=paste0("round",n,"_",level,"mean100.Rdata"))

}
######
# 各模型结果比较
library(corrplot)
library(paletteer)

level <- "module" # "pathway","module"
load(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/5.model/output/round7_",level,"mean100.Rdata"))

reM <- as.data.frame(reM)
reM[2:8] <- apply(reM[,2:8],2,as.numeric)
a <- which(is.nan(reM[,"F1 score"]))
reM[a,"F1 score"] <- reM[a,"Precision"]*reM[a,"Recall"]*2 / (reM[a,"Precision"]+reM[a,"Recall"])

reM_rf <- as.data.frame(reM_rf)
reM_rf[2:8] <- apply(reM_rf[,2:8],2,as.numeric)
a <- which(is.nan(reM_rf[,"F1 score"]))
reM_rf[a,"F1 score"] <- reM_rf[a,"Precision"]*reM_rf[a,"Recall"]*2 / (reM_rf[a,"Precision"]+reM_rf[a,"Recall"])

reM_tb <- as.data.frame(reM_tb)
reM_tb[2:8] <- apply(reM_tb[,2:8],2,as.numeric)
a <- which(is.nan(reM_tb[,"F1 score"]))
reM_tb[a,"F1 score"] <- reM_tb[a,"Precision"]*reM_tb[a,"Recall"]*2 / (reM_tb[a,"Precision"]+reM_tb[a,"Recall"])

# 文件结构
# reM_rf, reM_tb, reM，分成四个子列表，每个子列表为数据组合的结果矩阵
choCol <- c(4,5,6,7,8)

matr_rf <- reM_rf[,choCol]
matr_tb <- reM_tb[,choCol]
matr <- reM[,choCol]
for (i in c("matr_rf","matr_tb","matr")) {
  dat <- apply(get(i),2,as.numeric)
  row.names(dat) <- row.names(get(i))
  # row.names(dat) <- a
  dat <- dat[order(dat[,1],decreasing = TRUE),]
  pdf(paste0(paste0(i,"_Status.pdf")),width=18,height=10)
  corrplot(t(dat), is.corr = FALSE, 
           method = "color", outline = "white",
           col = colorRampPalette(colors = c("#0da9ce","white","#e74a32"))(100),
           addCoef.col = "#727D84",
           tl.col = "black",tl.cex = 0.9,tl.srt = 90,
           number.digits = 2,number.cex = 0.9,number.font = NULL,
           cl.pos = "n")
  dev.off()
}


######
# 使用tumor_module的数据的Fit21.Rdata(parRF)结果最好
library(caret)
library(dplyr)
load(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/5.model/round7/tumor_module_Fit1.Rdata"))
# 特征权重
reVar <- matrix(NA, nrow = 100, ncol = length(pro[[1]]),dimnames = list(c(1:100),pro[[1]]))
for (m in 1:100) {
  tem <- varImp(l1[[m]])
  reVar[m,] <- t(tem[[1]])
}
var <- as.data.frame(colMeans(reVar))
var$name <- row.names(var)
colnames(var) <- c("weight","name")


p2 <-var %>%
  # mutate(name = factor(name,levels=mixedsort(var$name,decreasing = TRUE)))%>%
  mutate(name = factor(name,levels=var$name[order(var$weight)]))%>%
  ggplot(aes(x=weight,y=name))+
    geom_col(fill="#009999")+
    theme_classic() + 
    theme(legend.position="none",
          axis.text.x=element_text(colour="black", size =12,angle=0,hjust=1,vjust=1),
          axis.text.y=element_text(colour="black", size = 12),
          axis.line=element_line(size=0.2,color="black"),
          axis.ticks = element_line(colour = "black",size=0.2),
          panel.border = element_blank(), panel.background = element_blank(),
          axis.ticks.length= unit(.5, "mm"))
pdf(paste0(paste0("feature","Weight.pdf")),width=6,height=10)
print(p2)
dev.off()

library(reshape2)
# accuracy,kappa,precision,recall的分布
# list1为4个数据组合的文件夹，list2为39种预测模型，list3为3种特征提取的结果，list4为100次循环的结果
load(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/5.model/output/round7_allOutputmodule.Rdata"))
mat <- re[["tumor_module_Fit1.Rdata"]][[1]][,-1]
mat <- as.data.frame(apply(mat,2,as.numeric))
mat$name <- 1:100
data <- melt(mat,
             id.vars = c("name"),
             measure.vars = c("Accuracy","Kappa","Precision","Recall","F1 score"),
             variable.name = "Measure",
             value.name = "value")
data$value[which(is.na(data$value))] <- 0 
data$value[which(data$value<0)] <- 0 

color <- c("#CAB2D6FF","#FB9A99FF","#18BC9CFF","#1F78B4FF","#B2DF8AFF")
# color <- c("#CD534C","#EFC000","#4A6990","#1F77B4","#9467BD")
p3 <- ggplot(data, aes(x=Measure, y=value, fill=Measure))+
# ggplot(data, aes(x=Measure, y=value, fill=Measure))+
  geom_boxplot(aes(fill=Measure),color="black",
               outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)+
  # stat_summary(fun.y = mean,geom="point",size=1,color="#0080FF")+
  theme_classic() + 
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size =12,angle=0,hjust=1,vjust=1),
        axis.text.y=element_text(colour="black", size = 12),
        axis.line=element_line(size=0.2,color="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"))+
  scale_fill_manual(name="Cell type", values=color)+
  scale_color_manual(name="Cell type", values=color)
pdf(paste0(paste0("measure",".pdf")),width=10,height=5)
print(p3)
dev.off()

######
# 所选特征在转移、非转移组的分布情况
library(reshape2)
library(dplyr)
library(ggplot2)
l<-read.table("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/status.meta")
load(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/5.model/round7/tumor_module_Fit1.Rdata"))
tumStage <- as.factor(l[,2])
tumSample <- l[,1]

scFEA_raw<-read_in_FC("module","status")
scFEA <- scFEA_raw[scFEA_raw$celltype %in% "Malignant cells",c(2,3,5)]
colnames(scFEA) <- c("sample","module","log2FC")

scFEA_stage <- scFEA[scFEA$sample %in% tumSample,]
scFEA_stage <- scFEA_stage[scFEA_stage$module %in% pro[[1]],]
dataStage0 <- acast(scFEA_stage,sample~module)
dataStage <- as.data.frame(cbind(tumStage,dataStage0))
dat <- melt(dataStage,
            id.vars = c("tumStage"),
            variable.name = "feature",
            value.name = "value")
p4 <- dat %>%
  mutate(tumStage = factor(tumStage,levels=c("1","2")))%>%
  # mutate(feature = factor(feature,levels=mixedsort(pro[[1]])))%>% 
  mutate(feature = factor(feature,levels=rev(var$name[order(var$weight)])))%>% 
  ggplot(aes(x=feature, y=value, fill=tumStage)) +
    geom_boxplot()+
    coord_cartesian(ylim = c(-1,3))+
    scale_fill_manual(name="Cell type", values=c("#e58c33","#1fc2dd"))+
    theme_classic() + 
    theme(legend.position="right",
          axis.text.x=element_text(colour="black", size =12,angle=0,hjust=1,vjust=1),
          axis.text.y=element_text(colour="black", size = 12),
          axis.line=element_line(size=0.2,color="black"),
          axis.ticks = element_line(colour = "black",size=0.2),
          panel.border = element_blank(), panel.background = element_blank(),
          axis.ticks.length= unit(.5, "mm"))
pdf(paste0(paste0("featureDistribution",".pdf")),width=15,height=6)
print(p4)
dev.off()

######
# 绘制模型的ROC曲线以及AUC的计算、
load(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/5.model/round7/tumor_module_Fit1.Rdata"))
library(pROC)
library(data.table)

val <- status_val("all_stage_sample")
tumStage <- as.factor(val[[1]])
tumSample <- val[[2]]
scFEA_raw<-read_in_FC("module","status")
scFEA <- scFEA_raw[scFEA_raw$celltype %in% "Malignant cells",c(2,3,5)]
colnames(scFEA) <- c("sample","module","log2FC")
scFEA_stage <- scFEA[scFEA$sample %in% tumSample,]
dataStage0 <- acast(scFEA_stage,sample~module)

sumAuc <- 0
sumSen <- 0
sumSpe <- 0
#在得到最终模型后，算出预测值
dat <- matrix(NA,nrow = 8,ncol = 2)
color <- colorRampPalette(colors = c("#0da9ce","white","#f49537"))(100)
for (i in 1:100) {
  # 还原测试集
  opt_col <- colnames(l1[[i]]$trainingData)[-length(l1[[i]]$trainingData)]
  opt_row <- rownames(l1[[i]]$trainingData)
  dataStage <- dataStage0[,opt_col]
  trainx = dataStage[opt_row,]
  testx = dataStage[-which(rownames(dataStage) %in% opt_row),]
  testx <- apply(testx,2,as.numeric)
  trainy = tumStage[which(rownames(dataStage) %in% opt_row)]
  testy = tumStage[-which(rownames(dataStage) %in% opt_row)]
  
  
  # 获得预测结果
  pred <- predict(l1[[i]], newdata = testx)
  # ROC, 输入真实的label和预测的值，ci是为了算95%CI（如果是做某个指标用于预测死亡率的ROC曲线，那就用是否死亡代替Y_test, 用这个指标的值代替pred）
  roc <- roc(testy, as.numeric(pred), ci = TRUE)
  auc <- auc(roc)
  
  # 把数据整理出来，用ggplot2再作图
  if(i==1){
    roc_data <- data.table(sensitivity = roc$sensitivities,
                           specificity = roc$specificities,
                           yorden = roc$sensitivities + roc$specificities - 1, # 约登指数
                           auc = roc$auc,
                           auc_low = as.numeric(roc$ci)[1],
                           auc_up = as.numeric(roc$ci)[3])
    roc_data$group <- i
    sen <- roc_data$sensitivity[which.max(roc_data$yorden)]
    spe <- roc_data$specificity[which.max(roc_data$yorden)]
  }else{
    tem <- data.table(sensitivity = roc$sensitivities,
                      specificity = roc$specificities,
                      yorden = roc$sensitivities + roc$specificities - 1, # 约登指数
                      auc = roc$auc,
                      auc_low = as.numeric(roc$ci)[1],
                      auc_up = as.numeric(roc$ci)[3])
    tem$group <- i
    sen <- tem$sensitivity[which.max(tem$yorden)]
    spe <- tem$specificity[which.max(tem$yorden)]
    roc_data <- rbind(roc_data,tem)
  }
  sumAuc <- sumAuc+auc
  sumSen <- sumSen+sen
  sumSpe <- sumSpe+spe
  
  # 获得预测概率
  
  dat[,1] <- ifelse(testy=="Primary",1,0)
  pred <- data.frame(predict(l1[[i]],type="prob", newdata=testx),
                     predicted=predict(l1[[i]], newdata=testx))
  dat[,2] <- pred$Primary
  ## 画在一起
  # if(i==1){
  #   x1<-plot.roc(dat[,1],dat[,2],direction = "<",smooth=F,lwd=2,
  #                ylim=c(0,1),xlim=c(1,0),legacy.axes=T,main="",col=color[i])
  #   sumAuc <- sumAuc+x1$auc
  # }else{
  #   x1<-plot.roc(dat[,1],dat[,2],direction = "<",smooth=F,lwd=2,
  #                ylim=c(0,1),xlim=c(1,0),legacy.axes=T,main="",col=color[i],add=TRUE)
  #   sumAuc <- sumAuc+x1$auc
  # }
  # 单独画
  pdf(paste0(paste0("ROC",i,".pdf")),width=7,height=7)
  jpeg(paste0("ROC",i,".png"),width=2000,height=2000,res=430);
  x1<-plot.roc(dat[,1],dat[,2],direction = "<",smooth=F,lwd=2,
               ylim=c(0,1),xlim=c(1,0),legacy.axes=T,main="",col=color[i])
  dev.off()
  
}
# 测试
# ROC曲线
# 方法1
# p5 <- ggplot()+
#   geom_segment(
#   aes(x=0, y=0, xend=1, yend=1),
#   linetype = "dotted",
#   color = "grey50"
# ) +
#   xlab("False Positive Rate") +
#   ylab("True Positive Rate") +
#   ggtitle("ROC Curve")+ 
#   theme_classic() + 
#   theme(legend.position="none",
#         axis.text.x=element_text(colour="black", size =12,angle=0,hjust=1,vjust=1),
#         axis.text.y=element_text(colour="black", size = 12),
#         axis.line=element_line(size=0.2,color="black"),
#         axis.ticks = element_line(colour = "black",size=0.2),
#         panel.border = element_blank(), panel.background = element_blank(),
#         axis.ticks.length= unit(.5, "mm"))+
#   geom_line(data=roc_data[which(roc_data$group==1),],
#             aes(x = 1-specificity, y=sensitivity),
#             linetype="solid",linewidth=1.2,
#             # color = color[as.numeric(1)],alpha = 0.3)
#             color = "black",alpha = 0.1)
# # p1 <- p+geom_line(data=roc_data[which(roc_data$group==2)],aes(x = 1-specificity, y=sensitivity),
# #                   color = color[as.numeric(roc_data$group[2])])
# for (i in 2:100) {
#   p5 <- p5+geom_line(data=roc_data[which(roc_data$group==i)],
#                      aes(x = 1-specificity, y=sensitivity),
#                      linetype="solid",linewidth=1.2,
#                      # color = color[as.numeric(i)],alpha = 0.1)
#                      color = "black",alpha = 0.3)
# }
# 
# p5 <- p5 +
#   annotate("text",
#            x = 0.7, y = 0.2,
#            label = paste0('mean(AUC) = ', round(sumAuc/100, 3))) +
#   annotate("text",
#            x = 0.7, y = 0.12,
#            label = paste0('mean(Sensitivity) = ', round(sumSen/100, 3))) +
#   annotate("text",
#            x = 0.7, y = 0.05,
#            label = paste0('mean(Specificity) = ', round(sumSpe/100, 3)))
# pdf(paste0(paste0("ROC",".pdf")),width=15,height=7)
# print(p5)
# dev.off()
