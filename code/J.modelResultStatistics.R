source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/5.model/output")
# 样本信息
l<-read.table("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/status.meta")
tumStage <- as.factor(l[,2])
tumSample <- l[,1]
# method <- c('treebag','bagEarth','bagEarthGCV','C5.0','rpart','rpart1SE','rpart2',
            # 'rpartScore','fda','avNNet','earth','gcvEarth','pcaNNet','pda',
            # 'multinom','ranger','C5.0Rules','C5.0Tree','evtree','AdaBag','parRF',
            # 'rFerns','ordinalRF','rf','rfRules','RRF','RRFglobal','wsrf','naive_bayes',
            # 'LogitBoost','bstTree','sparseLDA','sdwd','kernelpls','pls','glmnet',
            # 'svmLinear3','nnet','RFlda','gpls','mlp','mlpWeightDecay',
            # 'snn','gaussprRadial','svmRadialWeights')
method <- c("parRF","treebag","ranger","rf","gaussprRadial","RRFglobal")

library(caret)
library(reshape2)

level <- "module" # "pathway","module"
meta <- "high" # high, v1, v2, normal

# 数据读取
scFEA_raw<-read_in_FC(level, meta)
scFEA <- scFEA_raw[scFEA_raw$celltype %in% "Malignant cells",c(2,3,5)]
colnames(scFEA) <- c("sample","module","log2FC")
scFEA_stage <- scFEA[scFEA$sample %in% tumSample,]
dataStage0 <- acast(scFEA_stage,sample~module)

for(i in 1:8){

group = list.files(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/5.model/round",i))
dir = paste("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/5.model/round",i,"/",group,sep="")
n1 <- length(group)

re <- list()
for (j in 1:n1) { # 预测模型数量
  load(dir[j])
  # pro保存特征
  # l1, l2, l3保存rfFuncs, treebagFuncs和总体特征的100次循环结果
  matr <- list()
  for(k in 1:3){ # 3个特征层级
    l <- get(paste0("l",k))
    mat100 <- matrix(NA,nrow = 100,ncol = 8,
                     dimnames = list(1:100,c("parameter","Accuracy_train","Kappa_train","Accuracy","Kappa","Precision","Recall","F1 score")))
    if(length(l)<100){
      matr[[k]] <- "not enough"
    }else{
      for(m in 1:100){ # 100次
        print(paste0(group[j],";",k,";",m))
        # for(m in 16:18){
        if(is.null(l[[m]])) 
          next
        # 还原输入数据
        opt_col <- colnames(l[[m]]$trainingData)[-length(l[[m]]$trainingData)]
        opt_row <- rownames(l[[m]]$trainingData)
        dataStage <- dataStage0[,opt_col]
        trainx = dataStage[opt_row,]
        testx = dataStage[-which(rownames(dataStage) %in% opt_row),]
        testx <- apply(testx,2,as.numeric)
        trainy = tumStage[which(rownames(dataStage) %in% opt_row)]
        testy = tumStage[-which(rownames(dataStage) %in% opt_row)]
         
        # 训练参数及结果
        para <- paste(l[[m]]$bestTune,collapse = ";")
        mat100[m,1] <-  para
        for (b in 1:nrow(l[[m]]$results)) {
          paraM <- paste(l[[m]]$results[b,1:length(l[[m]]$bestTune)],collapse = ";")
          if(paraM==para){
            n3=b
          }
        }
        mat100[m,2] <- l[[m]]$results[n3,]$Accuracy
        mat100[m,3] <- l[[m]]$results[n3,]$Kappa
        
        # 预测结果
        pred <- predict(l[[m]],newdata = testx)
        
        cM <- confusionMatrix(pred, testy)
        mat100[m,4] <- cM$overall["Accuracy"]
        mat100[m,5] <- cM$overall["Kappa"]
        mat100[m,6] <- cM$byClass["Precision"]
        mat100[m,7] <- cM$byClass["Recall"]
        mat100[m,8] <- cM$byClass["F1"]
      }  
    }
	matr[[k]] <- mat100
  }
  re[[group[j]]] <- matr
}

save(re,file=paste0("round",i,"_allOutput",level,".Rdata"))

}
# # 统计100次结果的平均值
# reM_rf <- list()
# reM_tb <- list()
# reM <- list()
# for (i in 1:n1) {
  # mod <- list.files(dir[i])
  # reM_rf[[group[i]]] <- matrix(NA,nrow = length(mod),ncol = 8,
                            # dimnames = list(mod,c("method","Accuracy_train","Kappa_train","Accuracy","Kappa","Precision","Recall","F1 score")))
  # reM_tb[[group[i]]] <- matrix(NA,nrow = length(mod),ncol = 8,
                            # dimnames = list(mod,c("method","Accuracy_train","Kappa_train","Accuracy","Kappa","Precision","Recall","F1 score")))
  # reM[[group[i]]] <- matrix(NA,nrow = length(mod),ncol = 8,
                            # dimnames = list(mod,c("method","Accuracy_train","Kappa_train","Accuracy","Kappa","Precision","Recall","F1 score")))
  # for (j in 1:length(re[[group[i]]])) {
    # for (k in 1:3) {
      # re[[group[i]]][[mod[j]]][[k]][re[[group[i]]][[mod[j]]][[k]]<0] <- 0
      # re[[group[i]]][[mod[j]]][[k]][is.na(re[[group[i]]][[mod[j]]][[k]])] <- 0
    # }
    # reM_rf[[group[i]]][j,2:8] <- colMeans(apply(re[[group[i]]][[mod[j]]][[1]][,2:8],2,as.numeric))
    # reM_rf[[group[i]]][j,1] <- mod[j]
    # reM_tb[[group[i]]][j,2:8] <- colMeans(apply(re[[group[i]]][[mod[j]]][[2]][,2:8],2,as.numeric))
    # reM_tb[[group[i]]][j,1] <- mod[j]
    # reM[[group[i]]][j,2:8] <- colMeans(apply(re[[group[i]]][[mod[j]]][[3]][,2:8],2,as.numeric))
    # reM[[group[i]]][j,1] <- mod[j]
  # }
# }
