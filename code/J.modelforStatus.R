source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
# setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/4.downstream_analysis/moduleOut/group6")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/5.model")
args <- commandArgs(T)

library(caret)
library(reshape2)
library(corrplot)
library(dplyr)

seed <- 123
state <- "all_stage_sample" # "PDAC","all_stage_sample","CRC","LUAD"
ct <- args[1] # tumor
level <- args[2] # "pathway","module"
meta <- args[3] # "status","normal","tumor"
m <- args[4]
# ct <- "all"
# level <- "module"
# meta <- "status"
# m <- "2"
# method <- c('treebag','bagEarth','bagEarthGCV','C5.0','rpart','rpart1SE','rpart2',
            # 'rpartScore','fda','avNNet','earth','gcvEarth','pcaNNet','pda',
            # 'multinom','ranger','C5.0Rules','C5.0Tree','evtree','AdaBag','parRF',
            # 'rFerns','ordinalRF','rf','rfRules','RRF','RRFglobal','wsrf','naive_bayes',
            # 'LogitBoost','bstTree','sparseLDA','sdwd','kernelpls','pls','glmnet',
            # 'svmLinear3','nnet','RFlda','gpls','mlp','mlpWeightDecay',
            # 'snn','gaussprRadial','svmRadialWeights')
method <- c("parRF","treebag","ranger","rf","gaussprRadial","RRFglobal")

scFEA_raw<-read_in_FC(level,meta)
if(ct=="tumor"){
  scFEA <- scFEA_raw[scFEA_raw$celltype %in% "Malignant cells",c(2,3,5)]
# }else{ # 无法处理feature中的NA
  # tem <- scFEA_raw[,c(2,5)]
  # tem$feature <- paste(scFEA_raw[,3],scFEA_raw[,4],sep="-")
  # scFEA <- tem[,c(1,3,2)]
}
colnames(scFEA) <- c("sample","feature","log2FC")

if(meta=="status"){
  l<-read.table("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/status.meta")
}

tumStage <- as.factor(l[,2])
tumSample <- l[,1]

scFEA_stage <- scFEA[scFEA$sample %in% tumSample,]
dataStage <- acast(scFEA_stage,sample~feature)
# dataStage <- dataStage[!grepl("OV",rownames(dataStage)),]
######
# 特征选择
# # 删除常数自变量，或者是方差极小的自变量
# zerovar=nearZeroVar(dataStage)
# newdata1=dataStage[,-zerovar]
# 删除与其它自变量有很强相关性的变量
# descrCorr = cor(dataStage)
# highCorr = findCorrelation(descrCorr, 0.90)
# newdata2 = dataStage[, -highCorr]
# # 删除多重共线性自变量
# comboInfo = findLinearCombos(newdata2)
# if(!is.null(comboInfo$remove)){
  # newdata3=newdata2[, -comboInfo$remove]
  # dataStage=newdata3
# }else{
  # dataStage=newdata2
# }

# 然后定义控制参数，functions是确定用什么样的模型进行自变量排序，本例选择的模型是随机森林即rfFuncs，
# 可以选择的还有lmFuncs（线性回归），nbFuncs（朴素贝叶斯），treebagFuncs（装袋决策树），caretFuncs（自定义的训练模型）。lmFuncs rfFuncs nbFuncs treebagFuncs
# method是确定用什么样的抽样方法，本例使用cv即交叉检验, 还有提升boot以及留一交叉检验LOOCV
if(level=="module"){
  subsets = seq(5,100,by=5)
}
if(level=="pathway"){
  subsets = seq(3,30,by=3)
}
ctrl1= rfeControl(functions = rfFuncs, method = "LOOCV",verbose = FALSE, returnResamp = "final") # RRF,rfRules
Profile1 = rfe(dataStage, tumStage, sizes = subsets, rfeControl = ctrl1)
ctrl4= rfeControl(functions = treebagFuncs, method = "LOOCV",verbose = FALSE, returnResamp = "final") # RRF,rfRules
Profile4 = rfe(dataStage, tumStage, sizes = subsets, rfeControl = ctrl4)

print(Profile1)
print(Profile4)

######
# 生成模型
# featurePlot(trainx,trainy,plot='box')  # y需要是factor
print(paste0(m,":",method[as.numeric(m)]))

fitControl = trainControl(method = 'LOOCV')

# 3种特征选择结果的建模结果
l1 <- list()
l2 <- list()
l3 <- list()

pro <- list() # 选择的特征
train <- list() # 训练集
for(i in 1:3){
	if(i==1){
		matr <- dataStage[,Profile1$optVariables]
		pro[[i]] <- Profile1$optVariables
	}
	if(i==2){
		matr <- dataStage[,Profile4$optVariables]
		pro[[i]] <- Profile4$optVariables
	}
	if(i==3){
		matr <- dataStage
		pro[[i]] <- colnames(dataStage)
	}
	tem_list <- list()
	for(a in 1:100){
		# set.seed(seed)
		print(a)
	# 数据划分，分成75%的训练样本和25%检验样本
		inTrain = createDataPartition(tumStage, p = 0.75, list = FALSE)
		trainx = matr[inTrain,]
		testx = matr[-inTrain,]
		trainy = tumStage[inTrain]
		testy = tumStage[-inTrain]
		tem_list[[a]] <- inTrain
		if(i==1){
			try(l1[[a]]<-train(trainx,trainy,
				method = method[as.numeric(m)],
				trControl = fitControl))
		}
		if(i==2){
			try(l2[[a]]<-train(trainx,trainy,
				method = method[as.numeric(m)],
				trControl = fitControl))
		}
		if(i==3){
			try(l3[[a]]<-train(trainx,trainy,
				method = method[as.numeric(m)],
				trControl = fitControl))
		}
	}
	train[[i]] <- tem_list
}

save(l1,l2,l3,pro,train,file=paste0(ct,"_",level,"_Fit",m,".Rdata"))
