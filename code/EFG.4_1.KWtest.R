######
# 一、表示各细胞类型趋势一致的模块
# 读入数据
source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/4.module_pathway/KWtest")

level <- "module" # "pathway","sub_pathway","module"
meta <- "high" # high, v1, v2, normal
A<-read_in_FC(level, meta)
ct <- sort(unique(A$celltype)) # 细胞类型统计
# 细胞类型的样本数目
nSamCt <- c()
for (i in ct) {
  ct_sam <- A[which(A$celltype==i),]
  nSamCt <- c(nSamCt, length(unique(ct_sam$sample)))
}
names(nSamCt) <- ct
rem <- ct[which(as.numeric(nSamCt) < 5)]

a <- read.table(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/fluxandcelltype_",meta,".list"))[,1,drop=T]
n <- length(a)

# 每次处理一个样本
# 每个循环生成一个各模块的差异性p值向量
res <- c()
for (j in 1:n) {
  print(paste0(j,": ",a[j]))
  sam <- unlist(strsplit(a[j], ".", fixed = TRUE))[1]
  counts0 <- data.matrix(read.csv(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/flux/",sam,"_flux.csv"), header = T, row.names = 1))
  rownames(counts0) <- gsub('.',"-",rownames(counts0),fixed = T)
  rownames(counts0) <- gsub(':',"-",rownames(counts0),fixed = T)
  group <- read.table(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/celltype/",sam,".celltype"), header = F, sep = "\t")
  group[,1] <- gsub('.',"-",group[,1],fixed = T)
  group[,1] <- gsub(':',"-",group[,1],fixed = T)
  rownames(group) <- group[,1]
  group <- group[,-1,drop=F]
  colnames(group) <- c('celltype')
  # kw <- data.frame(counts0,group[,2])
  kw <- merge(group,counts0,by="row.names",all=F)
  rownames(kw) <- kw[,1]
  kw <- kw[,-1,drop=F]
  kw <- kw[which(!kw$celltype %in% c(rem,"Unassigned")),]
  re <- c()
  for (i in 2:(ncol(kw))) {
    kt <- kruskal.test(kw[,i]~celltype,data=kw)
    re[i] <- kt$p.value
  }
  re <- t(as.matrix(re))
  row.names(re) <- paste0(sam[1],sam[2],"-",sam[3])
  if(is.null(res)){
    res <- re
  }else{
    res <- rbind(res,re)
  }
}

res <- res[,-1]
colnames(res) <- colnames(counts0)
write.csv(res,file = "KWtest_flux.csv", quote = F)

# 读入通量预测结果
kw <- as.matrix(read.csv('KWtest_flux.csv'))
row.names(kw) <- kw[,1]
kw <- kw[,-1]
kw[is.na(kw)] <- 1

# # 对每个肿瘤类型进行研究
# tumor <- "scr"
# kw_t <- kw[grep(tumor,row.names(kw)),]
kw_t <- kw

# 阈值
p <- 0.05
# 保存每个样本高于阈值的module
l1 <- list()
# 保存每个样本低于阈值的module
l2 <- list()
for (i in 1:ncol(kw_t)) {
  l1[[i]] <- row.names(kw_t)[which(as.numeric(kw_t[,i]) > p)]
  l2[[i]] <- row.names(kw_t)[which(as.numeric(kw_t[,i]) < p)]
}
names(l1) <- colnames(kw_t)
names(l2) <- colnames(kw_t)

l <- c()
for (j in 1:length(l1)) {
  l[j] <- length(l1[[j]])
}
names(l)<-colnames(kw)

n_same <- 139
names(l)[which(l > n_same)]
l[which(l > n_same)]

# 找出不存在差异的样本
setdiff(row.names(kw),unlist(l1["M_59"]))

