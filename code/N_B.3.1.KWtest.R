######
# 一、表示各细胞类型趋势一致的模块
# 读入数据
source("/xtdisk/xiaojf_group/shangyf/a_commend/R_commend/3.TumorMetabolism/code/(function)utils.R")
setwd("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/2-3.normal/KWtest")

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

sam <- meta_data$sample[which(meta_data$class=="Normal")]
n_scFEA <- scFEA[scFEA$sample %in% sam,]
meta <- meta_data[meta_data$class %in% "Normal",]

for (i in tum){
  samT <- sam[grep(i, sam)]
  # 每次处理一个样本
  # 每个循环生成一个各模块的差异性p值向量
  res <- c()
  for (a in samT) {
    counts0 <- data.matrix(read.csv(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/flux/",a,"_flux.csv"), header = T, row.names = 1))
    rownames(counts0) <- gsub('.',"-",rownames(counts0),fixed = T)
    rownames(counts0) <- gsub(':',"-",rownames(counts0),fixed = T)
    group <- read.table(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/celltype/",a,".celltype"), header = F, sep = "\t")
    group[,1] <- gsub('.',"-",group[,1],fixed = T)
    group[,1] <- gsub(':',"-",group[,1],fixed = T)
    rownames(group) <- group[,1]
    group <- group[,-1,drop=F]
    colnames(group) <- c('celltype')
    # kw <- data.frame(counts0,group[,2])
    kw <- merge(group,counts0,by="row.names",all=F)
    rownames(kw) <- kw[,1]
    kw <- kw[,-1,drop=F]
    kw <- kw[which(kw$celltype %in% c(ct)),]
    re <- c()
    for (k in 2:(ncol(kw))) {
      kt <- kruskal.test(kw[,k]~celltype,data=kw)
      re[k] <- kt$p.value
    }
    re <- t(as.matrix(re))
    row.names(re) <- paste0(a)
    if(is.null(res)){
      res <- re
    }else{
      res <- rbind(res,re)
    }
  }

  res <- res[,-1,drop=F]
  colnames(res) <- colnames(counts0)
  write.csv(res,file = paste0("normal_",i,"_KWtest_flux.csv"), quote = F)
}

## 结果读取
# 读入通量预测结果
kw <- as.matrix(read.csv(paste0("normal_",t,"_KWtest_flux.csv"))) # i in group
row.names(kw) <- kw[,1]
kw <- kw[,-1]
kw[is.na(kw)] <- 1

# 阈值
p <- 0.05
# 保存每个样本高于阈值的module
l1 <- list()
# 保存每个样本低于阈值的module
l2 <- list()
for (i in 1:ncol(kw)) {
  l1[[i]] <- row.names(kw)[which(as.numeric(kw[,i]) > p)]
  l2[[i]] <- row.names(kw)[which(as.numeric(kw[,i]) < p)]
}
names(l1) <- colnames(kw)
names(l2) <- colnames(kw)

l <- c()
for (j in 1:length(l1)) {
  l[j] <- length(l1[[j]])
}
names(l)<-colnames(kw)

n_same <- ceiling(length(samT <- sam[grep(t, sam)])/2)
paste0(t,": ",length(samT),"-",n_same)
names(l)[which(l > n_same)]
sort(l[which(l > n_same)],decreasing=T)

# 找出不存在差异的样本
setdiff(row.names(kw),unlist(l1["M_59"]))

