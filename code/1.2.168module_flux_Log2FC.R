args=commandArgs(T)
ty <- args[1]
library(reshape2)

list<-read.table(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/fluxandcelltype_",ty,".list"),quote="",header=F)

bind<-c()
fc <- c()
for (i in 1:nrow(list)){
  print(paste0(i,": ",list[i,]))
  A<-read.table(paste("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/fluxandcelltype/",ty,"/",list[i,],sep=""),header=T,row.names=1,sep="\t")
  # me <- mean(apply(A[,-1],2,as.numeric))
  celltypes<-table(A$celltype)
  
  B <- A[,c(2:ncol(A))]
  C <- B
  if(min(C)<0){
    C <- C-min(C)
  }
  # B <- B/me
  
  # flux
  for (j in levels(as.factor(A[,1]))) {
    flux <- colMeans(C[which(A[,1] %in% j),])
    for (k in 1:length(flux)) {
      new <- cbind(names(flux)[k],unlist(strsplit(list[i,], ".", fixed = TRUE))[c(1)],
                   names(flux)[k],j,flux[k]) # 类似paste，连接数据形成单行的矩阵
      bind <- rbind(bind,new) # 将矩阵组合到一起
    }
  }
  # log2FC
  C <- B
  if(min(C)<0){
    C <- C-min(C)
  }
  aMean <- colMeans(C)
  ct<-A[,1,drop=F]
  D <- merge(ct,C,by="row.names",all=T)
  rownames(D)<-D[,1]
  D<-D[,-1]
  C <- melt(D)
  D <- aggregate(C$value,by=list(C$celltype, C$variable),FUN=mean)
  for(k in names(aMean)){
    D$log2FC[D$Group.2 %in% k] <- log2(D[which(D$Group.2 %in% k),"x"]/aMean[k])
  }
  df <- as.matrix(data.frame("module_1"=D$Group.2,"sample"=unlist(strsplit(list[i,], ".", fixed = TRUE))[c(1)],"module_2"=D$Group.2,"celltype"=D$Group.1,"log2FC"=D$log2FC))
  fc <- rbind(fc,df)
}
# 输出flux
colnames(bind) <- c("module_1","sample","module_2","celltype","flux")
bind <- bind[order(bind[,3],decreasing = F),]
rownames(bind) <- 1:nrow(bind)
write.table(bind,file=paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/input_df/samples.168modules_",ty,".flux"),
            sep="\t",quote=F,row.names = FALSE)

# 输出log2FC
fc <- fc[order(fc[,3],decreasing = F),]
rownames(fc) <- 1:nrow(fc)
write.table(fc,file=paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/input_df/samples.168modules_",ty,".log2FC"),
            sep="\t",quote=F,row.names = FALSE)