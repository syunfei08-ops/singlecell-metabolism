args=commandArgs(T)
sam <- args[1]
ds <- args[2]
ty <- args[3]
path <- args[4]

A <- read.csv(paste0(path,"/",sam,"_flux.csv"))
rownames(A) <- gsub('.',"-",A[,1],fixed = T)
A <- A[,-c(1)]
B<-read.table(paste0("/p300s/xiaojf_group/bigzengjy/cancerscem/",ds,"/seurat/all.celltype.forcellphonedb.txt"),header = T,sep="\t")
rownames(B) <- B[,1]
C <- merge(B,A,by="row.names",all=T)
colnames(C)[1] <- ""
colnames(C)[3] <- "celltype"
C <- C[,-c(2)]

write.table(C,file = paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/fluxandcelltype/",ty,"/",sam,".fluxandcelltype.edit.2"),
              sep = "\t", quote = F, row.names = FALSE)