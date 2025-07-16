meta <- "high" # v1, v2, high
dat <- read.csv(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/input_df/samples.34pathways_",meta,".log2FC"),
         sep = "\t", quote = "")

ct <- sort(unique(dat$celltype)) # 细胞类型统计
# 细胞类型的样本数目
nSamCt <- c()
for (i in ct) {
  ct_sam <- dat[which(dat$celltype==i),]
  nSamCt <- c(nSamCt, length(unique(ct_sam$sample)))
}
names(nSamCt) <- ct
sort(nSamCt,decreasing = T)
