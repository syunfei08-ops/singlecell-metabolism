# 成对的数据的环状条形图
source("C:/Users/augus/Desktop/metabolism/code/(function)utils.R")
setwd("C:/Users/augus/Desktop/metabolism/")

level <- "pathway" # "pathway","sub_pathway","module"
scFEA_raw<-read_in_FC(level)
scFEA <- ct_adj(scFEA_raw)

ct <- levels(as.factor(scFEA$celltype))
pw <- levels(as.factor(scFEA$pathway1))

sam <- c(54,60,64,84,83,83,83,82,81,90,60,67,90,95,35,57,90,57,84,7,66,64)
names(sam) <- ct
re_num_high <- matrix(NA,nrow = length(ct),ncol = length(pw),dimnames = list(c(ct),c(pw)))
re_num_low <- matrix(NA,nrow = length(ct),ncol = length(pw),dimnames = list(c(ct),c(pw)))
re_per_high <- matrix(NA,nrow = length(ct),ncol = length(pw),dimnames = list(c(ct),c(pw)))
re_per_low <- matrix(NA,nrow = length(ct),ncol = length(pw),dimnames = list(c(ct),c(pw)))
for (i in ct) {
  for (j in pw) {
    tem <- scFEA[which(scFEA$pathway1 %in% j & scFEA$celltype %in% i),]
    n_sam <- sam[i]
    n_high <- length(which(tem$log2FC>0)) 
    n_low <- length(which(tem$log2FC<0)) 
    per <- n_high/n_sam
    re_num_high[i,j] <- n_high
    re_per_high[i,j] <- per
    per <- n_low/n_sam
    re_num_low[i,j] <- n_low
    re_per_low[i,j] <- per
  }
}
library(ggplot2)
library(reshape2)
melt_re_high <- melt(re_per_high)
melt_re_low <- melt(re_per_low)
melt_re <- mutate(melt_re_high, high = melt_re_high$value*100,low = -(melt_re_low$value*100))
melt_re <- data.frame(celltype=melt_re$Var1,pathway=melt_re$Var2,high_value=melt_re$high,low_value=melt_re$low)

ggplot(data = melt_re,aes(x = pathway,fill=celltype))+
  geom_col(aes(y = high_value), fill = "#5d8402", position = 'dodge', width = 0.5)+
  geom_col(aes(y = low_value), fill = "#817d79", position = 'dodge', width = 0.5)+
  geom_text(aes(y = 100, label = pathway)) +
  coord_polar()+
  theme_void()

ggplot(data = melt_re,aes(x=pathway,fill=celltype))+
  geom_col(aes(y = high_value),position = 'dodge', width = 0.5)+
  geom_col(aes(y = low_value),position = 'dodge', width = 0.5)+
  geom_text(aes(y = 100, label = pathway)) +
  geom_text(aes(y = ifelse(high_value >= 15, 8, (high_value + 10)), color = ifelse(high_value >= 15, 'white', '#5d8402'), label = round(high_value, 2)), size = 2.5)+#ifelse(test, yes, no),所以这里代表Trees如果大于等于15时，y=8，颜色为白色，如果小于15则y为Trees+10，颜色为绿色
  geom_text(aes(y = ifelse(low_value <= -15, -8, (low_value - 10)), color = ifelse(low_value <= -15, 'white', '#817d79'), label = round(low_value, 1)), size = 2.5)+#Pop.10小于等于-15时，y=-8，颜色为白色，当Pop.10大于-15时，y=Pop.10-10，颜色为绿色
  theme(legend.position = 'none')+
  coord_polar()+
  theme_void()

# Hiplot的环状条形图
# paletteer_d("wesanderson::GrandBudapest1")