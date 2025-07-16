# 各细胞类型的总体代谢能力
load("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/2-2.tumor/meta_data.Rdata")
library(ggplot2)
library(magrittr)
library(dplyr)
colorTum=c("ATC"="#1F78B4FF","BCC"="#A6CEE3FF","BRCA"="#2FE4C0","ccRCC"="#BD26E3","CRC"="#CAB2D6FF",
         "DCIS"="#738BED","ESCC"="#F45B5BFF","HNSCC"="#FB9A99FF","LUAD"="#FF7F00FF",
         "LUSC"="#19BE4A","MCC"="#B2DF8AFF","MIUBC"="#18BC9CFF","NB"="#8CBE19","NBL"="#BE9D19",
         "NSCLC"="#CCBE93FF","PDAC"="#BE4A19","PRAD"="#F4E5A9","pRCC"="#E6C43D","STAD"="#93A2CC",
         "TNBC"="#93BECC","UCEC"="#9C8749","other"="#FFFFFF")

meta_data_count <- meta_data %>% 
  group_by(tumor, group2) %>%
  summarise(count = n())


p <- ggplot(meta_data_count, aes(x=group2, y=count,fill=tumor))+
  geom_col()+
  scale_fill_manual(name="Tumor", values=colorTum)+scale_color_manual(values=colorTum)+
  ylab("Group")+xlab("Count")+
  theme_classic()+
  theme(# legend.position="none",
        # panel.grid.major=element_line(color='gray', linetype='dotted')
        panel.grid.major=element_blank()
  )+
  geom_text(aes(label = count), position = position_stack(0.5), vjust = 0)

pdf(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/2-2.tumor/MalignantandMacrophages_sample.pdf"),width=12,height=4)
print(p)
dev.off()
