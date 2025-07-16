library(ggplot2)
list_tumor <- read.table("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/fluxandcelltype_high.list")
list_tumor[,2] <- "tumor"
list_normal <- read.table("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/fluxandcelltype_normal.list")
list_normal[,2] <- "tumor"
list_normal[1:33,2] <- "normal"
color <- c("Malignant cells"="#4F8AC1","CD4+ naive T cells"="#00FFFF","CD4+ effector T cells"="#87CEFA","CD4+ Tcm"="#000080","CD4+ Tem"="#4169E1","Tregs"="#7FFF00","CD8+ naive T cells"="#DAA520","CD8+ effector T cells"="#7BCD7B","CD8+ Tcm"="#218A21","CD8+ Tem"="#2F4F4F","NK cells"="#8B0000","Naive B cells"="#87481F","Activated B cells"="#FFB6C1","Memory B cells"="#7F7F7F","Plasma cells"="#FFD700","Macrophages"="#FAB17A","Monocytes"="#FF6246","Dendritic cells"="#FF1493","Mast cells"="#FF0000","Neutrophils"="#FFDEAD","Endothelial cells"="#0000FF","Epithelial cells"="#708090","Fibroblasts"="#440053","Erythrocytes"="#483D8B","Oligodendrocytes"="#D2691E","Microglial cells"="#000000","Neural progenitor cells"="#D8BFD8","Neurons"="#D2B48C","Astrocytes"="#7E5094","Acinar cells"="#800080","Progenitor cells"="#9370DB","HSCs"="#FF00FF","GMPs"="#E6E6FA","Unassigned"="#BEBEBE","Megakaryocytes"="#BC8F8F","Basophils"="#6A5ACD","Basal cells"="#F6D7BE","Smooth muscle cells"="#F660C6","Melanocytes"="#FFFACD","Myelocytes"="#F0FFFF","Neuroblasts"="#C1CDCD","Oligodendrocyte precursor cells"="#8B668B","Pericytes"="#CDB5CD","Plasmacytoid dendritic cells"="#008B8B","CLPs"="#96CDCD","CMPs"="#FAF0E6","Chondrocytes"="#53868B","Endocrine cells"="#00F5FF","Enterocytes"="#FFDAB9","Erythroid progenitor cells"="#8B3626","Hepatocytes"="#8B8378","Keratinocytes"="#E0EEE0","MEPs"="#8B1C62","Mesenchymal cells"="#836FFF","Muscle cells"="#27408B","Pancreatic ductal cells"="#6CA6CD","Platelets"="#CAE1FF","Retinal ganglion cells"="#B4CDCD","Astrocyte precursor cells"="#8B7B8B","Cajal-retzius cells"="#1C1C1C","Ciliated cells"="#FFE1FF","Club cells"="#5D478B","Eosinophils"="#4A708B","Glial cells"="#B9D3EE","Intermediated cells"="#4F4F4F","Interneurons"="#6E7B8B","Leukocytes"="#7A8B8B","Loop of Henle cells"="#838B8B","Mesothelial cells"="#FFE4E1","Mucosa cells"="#8B8386","Pinealocytes"="#CD8C95","Pneumocytes"="#8B3A62","Purkinje cells"="#EEE9BF","Pyramidal cells"="#CDC5BF","Renal principal cells"="#FFA500","Stem cells"="#CD8500","Adipocytes"="#778899","Plasmablasts"="#ADD8E6")

for(i in 1:nrow(list_tumor)){
  sam <- unlist(strsplit(list_tumor[i,1], ".", fixed = TRUE))[1]
  ct <- read.table(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/celltype/",sam,".celltype"),header = T,sep="\t")
  if(i==1){
    re <- as.data.frame(table(ct[,2]))
    re$sample <- sam
  }else{
    tem <- as.data.frame(table(ct[,2]))
    tem$sample <- sam
    re <- rbind(re,tem)
  }
}

p <- ggplot(re, aes(x=sample ,y=Freq ,fill=Var1))+
  geom_bar(position = "fill", stat = "identity", width = 0.6)+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "white" ),
        axis.text.x=element_text(colour="black",angle=90,hjust=0.5,size = 5),
        axis.text.y=element_text(colour="#222222",angle=0,hjust=1,size = 11),
        strip.text.x = element_text(size = 13, color = "#1E1E1E"),
        strip.background = element_rect(
          color="#CFD1D0", fill="white", size=1.5, linetype="solid")
  )+
  scale_fill_manual(name="Cell type", values=color)+
  labs(x="",y="Cell component")
  
pdf("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/3.metabolic_landscape/celltype/celltypeSummary.pdf",width=30,height=12)
print(p)
dev.off()

for (type in unique(list_normal$V2)) {
  list <- list_normal[which(list_normal$V2==type),]
  for(i in 1:nrow(list)){
    sam <- unlist(strsplit(list[i,1], ".", fixed = TRUE))[1]
    ct <- read.table(paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/celltype/",sam,".celltype"),header = T,sep="\t")
    if(i==1){
      re <- as.data.frame(table(ct[,2]))
      re$sample <- sam
    }else{
      tem <- as.data.frame(table(ct[,2]))
      tem$sample <- sam
      re <- rbind(re,tem)
    }
    p <- ggplot(re, aes(x=sample ,y=Freq ,fill=Var1))+
      geom_bar(position = "fill", stat = "identity", width = 0.6)+
      theme(panel.border = element_blank(),
            panel.background = element_rect(fill = "white" ),
            axis.text.x=element_text(colour="black",angle=90,hjust=0.5,size = 5),
            axis.text.y=element_text(colour="#222222",angle=0,hjust=1,size = 11),
            strip.text.x = element_text(size = 13, color = "#1E1E1E"),
            strip.background = element_rect(
              color="#CFD1D0", fill="white", size=1.5, linetype="solid")
      )+
      scale_fill_manual(name="Cell type", values=color)+
      labs(x="",y="Cell component")
    
    pdf(paste("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/3.metabolic_landscape/celltype/celltypeSummary",type,".pdf"),width=15,height=6)
    print(p)
    dev.off()
  }
}

