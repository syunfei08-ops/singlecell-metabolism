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
  A$hsa00010<-(A$M_1+A$M_2+A$M_3+A$M_4+A$M_5)/5;
  A$hsa00620<-A$M_6
  A$hsa00020<-(A$M_5+A$M_7+A$M_8+A$M_9+A$M_10+A$M_11+A$M_12+A$M_13+A$M_14)/9;
  A$hsa00260<-(A$M_15+A$M_16+A$M_17+A$M_18+A$M_19+A$M_20+A$M_21+A$M_22+A$M_23+A$M_32)/10;
  A$hsa00220<-(A$M_65+A$M_66)/2;
  A$hsa00330<-(A$M_61+A$M_62+A$M_63+A$M_64+A$M_67+A$M_68+A$M_69)/7;
  A$hsa00480<-(A$M_23+A$M_24+A$M_25+A$M_26+A$M_69)/5;
  A$hsa00270<-(A$M_27+A$M_28+A$M_29+A$M_30+A$M_31+A$M_70)/6;
  A$hsa00030<-A$M_33
  A$hsa00061<-A$M_34
  A$hsa00062<-A$M_34
  A$hsa00071<-A$M_35
  A$hsa00250<-(A$M_37+A$M_38+A$M_39+A$M_40+A$M_48+A$M_49+A$M_50+A$M_51)/8;
  A$hsa00410<-(A$M_41+A$M_42+A$M_43+A$M_44)/4;
  A$hsa00640<-(A$M_46+A$M_47)/2;
  A$hsa00340<-A$M_52
  A$hsa00280<-(A$M_53+A$M_54+A$M_55+A$M_56+A$M_164+A$M_166)/6;
  A$hsa00350<-(A$M_57+A$M_58)/2;
  A$hsa00360<-A$M_59
  A$hsa00310<-A$M_60
  A$TCDB<-(A$M_71+A$M_72+A$M_73+A$M_74+A$M_75+A$M_76+A$M_77+A$M_78+A$M_79+A$M_80+A$M_81+A$M_83+A$M_84+A$M_85+A$M_87+A$M_88+A$M_89+A$M_90+A$M_91+A$M_92+A$M_93+A$M_94+A$M_95+A$M_96+A$M_97+A$M_98+A$M_99+A$M_100+A$M_101+A$M_102+A$M_103+A$M_105)/32;
  A$hsa00520<-(A$M_106+A$M_107+A$M_108+A$M_109+A$M_110+A$M_112)/6;
  A$hsa04922<-A$M_111
  A$hsa00900<-(A$M_113+A$M_114+A$M_115)/3;
  A$hsa00130<-A$M_116
  A$hsa00510<-(A$M_117+A$M_118+A$M_119+A$M_120+A$M_121+A$M_122+A$M_123+A$M_124)/8;
  A$hsa00515<-(A$M_125+A$M_126+A$M_127+A$M_128)/4;
  A$hsa00532<-(A$M_129+A$M_130+A$M_131)/3;
  A$hsa00534<-A$M_132
  A$hsa00230<-(A$M_133+A$M_134+A$M_135+A$M_136+A$M_137+A$M_138+A$M_139+A$M_140+A$M_141+A$M_142+A$M_143+A$M_144+A$M_145+A$M_146+A$M_147+A$M_148+A$M_149+A$M_170)/18;
  A$hsa00240<-(A$M_36+A$M_150+A$M_151+A$M_152+A$M_153+A$M_154+A$M_155+A$M_156+A$M_157+A$M_158+A$M_159+A$M_160+A$M_161+A$M_162+A$M_163+A$M_164+A$M_165+A$M_166+A$M_171+A$M_45)/20;
  A$hsa00100<-A$M_167
  A$hsa00120<-A$M_168
  A$hsa00140<-A$M_169
  
  # B<-A[,c(1,170:203)]
  # file<-paste("/asnas/pmod/zengjy/cancer-single-cell-database/scFEA/log2FC_pathways_edit/PCA/",list[i,],".34pathways",sep="")
  # write.table(B,file=file,sep="\t",quote=F)
  
  B <- A[,c(170:ncol(A))]
  # B <- B/me
  # rbind(A[,1],B)
  for (j in levels(as.factor(A[,1]))) {
    flux <- colMeans(B[which(A[,1] %in% j),])
    for (k in 1:length(flux)) {
      new <- cbind(names(flux)[k],unlist(strsplit(list[i,], ".", fixed = TRUE))[c(1)],
                   paste0("pathway",k),j,flux[k])
      bind <- rbind(bind,new)
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
  df <- as.matrix(data.frame("pathway_1"=D$Group.2,"sample"=unlist(strsplit(list[i,], ".", fixed = TRUE))[c(1)],"pathway_2"=D$Group.2,"celltype"=D$Group.1,"log2FC"=D$log2FC)) 
  fc <- rbind(fc,df)
}
# 输出flux
colnames(bind) <- c("pathway_1","sample","pathway_2","celltype","flux")
bind <- bind[order(bind[,3],decreasing = F),]
rownames(bind) <- 1:nrow(bind)
write.table(bind,file=paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/input_df/samples.34pathways_",ty,".flux"),
            sep="\t",quote=F,row.names = FALSE)

# 输出log2FC
fc <- fc[order(fc[,3],decreasing = F),]
rownames(fc) <- 1:nrow(fc)
write.table(fc,file=paste0("/xtdisk/xiaojf_group/shangyf/a_result/3.metabolism/1.input/input_df/samples.34pathways_",ty,".log2FC"),
            sep="\t",quote=F,row.names = FALSE)