library(readxl)
library(writexl)
# load( "CRG/extdata/featureselect_data.Rdata")
# exp <- round(exp,0)
# colnames(exp) <- gsub('\\.1','',colnames(exp))
# write.table(exp,'~/Desktop/tcgaexp.txt',sep = '\t',quote = F)
# group <- ifelse(as.numeric(substr(colnames(exp),14,15))<10,'Tumor','Normal')
# group <- as.data.frame(group)
# group <- cbind(colnames(exp),group )
# write_xlsx(group ,'~/Desktop/tcgaexp.xlsx')


CRGs <- read_xlsx('CRG/result/files/CRGs_209.xlsx')
A <- read_xlsx('~/Downloads/GSE48000_degs_1.2_0.05_limma.xlsx')
B <- read_xlsx('~/Downloads/GSE19151_degs_1.2_0.05_limma.xlsx')
load('CRG/extdata/DEG_TN_TCGA.Rdata')
C <- rownames(DEG_limma_voom)[DEG_limma_voom$adj.P.Val<0.05&abs(DEG_limma_voom$logFC)>log2(1)]
load('../project/R/my_luad/Rdata/DEG_Normal-Tumor.Rdata')
D <- rownames(DEG_limma_voom)[DEG_limma_voom$adj.P.Val<0.05&abs(DEG_limma_voom$logFC)>log2(1)]
hubgenes <-intersect(intersect(intersect(intersect(A$Tag,B$Tag),C),D),CRGs$gene)
length(intersect(D,CRGs$gene))
length(intersect(intersect(C,D),CRGs$gene))
length(intersect(intersect(intersect(B$Tag,C),D),CRGs$gene))
length(intersect(intersect(intersect(intersect(A$Tag,B$Tag),C),D),CRGs$gene))
########################################
library(grid)
library(futile.logger)
library(VennDiagram)
venn.plot <- draw.quintuple.venn(
  area1 = 2972,
  area2 = 9002,
  area3 = 7113,
  area4 =25822,
  area5 = 209,
  n12 = 810,
  n13 = 1444,
  n14 = 2108,
  n15 = 53,
  n23 = 1927,
  n24 = 3928,
  n25 = 50,
  n34 = 5412,
  n35 = 141,
  n45 = 157,
  n123 = 451,
  n124 = 633,
  n125 = 18,
  n134 = 1188,
  n135 = 43,
  n145 = 47,
  n234 = 1557,
  n235 = 35,
  n245 = 41,
  n345 = 116,
  n1234 = 373,
  n1235 = 16,
  n1245 = 15,
  n1345 = 37,
  n2345 = 31,
  n12345 = 13,
  category = c("GSE48000 DEGs", "GSE19151 DEGs", "TCGA DEGs","FUSCC DEGs", "CRGs"),
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  lwd = rep(2, 5), lty = rep("solid",5),
  cex = 2,
  fontface = 'bold',
  cat.cex = 2,
  cat.fontface = 'bold'
)
png(filename = "CRG/result/figs/02_Venn_diagram.png",width = 3000,height = 3000,res=200)
grid.draw(venn.plot)
dev.off()

#####################################
MSI <- read_xlsx('CRG/result/files/02_tmb_msi.xlsx',sheet = 'MSI')
TMB <- read_xlsx('CRG/result/files/02_tmb_msi.xlsx',sheet = 'TMB')
HRD <- read_xlsx('CRG/result/files/02_tmb_msi.xlsx',sheet = 'HRD')
MATH <- read_xlsx('CRG/result/files/02_tmb_msi.xlsx',sheet = 'MATH')
LOH <- read_xlsx('CRG/result/files/02_tmb_msi.xlsx',sheet = 'LOH')
ESTIMATE <- read_xlsx('CRG/result/files/02_tmb_msi.xlsx',sheet = 'ESTIMATE')
purity <- read_xlsx('CRG/result/files/02_tmb_msi.xlsx',sheet = 'purity')
ploidy <- read_xlsx('CRG/result/files/02_tmb_msi.xlsx',sheet = 'ploidy')
CIBERSORT <- read_xlsx('CRG/result/files/02_CR2_CIBERSORT.xlsx')
clin_tumor$SampleName <- clin_tumor$sampleID

hubdat$SampleName<- hubdat$sample
TMB <- merge(hubdat,TMB,by.x='SampleName')
TMB$group <- ifelse(TMB$RASGRP2>median(TMB$RASGRP2),'High RASGRP2','Low RASGRP2')
pp <- ggboxplot(TMB,x = 'group',y = 'TMB',fill='group')+theme_bw()+ylab('TMB')+xlab('')+ggtitle('TCGA')+ylim(min(TMB$TMB)-3,max(TMB$TMB)+10)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(size=1,textsize=12,comparisons = list(c('High RASGRP2','Low RASGRP2')),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/02_TMB_TCGA.tiff',pp,width = 5,height = 5,dpi=200)

hubdat$SampleName<- hubdat$sample
purity <- merge(hubdat,purity,by.x='SampleName')
purity$group <- ifelse(purity$RASGRP2>median(purity$RASGRP2),'High RASGRP2','Low RASGRP2')
pp <- ggboxplot(purity,x = 'group',y = 'purity',fill='group')+theme_bw()+ylab('purity')+xlab('')+ggtitle('TCGA')+ylim(min(purity$purity)-0.1,max(purity$purity)+0.5)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(size=1,textsize=12,comparisons = list(c('High RASGRP2','Low RASGRP2')),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/02_purity_TCGA.tiff',pp,width = 5,height = 5,dpi=200)


CIBERSORT$SampleName<- CIBERSORT$ID
CIBERSORT <- merge(hubdat,CIBERSORT,by.x='SampleName')
CIBERSORT$group <- ifelse(CIBERSORT$RASGRP2>median(CIBERSORT$RASGRP2),'High RASGRP2','Low RASGRP2')
rownames(CIBERSORT) <- CIBERSORT$ID
pd <-CIBERSORT%>%
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(CIBERSORT)[-c(1:5,ncol(CIBERSORT))]),names_to = 'Gene',values_to = 'Composition')
pd$Gene <- factor(pd$Gene)
df_p_val1 <- pd %>%
  group_by(Gene) %>%
  t_test(formula = Composition ~ group) %>%
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>%
  add_xy_position(x='Gene')

pp <- ggboxplot(pd,x = 'Gene', y='Composition',fill='group')+theme_bw()+ylab('CIBERSORT immune infiltration score')+xlab('TCGA')+ylim(min(pd$Composition),max(pd$Composition)+0.2)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  stat_pvalue_manual(size=10,df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(axis.text.x = element_text(size=18,face='bold',color='black',angle = 90,hjust = 0.8,vjust = 0.9),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/02_CIBERSORT_TCGA.tiff',pp,width = 15,height = 10,dpi=200)



ESTIMATE<- merge(hubdat,ESTIMATE,by.x='SampleName')
ESTIMATE$group <- ifelse(ESTIMATE$RASGRP2>median(ESTIMATE$RASGRP2),'High RASGRP2','Low RASGRP2')

imm <- ggboxplot(ESTIMATE,x = 'group',y = 'ImmuneScore',fill='group')+theme_bw()+ylab('ImmuneScore')+xlab('')+ggtitle('TCGA')+ylim(min(ESTIMATE$ImmuneScore)-1000,max(ESTIMATE$ImmuneScore)+1000)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(size=1,textsize=12,comparisons = list(c('High RASGRP2','Low RASGRP2')),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/02_imm_TCGA.tiff',imm,width = 5.5,height = 5,dpi=200)

strom <- ggboxplot(ESTIMATE,x = 'group',y = 'StromalScore',fill='group')+theme_bw()+ylab('StromalScore')+xlab('')+ggtitle('TCGA')+ylim(min(ESTIMATE$StromalScore)-1000,max(ESTIMATE$StromalScore)+1000)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(size=1,textsize=12,comparisons = list(c('High RASGRP2','Low RASGRP2')),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/02_strom_TCGA.tiff',strom,width = 5.5,height = 5,dpi=200)

est <- ggboxplot(ESTIMATE,x = 'group',y = 'ESTIMATEScore',fill='group')+theme_bw()+ylab('ESTIMATEScore')+xlab('')+ggtitle('TCGA')+ylim(min(ESTIMATE$ESTIMATEScore)-1000,max(ESTIMATE$ESTIMATEScore)+1000)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(size=1,textsize=12,comparisons = list(c('High RASGRP2','Low RASGRP2')),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/02_est_TCGA.tiff',est,width = 5.5,height = 5,dpi=200)




# library(pROC)
# roc1 <- roc(a$status,a$Expression)
# roc2 <- roc(a$status,a$MSI)
# png("CRG/result/figs/02_roc_compare.png",width = 1500,height = 1500,res=200)
# plot(roc1,col="red",lwd=4,cex.axis=1.5,ann = F, xlim=c(1,0), ylim=c(0,1))
# title(xlab="Specificity", ylab="Sensitivity",cex.lab=1.6,font.lab=2)
# lines(roc2$specificities,roc2$sensitivities, type="l",lwd=4,col="orange")
# legend(0.55,0.2,c(paste('CR2 =',round(roc1$auc,3)) ,
#                   paste('MSI =',round(roc2$auc,3)) ),
#        x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
#        bty='n',seg.len = 1.5,cex = 1.2)
# dev.off()
# 
# MSI$group <- ifelse(MSI$Expression>median(MSI$Expression),'High CR2','Low CR2')
# pp <- ggboxplot(MSI,x = 'group',y = 'MSI',fill='group')+theme_bw()+ylab('MSI')+xlab('')+ggtitle('TCGA')+ylim(min(MSI$MSI)-0.1,max(MSI$MSI)+0.1)+
#   scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
#   geom_signif(size=1,textsize=12,comparisons = list(c('High CR2','Low CR2')),test = wilcox.test,map_signif_level = T)+
#   theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
#         axis.text.x = element_text(size=18,face='bold',color='black'),
#         axis.text.y = element_text(size=18,face='bold',color='black'),
#         text=element_text(size=18,face='bold'),
#         legend.title = element_blank(),legend.position = "none",
#         axis.line = element_line(size=1),
#         panel.grid = element_blank(),panel.border = element_rect(size=2))
# ggsave('CRG/result/figs/02_MSI_TCGA.tiff',pp,width = 5,height = 5,dpi=200)
# 
# LOH$group <- ifelse(LOH$Expression>median(LOH$Expression),'High CR2','Low CR2')
# pp <- ggboxplot(LOH,x = 'group',y = 'LOH',fill='group')+theme_bw()+ylab('LOH')+xlab('')+ggtitle('TCGA')+ylim(min(LOH$LOH)-0.1,max(LOH$LOH)+0.1)+
#   scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
#   geom_signif(size=1,textsize=12,comparisons = list(c('High CR2','Low CR2')),test = wilcox.test,map_signif_level = T)+
#   theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
#         axis.text.x = element_text(size=18,face='bold',color='black'),
#         axis.text.y = element_text(size=18,face='bold',color='black'),
#         text=element_text(size=18,face='bold'),
#         legend.title = element_blank(),legend.position = "none",
#         axis.line = element_line(size=1),
#         panel.grid = element_blank(),panel.border = element_rect(size=2))
# ggsave('CRG/result/figs/02_LOH_TCGA.tiff',pp,width = 5,height = 5,dpi=200)
# 
# 
# ploidy$group <- ifelse(ploidy$Expression>median(ploidy$Expression),'High CR2','Low CR2')
# pp <- ggboxplot(ploidy,x = 'group',y = 'ploidy',fill='group')+theme_bw()+ylab('ploidy')+xlab('')+ggtitle('TCGA')+ylim(min(ploidy$ploidy)-0.1,max(ploidy$ploidy)+1)+
#   scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
#   geom_signif(size=1,textsize=12,comparisons = list(c('High CR2','Low CR2')),test = wilcox.test,map_signif_level = T)+
#   theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
#         axis.text.x = element_text(size=18,face='bold',color='black'),
#         axis.text.y = element_text(size=18,face='bold',color='black'),
#         text=element_text(size=18,face='bold'),
#         legend.title = element_blank(),legend.position = "none",
#         axis.line = element_line(size=1),
#         panel.grid = element_blank(),panel.border = element_rect(size=2))
# ggsave('CRG/result/figs/02_ploidy_TCGA.tiff',pp,width = 5,height = 5,dpi=200)
# 
# 
# MATH$group <- ifelse(MATH$Expression>median(MATH$Expression),'High CR2','Low CR2')
# pp <- ggboxplot(MATH,x = 'group',y = 'MATH',fill='group')+theme_bw()+ylab('MATH')+xlab('')+ggtitle('TCGA')+ylim(min(MATH$MATH)-5,max(MATH$MATH)+10)+
#   scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
#   geom_signif(size=1,textsize=12,comparisons = list(c('High CR2','Low CR2')),test = wilcox.test,map_signif_level = T)+
#   theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
#         axis.text.x = element_text(size=18,face='bold',color='black'),
#         axis.text.y = element_text(size=18,face='bold',color='black'),
#         text=element_text(size=18,face='bold'),
#         legend.title = element_blank(),legend.position = "none",
#         axis.line = element_line(size=1),
#         panel.grid = element_blank(),panel.border = element_rect(size=2))
# ggsave('CRG/result/figs/02_MATH_TCGA.tiff',pp,width = 5,height = 5,dpi=200)
# 
# 
# 
# HRD$group <- ifelse(HRD$Expression>median(HRD$Expression),'High CR2','Low CR2')
# pp <- ggboxplot(HRD,x = 'group',y = 'HRD',fill='group')+theme_bw()+ylab('HRD')+xlab('')+ggtitle('TCGA')+ylim(min(HRD$HRD)-3,max(HRD$HRD)+10)+
#   scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
#   geom_signif(size=1,textsize=12,comparisons = list(c('High CR2','Low CR2')),test = wilcox.test,map_signif_level = T)+
#   theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
#         axis.text.x = element_text(size=18,face='bold',color='black'),
#         axis.text.y = element_text(size=18,face='bold',color='black'),
#         text=element_text(size=18,face='bold'),
#         legend.title = element_blank(),legend.position = "none",
#         axis.line = element_line(size=1),
#         panel.grid = element_blank(),panel.border = element_rect(size=2))
# ggsave('CRG/result/figs/02_HRD_TCGA.tiff',pp,width = 5,height = 5,dpi=200)


#########################hub基因正常和肿瘤boxplot
hubgenes='RASGRP2'
library(tidyverse)
library(rstatix)
library(ggtext)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
library(tidyr)
load( "CRG/extdata/featureselect_data.Rdata")
group <- ifelse(as.numeric(substr(colnames(exp_CRGs),14,15))<10,'Tumor','Normal')
hubdat <- as.data.frame(t(exp_CRGs[rownames(exp_CRGs)%in%hubgenes,]))
hubdat <- do.call(data.frame,hubdat )
hubdat <- cbind(group,hubdat )
rownames(hubdat) <- colnames(exp_mnc)
hubdat$sample <- rownames(hubdat)
####TCGA汇总
# pd <-hubdat%>% 
#   rownames_to_column('Sample')%>%
#   pivot_longer(cols=c(colnames(hubdat)[-c(1,ncol(hubdat))]),names_to = 'Gene',values_to = 'Composition')
# pd$Gene <- factor(pd$Gene)
# df_p_val1 <- pd %>% 
#   group_by(Gene) %>% 
#   wilcox_test(formula = Composition ~ group) %>% 
#   add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
#   add_xy_position(x='Gene')
# dev.new()
pp <- ggboxplot(hubdat,x = 'group',y = hubgenes,fill='group')+theme_bw()+ylab('Expression')+xlab('')+ggtitle(paste0('TCGA ',hubgenes))+ylim(3,14)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(size=1,textsize=12,comparisons = list(c("Tumor", "Normal")),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/02_hub_TCGA.tiff',pp,width = 5,height = 5,dpi=200)

############TCGA 甲基化
meth <- read.table('LUADdata/TCGA/TCGA.LUAD.sampleMap_HumanMethylation450.txt',header = T,row.names = 1)
colnames(meth) <- gsub('\\.','-',colnames(meth) )
probe <- read.table('LUADdata/TCGA/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy',header = T,sep = '\t')
probe <- probe[probe$gene!='.',]
probe <- probe[probe$id%in%rownames(meth),]
probe <- probe[probe$gene%in%hubgenes,]
hub_meth <- data.frame()
i=0
for(gene in hubgenes){
  i=i+1
  hub_meth <- rbind(hub_meth,colSums(meth[rownames(meth)%in%probe$id[probe$gene==gene],]))
  rownames(hub_meth)[i] <- gene
  colnames(hub_meth) <- colnames(meth)
}
hub_meth <- as.data.frame(t(hub_meth))
hub_meth <- na.omit(hub_meth)
group <- ifelse(as.numeric(substr(rownames(hub_meth),14,15))<10,'Tumor','Normal')
# Normal  Tumor 
# 30    435 
hub_meth1 <- do.call(data.frame,hub_meth )
hub_meth1$sample <- rownames(hub_meth)
hub_meth1 <- cbind(group,hub_meth1 )
rownames(hub_meth1) <- hub_meth1$sample
# pd <-hub_meth1 %>% 
#   rownames_to_column('Sample')%>%
#   pivot_longer(cols=c(colnames(hub_meth1)[-c(1,ncol(hub_meth1))]),names_to = 'Gene',values_to = 'Composition')
# pd$Gene <- factor(pd$Gene)
# df_p_val1 <- pd %>% 
#   group_by(Gene) %>% 
#   wilcox_test(formula = Composition ~ group) %>% 
#   add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
#   add_xy_position(x='Gene')
pp <- ggboxplot(hub_meth1 ,x = 'group',y = hubgenes,fill='group')+theme_bw()+ylab('Methylation')+xlab('')+ggtitle(paste0('TCGA ',hubgenes))+ylim(3,18)+
  scale_fill_manual(values=c("orange", "darkgreen"))+
  geom_signif(size=1,textsize=10,comparisons = list(c("Tumor", "Normal")),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))

ggsave('CRG/result/figs/02_hubMETH_TCGA.tiff',pp,width = 5,height =5,dpi=200)

#############TCGA cox
load( "CRG/extdata/featureselect_data.Rdata")
load('LUADdata/TCGA/clinic.Rdata')
clin_tcga <- CLINIC[,c("sampleID","gender","age_at_initial_pathologic_diagnosis","pathologic_stage","pathologic_T",
                       "pathologic_N","pathologic_M","vital_status","tobacco_smoking_history_indicator","new_neoplasm_event_type","radiation_therapy",
                       "additional_pharmaceutical_therapy",'history_of_neoadjuvant_treatment','targeted_molecular_therapy','days_to_additional_surgery_locoregional_procedure',
                       'days_to_additional_surgery_metastatic_procedure')]
colnames(clin_tcga) <- c("sampleID",'Gender','Age','Stage',"Tumor","Lymph_Node","Metastasis",'Status','Smoking','New_tumor_type',"radiation_therapy",
                         'pharmaceutical_therapy', 'neoadjuvant_treatment','targeted_molecular_therapy','additional_surgery_locoregional_days','additional_surgery_metastatic_days')
clin_tcga <- as.data.frame(clin_tcga)
rownames(clin_tcga) <- clin_tcga$sampleID
clin_tumor<- clin_tcga[clin_tcga$sampleID%in%exp_CRGs_tumor$sampleID,]
clin_tumor$Gender <- ifelse(clin_tumor$Gender=='MALE','Male','Female')
clin_tumor$Stage <- ifelse(clin_tumor$Stage%in%c('Stage I','Stage IA','Stage IB'),'I',
                           ifelse(clin_tumor$Stage%in%c('Stage II','Stage IIA','Stage IIB'),"II",
                                  ifelse(clin_tumor$Stage%in%c('Stage IIIA','Stage IIIB'),"III",
                                         ifelse(clin_tumor$Stage%in%c('Stage IV'),'IV','UNKOWN'))))
clin_tumor$Tumor <- ifelse(clin_tumor$Tumor%in%c('T1', 'T1a', 'T1b'),'T1',ifelse(clin_tumor$Tumor%in%c('T2', 'T2a', 'T2b'),'T2',"T3-4"))
clin_tumor$`Lymph_Node` <- ifelse(clin_tumor$`Lymph_Node`%in%c('N0'),'N0',ifelse(clin_tumor$`Lymph_Node`%in%c('N1'),'N1',"N2-3"))
clin_tumor$Metastasis <- ifelse(clin_tumor$Metastasis%in%c('M0'),'M0','M1')
clin_tumor$Smoking <- ifelse(clin_tumor$Smoking%in%c(''),'UNKOWN',ifelse(clin_tumor$Smoking%in%c('Lifelong Non-smoker'),'Never smoker','Current smoker'))
clin_tumor$Status <- ifelse(clin_tumor$Status%in%c('LIVING'),'Alive','Dead')
clin_tumor$Age <- ifelse(clin_tumor$Age>70,'old','young')
clin_tumor <- merge(clin_tumor,clin_rna,by.x='sampleID')
clin_tumor <- clin_tumor[clin_tumor$time!=0,]
# clin_tumor <- clin_tumor[,-c(10:12)]
clin_tumor$RFS <- clin_tumor$additional_surgery_locoregional_days
clin_tumor$RFS[which(is.na(clin_tumor$RFS))] <- clin_tumor$days_to_death[which(is.na(clin_tumor$RFS))]
clin_tumor$RFS[which(is.na(clin_tumor$RFS))] <- clin_tumor$days_to_last_followup[which(is.na(clin_tumor$RFS))]
clin_tumor$RFS <- round(clin_tumor$RFS/30,3)
clin_tumor$RFS.E <- clin_tumor$status
clin_tumor$RFS.E[grep('Recurrence',clin_tumor$New_tumor_type)] <- 1

hubdat <-  exp_CRGs_tumor[,c('sampleID',hubgenes)]
modeldat<- merge(hubdat,clin_tumor,by.x = 'sampleID')
lxdata <- modeldat
lxdata1 <- modeldat
lxdata1$Age <- ifelse(lxdata$Age=='young',1,2 )
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$Stage <- ifelse(lxdata$Stage=='I',1,ifelse(lxdata$Stage=='II',2 ,ifelse(lxdata$Stage=='III',3,ifelse(lxdata$Stage=='IV',4,0))))
lxdata1$Stage <- factor(lxdata1$Stage,levels=c(0,1,2,3,4),labels=c('UNKOWN','I','II','III','IV'))

library(survival)
library(survminer)
library(rms)
library(survivalROC)
library(dplyr)
covariates <- c(hubgenes,'Age','Stage','Tumor','Lymph_Node','Metastasis')
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = lxdata1)})
#提取HR，95%置信区间和p值
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         #获取HR
                         HR <-signif(x$coef[2], digits=2);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })

res <- t(as.data.frame(univ_results, check.names = FALSE))
res <-as.data.frame(res,stringsAsFactors=F)
HR=gsub("[\\(\\)]","",res$`HR (95% CI for HR)`)
HR=gsub("-"," ",HR)
HR=as.data.frame(do.call(cbind,strsplit(HR," ")),stringsAsFactors=F)
names(HR)=rownames(res)
HR <- t(HR)
mulcox<- coxph(Surv(time, status) ~ RASGRP2+Age+Tumor+Lymph_Node+Metastasis, data = lxdata);mulcox


forest <- ggforest(mulcox,  #coxph得到的Cox回归结果
                   data = lxdata1,  #数据集
                   main = 'Hazard ratio of all data',  #标题
                   cpositions = c(0.05, 0.15, 0.35),  #前三列数值的距离
                   fontsize = 1.4, #字体大小
                   refLabel = 'Reference', #相对变量的标签
                   noDigits = 3 #HR、95%CI的小数位
)
ggsave(filename="CRG/result/figs/02_RASGRP2muti_forest.tiff",forest, width=17, height=10)

lxdata1$Group <- ifelse(modeldat$RASGRP2>median(modeldat$RASGRP2),'High RASGRP2_expr','Low RASGRP2_expr')
png("CRG/result/figs/02_OSKM_RASGRP2_expr_TCGA.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(time, status) ~ Group, data=lxdata1),
                surv.median.line = "hv",size=3,
                title = "TCGA RASGRP2",
                legend.labs = c('High RASGRP2','Low RASGRP2'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("CRG/result/figs/02_RFSKM_RASGRP2_expr_TCGA.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(RFS, RFS.E) ~ Group, data=lxdata1),size=3,
                surv.median.line = "hv",
                title = "TCGA RASGRP2",
                legend.labs = c('High RASGRP2','Low RASGRP2'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("red", "blue",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of RFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()




# Female   Male 
#  270      232
# 53.78%   46.22%
pp <- ggboxplot(lxdata1,x = 'Gender',y = 'RASGRP2',fill='Gender')+theme_bw()+ylab('Expression')+xlab('')+ggtitle('TCGA RASGRP2')+ylim(min(lxdata1$RASGRP2)-1,max(lxdata1$RASGRP2)+2)+
  scale_fill_manual(values=c("#354e97", "#df5b3f"))+
  geom_signif(size=1,textsize=12,comparisons = list(c("Male", "Female")),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/02_gender_TCGA.tiff',pp,width = 5,height = 5,dpi=200)

#total 502
# UNKOWN      I     II    III     IV 
#      8    270    119     80     25 
#1.59%    53.78% 23.71% 15.94%  4.98%
pp <- ggboxplot(lxdata1[lxdata1$Stage!='UNKOWN',],x = 'Stage',y = 'RASGRP2',fill='Stage')+theme_bw()+ylab('Expression')+xlab('')+ggtitle('TCGA RASGRP2')+ylim(min(lxdata1$RASGRP2)-1,max(lxdata1$RASGRP2)+16)+
  scale_fill_manual(values=c("#354e97","70a3c4",'#f5b46f', "#df5b3f"))+
  geom_signif(size=1,textsize=10,comparisons = list(c("I", "II"),c("II", "III"),c("III",'IV'),c("II", "IV"),c("I", "III"),c("I", "IV")),y_position = c(13,14,15,18,21,24,27),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/02_Stage_TCGA.tiff',pp,width = 5,height = 5,dpi=200)
# young   old 
#   330   162
#67.07% 32.93%
pp <- ggboxplot(lxdata1[lxdata1$Age%in%c("young", "old"),],x = 'Age',y = 'RASGRP2',fill='Age')+theme_bw()+ylab('Expression')+xlab('')+ggtitle('TCGA RASGRP2')+ylim(min(lxdata1$RASGRP2)-1,max(lxdata1$RASGRP2)+3)+
  scale_fill_manual(values=c("#354e97","70a3c4"))+
  geom_signif(size=1,textsize=12,comparisons = list(c("young", "old")),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/02_Age_TCGA.tiff',pp,width = 5,height = 5,dpi=200)



###########################################
load("LUADdata/GEO/GSE19151.Rdata")
group <- clindata133[colnames(exp133),4]

load("LUADdata/GEO/GSE48000.Rdata")
group <- clindata133[colnames(exp133),2]

exp <- as.data.frame(t(exp133[hubgenes,]))
exp1 <- do.call(data.frame,exp)
exp1$sample <- rownames(exp)
exp1 <- cbind(group,exp1 )
rownames(exp1) <- exp1$sample
# pd <-exp1 %>% 
#   rownames_to_column('Sample')%>%
#   pivot_longer(cols=c(colnames(exp1)[-c(1,ncol(exp1))]),names_to = 'Gene',values_to = 'Composition')
# pd$Gene <- factor(pd$Gene)
# df_p_val1 <- pd %>% 
#   group_by(Gene) %>% 
#   wilcox_test(formula = Composition ~ group) %>% 
#   add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
#   add_xy_position(x='Gene')

pp <- ggboxplot(exp1,x = 'group',y = hubgenes,fill='group')+theme_bw()+ylab('Expression')+xlab('')+ggtitle(paste0('GSE19151 ',hubgenes))+ylim(9,12)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(size=1,textsize=12,comparisons = list(c("VTE", "Control")),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))

pp <- ggboxplot(exp1,x = 'group',y = hubgenes,fill='group')+theme_bw()+ylab('Expression')+xlab('')+ggtitle(paste0('GSE48000 ',hubgenes))+ylim(10,14)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(size=1,textsize=12,comparisons = list(c("VTE", "Control")),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))

ggsave('CRG/result/figs/02_VTEexp_GSE19151.tiff',pp,width = 5,height = 5,dpi=200)
ggsave('CRG/result/figs/02_VTEexp_GSE48000.tiff',pp,width = 5,height = 5,dpi=200)

#########################################fuscc
rna_data<-read.table('LUADdata/FUSCC/all.counts.id1.txt',row.names = 1,header = T)
rna_data<-rna_data[,c(6:409)]
names<-c("X1684LC","X1890N","X2310N","X3238LC", "X3244LC","X3246LC","X3255LC","X3284LC","X3338LC" 
         ,"X3439LC","X4873LC" ,"X4894LC", "X4912LC","X4933LC","X4941LC","X4942LC"  
         ,"X4948LC","X5069LC","X5081LC","X5091LC")
colnames(rna_data)<-c(sub("X","",colnames(rna_data)[colnames(rna_data)%in%names]),colnames(rna_data)[which(!colnames(rna_data)%in%names)])

load('LUADdata/FUSCC/clin_data.Rdata')

Barcode.T<-clin_data$RNA.LC
Barcode.N<-clin_data$RNA.N
Barcodes<-c(as.character(Barcode.T),as.character(Barcode.N))
Data.T<-rna_data[,which(colnames(rna_data)%in%Barcode.T)]
Data.N<-rna_data[,which(colnames(rna_data)%in%Barcode.N)]
rna_data1<-rna_data[,colnames(rna_data)%in%Barcodes]
data<-rna_data1[,colnames(rna_data1)[which(colnames(rna_data1)%in%Barcodes)]]
group_list<-ifelse(colnames(data)%in%as.character(Barcode.T),"Tumor","Normal")
hubdat <- as.data.frame(t(data[hubgenes,]))
hubdat <- do.call(data.frame,hubdat )
hubdat <- log2(hubdat+1)
hubdat <- cbind(group_list,hubdat )
rownames(hubdat) <- colnames(data)
hubdat$sample <- rownames(hubdat)

# pd <-hubdat%>% 
#   rownames_to_column('Sample')%>%
#   pivot_longer(cols=c(colnames(hubdat)[-c(1,ncol(hubdat))]),names_to = 'Gene',values_to = 'Composition')
# pd$Gene <- factor(pd$Gene)
# df_p_val1 <- pd %>% 
#   group_by(Gene) %>% 
#   wilcox_test(formula = Composition ~ group_list) %>% 
#   add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
#   add_xy_position(x='Gene')

pp <- ggboxplot(hubdat,x = 'group_list',y = hubgenes,fill='group_list')+theme_bw()+ylab('Expression')+xlab('')+ggtitle(paste0('FUSCC ',hubgenes))+ylim(6,15)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(size=1,textsize=12,comparisons = list(c("Tumor", "Normal")),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/02_hub_FUSCC.tiff',pp,width = 5,height = 5,dpi=200)

modeldat <- as.data.frame(t(data[hubgenes,]))
modeldat<- do.call(data.frame,modeldat )
modeldat<- log2(modeldat+1)
rownames(modeldat) <- colnames(data)
modeldat$RNA.LC<- rownames(modeldat)
modeldat <- merge(modeldat,clin_data,by.y='RNA.LC')
modeldat <-  modeldat[!is.na(modeldat$RFS.E),]
modeldat$group<- ifelse(modeldat$RASGRP2>median(modeldat$RASGRP2),'High RASGRP2_expr','Low RASGRP2_expr')
png("CRG/result/figs/02_OSKM_RASGRP2_expr_FUSCC.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(OS, OS.E) ~ group, data=modeldat),
                surv.median.line = "hv",size=3,
                title = "FUSCC RASGRP2",
                legend.labs = c('High RASGRP2','Low RASGRP2'),legend=c(0.8,0.3),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("CRG/result/figs/02_RFSKM_RASGRP2_expr_FUSCC.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(RFS, RFS.E) ~ group, data=modeldat),size=3,
                surv.median.line = "hv",
                title = "FUSCC RASGRP2",
                legend.labs = c('High RASGRP2','Low RASGRP2'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("red", "blue",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of RFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("CRG/result/figs/02_RFSKM_RASGRP2_expr_later_FUSCC.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(RFS, RFS.E) ~ group, data=modeldat[modeldat$Pathologic_stage=='IIIA',]),size=3,
                surv.median.line = "hv",
                title = "FUSCC RASGRP2 Stage III",
                legend.labs = c('High RASGRP2','Low RASGRP2'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("red", "blue",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of RFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

covariates <- c(hubgenes)
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS, OS.E)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = modeldat)})
#提取HR，95%置信区间和p值
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         #获取HR
                         HR <-signif(x$coef[2], digits=2);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })


res <- t(as.data.frame(univ_results, check.names = FALSE))
res <-as.data.frame(res,stringsAsFactors=F)
HR=gsub("[\\(\\)]","",res$`HR (95% CI for HR)`)
HR=gsub("-"," ",HR)
HR=as.data.frame(do.call(cbind,strsplit(HR," ")),stringsAsFactors=F)
names(HR)=rownames(res)
HR <- t(HR)

