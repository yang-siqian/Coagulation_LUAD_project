library(dplyr)
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

bsrgenes <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7_PLA2G4B','CR2')
JIN_genes <- c('SERPINA1','CFHR3','PPP1CB','P2RX1','PLCB3','PLCB4','PIK3R6','GP6','PIK3R1','GP1BA','PLA2G4F')
LI_genes <- c('CASP6','NLRP7','NOD1','NLRP1','NLRP2')
GENG_genes <- c('NLRP7','CASP1','PLCG1','NLRP1')
exp_tumor <- exp[,colnames(exp)%in%clin_rna$sampleID[clin_rna$group=='Tumor']]
exp_BJ_tumor <- exp_tumor[rownames(exp_tumor)%in%c(bsrgenes[-10],gsub('_','-',bsrgenes[10]),JIN_genes,LI_genes,GENG_genes),]%>%t%>%as.data.frame()
exp_BJ_tumor$sampleID <- rownames(exp_BJ_tumor)
exp_BJ_tumor <- do.call(data.frame,exp_BJ_tumor)
colnames(exp_BJ_tumor)[23] <- 'JMJD7_PLA2G4B'
lxdata <- merge(exp_BJ_tumor,clin_tumor,by.x = 'sampleID')
lxdata1 <- merge(exp_BJ_tumor,clin_tumor,by.x = 'sampleID')
lxdata1$time <- round(lxdata1$time,3)
lxdata1$Age <- ifelse(lxdata$Age=='young',1,2 )
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$Stage <- ifelse(lxdata$Stage=='I',1,ifelse(lxdata$Stage=='II',2 ,ifelse(lxdata$Stage=='III',3,ifelse(lxdata$Stage=='IV',4,0))))
lxdata1$Stage <- factor(lxdata1$Stage,levels=c(0,1,2,3,4),labels=c('UNKOWN','I','II','III','IV'))
lxdata1$Risk=0.09*lxdata1$F2-0.299*lxdata1$P2RX1-0.106*lxdata1$GUCY1A1+0.107*lxdata1$PLAUR-0.124*lxdata1$PRKCZ+0.093*lxdata1$F12+0.056*lxdata1$SERPINE1+
  0.137*lxdata1$COL1A2-0.022*lxdata1$C4BPA-0.031*lxdata1$JMJD7_PLA2G4B-0.017*lxdata1$CR2
lxdata1$Group <- ifelse(lxdata1$Risk>median(lxdata1$Risk),'High risk','Low risk')
############gsea
groupFile <- as.data.frame(lxdata1[,c(1,49)])
write_xlsx(groupFile,path = 'CRG/extdata/GSEA_groupfile.xlsx')

load('CRG/extdata/DEG_TN_TCGA.Rdata')
gseaexp <- exp[rownames(exp)%in%rownames(DEG_limma_voom)[DEG_limma_voom$adj.P.Val<0.05&abs(DEG_limma_voom$logFC)>log2(1)],]
write.table(gesaexp ,'~/Desktop/GSEA_exp.txt',sep = '\t',quote = F)
##############TMB
Mut_tcga=read.table('LUADdata/TCGA/LUAD_mutation.txt.gz',header = T,fill = TRUE)
Mut_tcga$Hugo_Symbol <- Mut_tcga$gene
Mut_tcga$Chromosome <- Mut_tcga$chr
Mut_tcga$Start_Position <- Mut_tcga$start
Mut_tcga$End_Position <- Mut_tcga$end
Mut_tcga$Reference_Allele <- Mut_tcga$reference
Mut_tcga$Tumor_Seq_Allele2 <- Mut_tcga$alt
Mut_tcga$Variant_Classification<- Mut_tcga$effect
Mut_tcga$Variant_Type<- Mut_tcga$DNA_VAF
Mut_tcga$Tumor_Sample_Barcode<- Mut_tcga$sample

library(maftools)
laml = read.maf(maf = Mut_tcga)
tmb <- tmb(maf = laml,logScale=F)
tmb <- tmb[tmb$Tumor_Sample_Barcode%in%lxdata1$sampleID,]
d <- c()
for(i in tmb$Tumor_Sample_Barcode){
  a <- lxdata1$Group[which(lxdata1$sampleID==i)]
  d<- c(d,a)
}
tmb$group <- d
tmb$total_perMB <- as.numeric(tmb$total_perMB)
library(ggpubr)
library(rstatix)

tb <- ggplot(tmb,aes(group,total_perMB,fill=group))+geom_boxplot()+theme_classic()+ylab('TMB')+xlab('')+ggtitle('TCGA')+ylim(0,40)+
  geom_signif(size=1,textsize=12,comparisons = list(c('High risk','Low risk')),test = wilcox.test,map_signif_level = T)+
  theme( plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
         axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1))
ggsave('CRG/result/figs/06_TMB.tiff',tb,width = 5,height =5,dpi=200)
####################CIBERSORT
CIBERSORT <- read_xlsx('CRG/result/files/02_CR2_CIBERSORT.xlsx')
colnames(CIBERSORT)[1] <- 'sampleID'
a <- merge(lxdata1,CIBERSORT,by.x='sampleID')
b <- a[,c(1,49:71)]
rownames(b) <- b$sampleID
pd <-b %>%
  rownames_to_column('sample')%>%
  pivot_longer(cols=c(colnames(b)[-c(1:2)]),names_to = 'Gene',values_to = 'Composition')
pd$Gene <- factor(pd$Gene)
df_p_val1 <- pd %>%
  group_by(Gene) %>%
  t_test(formula = Composition ~ Group) %>%
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>%
  add_xy_position(x='Gene')

pp <- ggboxplot(pd,x = 'Gene', y='Composition',fill='Group')+theme_bw()+ylab('CIBERSORT immune infiltration score')+xlab('')+ylim(min(pd$Composition),max(pd$Composition)+0.2)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  stat_pvalue_manual(size=10,remove.bracket=T,df_p_val1,label = '{p.signif}',tip.length = 0,y.position = 0.8)+
  theme(axis.text.x = element_text(size=18,face='bold',color='black',angle = 45,hjust = 0.8,vjust = 0.9),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/06_CIBERSORT_TCGA.tiff',pp,width = 15,height = 10,dpi=200)
##############Estimate
library(estimate)
library(GSVA)
library(GSEABase)
library(limma)
library(readxl)

gene_set<-read_xlsx('~/IMOA/LUAD/data/immueGenelist.xlsx')
gene_set <- gene_set[,c(1,2)]
colnames(gene_set) <- c('gene','type') 

# length(intersect(var_195_names,gene_set$gene))
exp_cluster1 <- exp[rownames(exp)%in%gene_set$gene,colnames(exp)%in%lxdata1$sampleID[lxdata1$Group=='High risk']]
exp_cluster2 <- exp[rownames(exp)%in%gene_set$gene,colnames(exp)%in%lxdata1$sampleID[lxdata1$Group=='Low risk']]
immueExpr <- exp_cluster1
immueExpr <- exp_cluster2
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
ssGSEA_Score <- gsva(as.matrix(immueExpr),list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}

norm_ssGSEA_Score1 <- normalization(ssGSEA_Score)%>%t%>%as.data.frame()
norm_ssGSEA_Score1$group <- rep('High risk',nrow(norm_ssGSEA_Score1))
norm_ssGSEA_Score2 <- normalization(ssGSEA_Score)%>%t%>%as.data.frame()
norm_ssGSEA_Score2$group <- rep('Low risk',nrow(norm_ssGSEA_Score2))
norm_ssGSEA_Score <- rbind(norm_ssGSEA_Score1,norm_ssGSEA_Score2)
colnames(norm_ssGSEA_Score ) <- c('Act_B','Act_CD4','Act_CD8','Act_DC','CD56bright',
                                  'CD56dim','Tcm_CD4','Tcm_CD8','Tem_CD4','Tem_CD8',
                                  "Eosinophil",'Tgd','Imm_B','iDC',"Macrophage" ,
                                  'Mast',"MDSC",'Mem_B',"Monocyte",'NK',
                                  'NKT',"Neutrophil" ,'pDC','Treg','Tfh',
                                  'Th1','Th17','Th2',"group" )

library(tidyverse)
library(rstatix)
library(ggtext)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
library(tidyr)
fd <-norm_ssGSEA_Score %>% 
  rownames_to_column('sample')%>%
  pivot_longer(cols=c(colnames(norm_ssGSEA_Score)[-ncol(norm_ssGSEA_Score)]),names_to = 'Celltype',values_to = 'Composition')
fd$Celltype <- factor(fd$Celltype)
df_p_val1 <- fd %>% 
  group_by(Celltype) %>% 
  wilcox_test(formula = Composition ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='Celltype')

im<- ggboxplot(fd,x = 'Celltype', y='Composition',fill='group')+theme_bw()+ylab('')+xlab('')+
  scale_fill_manual(values=c("#DC143C", "#0000FF"))+
  stat_pvalue_manual(size=8,remove.bracket=T,y.position = 1.2,vjust=1,df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(axis.text.x = element_text(size=18,face='bold',color='black',angle = 90,hjust = 0.8,vjust = 0.9),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))

ggsave('CRG/result/figs/06_immue.tiff',im,width = 15,height = 10,dpi=200)



exp_cluster1 <- exp[rownames(exp)%in%gene_set$gene,colnames(exp)%in%lxdata1$sampleID[lxdata1$Group=='High risk']]
write.table(exp_cluster1 ,file="CRG/extdata/exp_cluster1.txt",sep = '\t',quote = F)
filterCommonGenes(input.f="CRG/extdata/exp_cluster1.txt",
                  output.f="CRG/extdata/exp_cluster1.gct",
                  id="GeneSymbol")

estimateScore(input.ds = "CRG/extdata/exp_cluster1.gct", 
              output.ds="CRG/extdata/estimateScore_cluster1.gct",
              platform="illumina")
scores1 <- read.table("CRG/extdata/estimateScore_cluster1.gct",skip = 2,header = T)
scores1 <- scores1 %>% t %>% as.data.frame()
colnames(scores1) <- scores1[1,]
scores1 <- scores1[-c(1,2),]
scores1$sample <- gsub('\\.','-',rownames(scores1))
scores1$group <- rep('High risk',nrow(scores1))
exp_cluster2 <- exp[rownames(exp)%in%gene_set$gene,colnames(exp)%in%lxdata1$sampleID[lxdata1$Group=='Low risk']]
write.table(exp_cluster2 ,file="CRG/extdata/exp_cluster2.txt",sep = '\t',quote = F)
filterCommonGenes(input.f="CRG/extdata/exp_cluster2.txt",
                  output.f="CRG/extdata/exp_cluster2.gct",
                  id="GeneSymbol")

estimateScore(input.ds = "CRG/extdata/exp_cluster2.gct", 
              output.ds="CRG/extdata/estimateScore_cluster2.gct",
              platform="illumina")
scores2 <- read.table("CRG/extdata/estimateScore_cluster2.gct",skip = 2,header = T)
scores2 <- scores2 %>% t %>% as.data.frame()
colnames(scores2) <- scores2[1,]
scores2 <- scores2[-c(1,2),]
scores2$sample <- gsub('\\.','-',rownames(scores2))
scores2$group <- rep('Low risk',nrow(scores2))

scores <- rbind(scores1,scores2)
scores$StromalScore <- as.numeric(scores$StromalScore)
scores$ImmuneScore <- as.numeric(scores$ImmuneScore)
scores$ESTIMATEScore <- as.numeric(scores$ESTIMATEScore)
ss <- ggboxplot(scores,x="group",y="StromalScore",fill = "group")+theme_classic()+ylab('StromalScore')+xlab('')+ggtitle('TCGA')+ylim(-100,130)+
  geom_signif(size=1,textsize=12,comparisons = list(c('High risk','Low risk')),test = wilcox.test,map_signif_level = T)+
  theme( plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
         axis.text.x = element_text(size=18,face='bold',color='black'),
         axis.text.y = element_text(size=18,face='bold',color='black'),
         text=element_text(size=25,face='bold'),
         legend.title = element_blank(),legend.position = "none",
         axis.line = element_line(size=1))
is <- ggboxplot(scores,x='group',y='ImmuneScore',fill = "group")+theme_classic()+ylab('ImmuneScore')+xlab('')+ggtitle('TCGA')+ylim(-60,140)+
  geom_signif(size=1,textsize=12,comparisons = list(c('High risk','Low risk')),test = wilcox.test,map_signif_level = T)+
  theme( plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
         axis.text.x = element_text(size=18,face='bold',color='black'),
         axis.text.y = element_text(size=18,face='bold',color='black'),
         text=element_text(size=25,face='bold'),
         legend.title = element_blank(),legend.position = "none",
         axis.line = element_line(size=1))
es <- ggboxplot(scores,x='group',y='ESTIMATEScore',fill = "group")+theme_classic()+ylab('ESTIMATEScore')+xlab('')+ggtitle('TCGA')+ylim(-150,220)+
  geom_signif(size=1,textsize=12,comparisons = list(c('High risk','Low risk')),test = wilcox.test,map_signif_level = T)+
  theme( plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
         axis.text.x = element_text(size=18,face='bold',color='black'),
         axis.text.y = element_text(size=18,face='bold',color='black'),
         text=element_text(size=25,face='bold'),
         legend.title = element_blank(),legend.position = "none",
         axis.line = element_line(size=1))
ggsave('CRG/result/figs/06_stromalscore.tiff',ss,width = 5,height = 5,dpi=200)
ggsave('CRG/result/figs/06_immuescore.tiff',is,width = 5,height = 5,dpi=200)
ggsave('CRG/result/figs/06_estimate.tiff',es,width = 5,height = 5,dpi=200)
##############################MHC
HLA <- c("HLA-F","HLA-E","HLA-G","HLA-DRB5","HLA-DRB1",
         "HLA-DRA","HLA-DQB2","HLA-DQB1","HLA-DQA2","HLA-DQA1",
         "HLA-DPB1","HLA-DPA1","HLA-DOB","HLA-DOA","HLA-DMB","HLA-DMA","HLA-C","HLA-B","HLA-A")
HLA_c1 <- exp[rownames(exp)%in%HLA,colnames(exp)%in%lxdata1$sampleID[lxdata1$Group=='High risk']]
HLA_c2 <- exp[rownames(exp)%in%HLA,colnames(exp)%in%lxdata1$sampleID[lxdata1$Group=='Low risk']]
HLA_c1 <- HLA_c1 %>% t %>% as.data.frame()
HLA_c1$sample <- rownames(HLA_c1)
HLA_c1$group <- rep("High risk",nrow(HLA_c1))
HLA_c2 <- HLA_c2 %>% t %>% as.data.frame()
HLA_c2$sample <- rownames(HLA_c2)
HLA_c2$group <- rep("Low risk",nrow(HLA_c2))
HLA_all <- rbind(HLA_c1,HLA_c2)

library(tidyr)
pd <-HLA_all %>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(HLA_all)[-c(ncol(HLA_all)-1,ncol(HLA_all))]),names_to = 'Celltype',values_to = 'Composition')
pd $Celltype <- factor(pd $Celltype)
df_p_val1 <- pd  %>% 
  group_by(Celltype) %>% 
  wilcox_test(formula = Composition ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='Celltype')

PP<- ggboxplot(pd ,x = 'Celltype', y='Composition',fill='group')+theme_bw()+ylab('')+xlab('')+
  scale_fill_manual(values=c("#DC143C", "#0000FF"))+
  stat_pvalue_manual(size=14,remove.bracket=T,y.position = 25,vjust=1,df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))+coord_flip()

ggsave('CRG/result/figs/06_HLA.tiff',PP,width = 10,height = 12,dpi=200)

##################################TCAF
TCAF <- c('TNFRSF8','TNFRSF25','TNFRSF17','TNFRSF14',
          'ICOS','CD40LG','CD28','CD27','CD226','CD2')
TCAF_c1 <- exp[rownames(exp)%in%TCAF,colnames(exp)%in%lxdata1$sampleID[lxdata1$Group=='High risk']]
TCAF_c2 <- exp[rownames(exp)%in%TCAF,colnames(exp)%in%lxdata1$sampleID[lxdata1$Group=='Low risk']]
TCAF_c1 <- TCAF_c1 %>% t %>% as.data.frame()
TCAF_c1$sample <- rownames(TCAF_c1)

TCAF_c1$group <- rep("High risk",nrow(TCAF_c1))
TCAF_c2 <- TCAF_c2 %>% t %>% as.data.frame()
TCAF_c2$sample <- rownames(TCAF_c2)
TCAF_c2$group <- rep("Low risk",nrow(TCAF_c2))
TCAF_all <- rbind(TCAF_c1,TCAF_c2)
library(tidyr)
library(ggplot2)
td <-TCAF_all %>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(TCAF_all)[-c(ncol(TCAF_all)-1,ncol(TCAF_all))]),names_to = 'Celltype',values_to = 'Composition')
td$Celltype <- factor(td$Celltype)
df_p_val1 <- td  %>% 
  group_by(Celltype) %>% 
  wilcox_test(formula = Composition ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='Celltype')

tP<- ggboxplot(td ,x = 'Celltype', y='Composition',fill='group')+theme_bw()+ylab('')+xlab('')+
  scale_fill_manual(values=c("#DC143C", "#0000FF"))+
  stat_pvalue_manual(size=14,remove.bracket=T,y.position = 15,vjust=1,df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))+coord_flip()
ggsave('CRG/result/figs/06_TCAF.tiff',tP,width = 10,height = 12,dpi=200)

