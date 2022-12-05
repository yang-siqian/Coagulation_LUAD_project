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
bsrgenes[10] <- gsub('_','-',bsrgenes[10])
hubdat <-  exp_CRGs_tumor[,c('sampleID',bsrgenes)]
modeldat<- merge(hubdat,clin_tumor,by.x = 'sampleID')
lxdata <- modeldat
lxdata1 <- merge(exp_CRGs_tumor,clin_tumor,by.x = 'sampleID')
lxdata1$time <- round(lxdata1$time,3)
lxdata1$Age <- ifelse(lxdata$Age=='young',1,2 )
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$Stage <- ifelse(lxdata$Stage=='I',1,ifelse(lxdata$Stage=='II',2 ,ifelse(lxdata$Stage=='III',3,ifelse(lxdata$Stage=='IV',4,0))))
lxdata1$Stage <- factor(lxdata1$Stage,levels=c(0,1,2,3,4),labels=c('UNKOWN','I','II','III','IV'))
bsrgenes[10] <- gsub('-','_',bsrgenes[10])
colnames(lxdata1)[11] <- gsub('-','_',colnames(lxdata1)[11])
exp_tumor <- exp[,colnames(exp)%in%clin_rna$sampleID[clin_rna$group=='Tumor']]
PD <- c('PDCD1','CTLA4','HAVCR2','LAG3','CD274')
exp_BJ_tumor <- exp_tumor[rownames(exp_tumor)%in%c(bsrgenes[-10],gsub('_','-',bsrgenes[10]),PD),]%>%t%>%as.data.frame()
exp_BJ_tumor$sampleID <- rownames(exp_BJ_tumor)

exp <- read.table('LUADdata/TCGA/TCGA-LUAD.htseq_counts.tsv.gz',header=T,row.names = 1)
colnames(exp) <- gsub('\\.','-',colnames(exp))
colnames(exp) <- substr(colnames(exp),1,15)
rownames(exp) <- str_sub(rownames(exp), start = 1, end = 15)
exp$ENSEMBL <- rownames(exp)
HAVCR2 <- exp[rownames(exp)%in%'ENSG00000135077',] %>% t %>% as.data.frame()
HAVCR2$sampleID <- gsub('\\.','-',rownames(HAVCR2))
colnames(HAVCR2)[1] <- 'HAVCR2'
exp_BJ_tumor <- merge(exp_BJ_tumor,HAVCR2,by.x='sampleID')
exp_BJ_tumor <- do.call(data.frame,exp_BJ_tumor)
colnames(exp_BJ_tumor)[16] <- 'JMJD7_PLA2G4B'
colnames(exp_BJ_tumor)[7] <- 'PD-L1'
colnames(exp_BJ_tumor)[12] <- 'PD-1'
lxdata1 <- merge(exp_BJ_tumor,clin_tumor,by.x = 'sampleID')
lxdata1$time <- round(lxdata1$time,3)
lxdata1$Age <- ifelse(lxdata$Age=='young',1,2 )
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$Stage <- ifelse(lxdata$Stage=='I',1,ifelse(lxdata$Stage=='II',2 ,ifelse(lxdata$Stage=='III',3,ifelse(lxdata$Stage=='IV',4,0))))
lxdata1$Stage <- factor(lxdata1$Stage,levels=c(0,1,2,3,4),labels=c('UNKOWN','I','II','III','IV'))

lxdata1$Risk=0.09*lxdata1$F2-0.299*lxdata1$P2RX1-0.106*lxdata1$GUCY1A1+0.107*lxdata1$PLAUR-0.124*lxdata1$PRKCZ+0.093*lxdata1$F12+0.056*lxdata1$SERPINE1+
  0.137*lxdata1$COL1A2-0.022*lxdata1$C4BPA-0.031*lxdata1$JMJD7_PLA2G4B-0.017*lxdata1$CR2
lxdata1$Group <- ifelse(lxdata1$Risk>median(lxdata1$Risk),'High risk','Low risk')
library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
library(tidyr)
library(rstatix)
##############################PD
PDG <- c('PD-1','CTLA4','HAVCR2','LAG3','PD-L1')
rownames(lxdata1) <- lxdata1$sampleID
PDG_c1 <- lxdata1[,colnames(lxdata1)%in%PDG ]
PDG_all <- PDG_c1
PDG_all$sample <- rownames(lxdata1)
PDG_all$group <- lxdata1$Group
PDG_all$HAVCR2 <- as.numeric(PDG_all$HAVCR2)


pd <-PDG_all %>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(PDG_all)[-c(ncol(PDG_all)-1,ncol(PDG_all))]),names_to = 'Celltype',values_to = 'Composition')
pd $Celltype <- factor(pd $Celltype)
df_p_val1 <- pd  %>% 
  group_by(Celltype) %>% 
  wilcox_test(formula = Composition ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='Celltype')

PP<- ggboxplot(pd ,x = 'Celltype', y='Composition',fill='group')+theme_bw()+ylab('Expression')+xlab('')+
  scale_fill_manual(values=c("#DC143C", "#0000FF"))+
  stat_pvalue_manual(size=10,remove.bracket=T,y.position = 16,vjust=1,hjust=1,df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))

ggsave('CRG/result/figs/07_PDG.tiff',PP,width = 8,height =6,dpi=200)

###############################################
if (!requireNamespace("TMEscore", quietly = TRUE))
  devtools::install_github("DongqiangZeng0808/TMEscore")
library('TMEscore')
pd_stad <- lxdata1[,c(1,42,36:37)]
colnames(pd_stad)[1:2] <- c('sampleID','subtype')
eset_stad <- exp_tumor[,colnames( exp_tumor)%in%lxdata1$sampleID]
tmescore<-tmescore(eset = eset_stad ,pdata = pd_stad,column_of_sample = "sampleID",classify = T)
library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
cb <- ggboxplot(tmescore,x = 'subtype',y = 'TMEscore',fill='subtype')+theme_bw()+ylab('TMEscore')+xlab('')+ggtitle('TCGA')+ylim(-20,20)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_signif(size=1,textsize=12,comparisons = list(c("High risk", "Low risk")),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('CRG/result/figs/07_tmescore.tiff',cb,width = 5,height = 5,dpi=200)

cor(tmescore$TMEscore,lxdata1$Risk)

#############################################
IPS <- read_tsv('CRG/extdata/TCIA-ClinicalData.tsv')
IPS$sample <- IPS$barcode
pd_stad$sample <- substr(pd_stad$sampleID,1,12)
IPS <- merge(IPS,pd_stad,by.x='sample')
IPS <- IPS [,c(30:31,25:28)]
colnames(IPS)[2:6] <- c('group','IPS','IPS_PD1','IPS_CTLA4','IPS_CTLA4_PD1')
rownames(IPS) <- IPS$sampleID
IPS <- na.omit(IPS)
library(tidyr)
pd <-IPS %>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(IPS)[-c(1,2)]),names_to = 'Celltype',values_to = 'Composition')
pd $Celltype <- factor(pd $Celltype)
df_p_val1 <- pd  %>% 
  group_by(Celltype) %>% 
  wilcox_test(formula = Composition ~ group) %>% 
  #kruskal_test(formula = Composition ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='Celltype')

PP<- ggboxplot(pd ,x = 'Celltype', y='Composition',fill='group')+theme_bw()+ylab('IPS')+xlab('')+
  scale_fill_manual(values=c("#DC143C", "#0000FF"))+
  stat_pvalue_manual(size=10,remove.bracket=T,y.position =12,vjust=1,hjust=1,df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))

ggsave('CRG/result/figs/07_IPS.tiff',PP,width = 10,height =6,dpi=200)

##########################################
tide <- exp_tumor%>%t%>%as.data.frame()
tide1 <- tide-colMeans(tide)
tide1 <- tide1%>%t%>%as.data.frame()
write.table(tide1[1:2,1:5],'~/Desktop/tide.txt',sep='\t',quote = F)


















