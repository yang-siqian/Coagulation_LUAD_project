library(survival)
library(survminer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(readxl)

exp <- read.table('rawdata/TCGA_LUAD/TCGA-LUAD.htseq_counts.tsv.gz',header=T,row.names = 1)
colnames(exp) <- gsub('\\.','-',colnames(exp))
colnames(exp) <- substr(colnames(exp),1,15)
rownames(exp) <- str_sub(rownames(exp), start = 1, end = 15)
exp$ENSEMBL <- rownames(exp)
HAVCR2 <- exp[rownames(exp)%in%'ENSG00000135077',] %>% t %>% as.data.frame()
HAVCR2$sampleID <- gsub('\\.','-',rownames(HAVCR2))
colnames(HAVCR2)[1] <- 'HAVCR2'
df <- bitr( rownames(exp), fromType = "ENSEMBL", toType = c( "ENTREZID" ), OrgDb = org.Hs.eg.db )
exp <- merge(exp, df, by='ENSEMBL')
gf=bitr(rownames(exp), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
exp <- merge(exp,gf, by='ENTREZID')
exp<-exp[!duplicated(exp$SYMBOL),]
rownames(exp)<-exp$SYMBOL[!duplicated(exp$SYMBOL)]
exp <- exp[,-which(colnames(exp)%in%c( "ENSEMBL", "ENTREZID","SYMBOL"))]
group <- ifelse(as.numeric(substr(colnames(exp),14,15))<10,'Tumor','Normal')
bsrgenes <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7_PLA2G4B','CR2')
bsrgenes[10] <- gsub('_','-',bsrgenes[10])
PD <- c('PDCD1','CTLA4','HAVCR2','LAG3','CD274')
exp_immune <- exp[rownames(exp)%in%c(bsrgenes,PD),] %>% t %>% as.data.frame()
exp_immune$sampleID <- rownames(exp_immune)
exp_immune <- merge(exp_immune,HAVCR2,by.x='sampleID')
exp_immune <- exp_immune %>% t %>% as.data.frame()
colnames(exp_immune) <- exp_immune[1,]

load('Rdata/clinic.Rdata')
clin_tcga <- CLINIC[,c("sampleID","gender","age_at_initial_pathologic_diagnosis","pathologic_stage","pathologic_T",
                       "pathologic_N","pathologic_M","days_to_death" ,"days_to_last_followup","vital_status","tobacco_smoking_history_indicator","new_neoplasm_event_type","radiation_therapy",
                       "additional_pharmaceutical_therapy",'history_of_neoadjuvant_treatment','targeted_molecular_therapy','days_to_additional_surgery_locoregional_procedure',
                       'days_to_additional_surgery_metastatic_procedure')]
colnames(clin_tcga) <- c("sampleID",'Gender','Age','Stage',"Tumor","Lymph_Node","Metastasis","days_to_death" ,"days_to_last_followup",'Status','Smoking','New_tumor_type',"radiation_therapy",
                         'pharmaceutical_therapy', 'neoadjuvant_treatment','targeted_molecular_therapy','additional_surgery_locoregional_days','additional_surgery_metastatic_days')
clin_tcga <- as.data.frame(clin_tcga)
rownames(clin_tcga) <- clin_tcga$sampleID
clin_tcga$group <-  ifelse(as.numeric(substr(clin_tcga$sampleID,14,15))<10,'Tumor','Normal')
clin_tcga$OS <- ifelse(clin_tcga$Status=='LIVING',clin_tcga$days_to_last_followup/30,clin_tcga$days_to_death/30)
clin_tcga$OS.E <- ifelse(clin_tcga$Status=='LIVING',0,1)
exp_imm_tumor <- exp_immune[,colnames(exp_immune)%in%clin_tcga$sampleID[clin_tcga$group=='Tumor']] #515
exp_imm_tumor <- exp_imm_tumor %>% t%>%as.data.frame()

clin_tumor<- clin_tcga[clin_tcga$sampleID%in%exp_imm_tumor$sampleID,]
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
clin_tumor$RFS <- clin_tumor$additional_surgery_locoregional_days
clin_tumor$RFS[which(is.na(clin_tumor$RFS))] <- clin_tumor$days_to_death[which(is.na(clin_tumor$RFS))]
clin_tumor$RFS[which(is.na(clin_tumor$RFS))] <- clin_tumor$days_to_last_followup[which(is.na(clin_tumor$RFS))]
clin_tumor$RFS <- round(clin_tumor$RFS/30,3)
clin_tumor$RFS.E <- clin_tumor$OS.E
clin_tumor$RFS.E[grep('Recurrence',clin_tumor$New_tumor_type)] <- 1
clin_tumor <- clin_tumor[which(clin_tumor$OS!=0),] #502

modeldat<- merge(exp_imm_tumor,clin_tumor,by.x = 'sampleID')
lxdata1 <- modeldat
lxdata1$OS <- round(lxdata1$OS,3)
lxdata1$Age <- ifelse(lxdata$Age=='young',1,2 )
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$Stage <- ifelse(lxdata$Stage=='I',1,ifelse(lxdata$Stage=='II',2 ,ifelse(lxdata$Stage=='III',3,ifelse(lxdata$Stage=='IV',4,0))))
lxdata1$Stage <- factor(lxdata1$Stage,levels=c(0,1,2,3,4),labels=c('UNKOWN','I','II','III','IV'))
colnames(lxdata1) <- gsub('-','_',colnames(lxdata1))
lxdata1$Risk=0.09*as.numeric(lxdata1$F2)-0.299*as.numeric(lxdata1$P2RX1)-0.106*as.numeric(lxdata1$GUCY1A1)+0.107*as.numeric(lxdata1$PLAUR)-0.124*as.numeric(lxdata1$PRKCZ)+
  0.093*as.numeric(lxdata1$F12)+0.056*as.numeric(lxdata1$SERPINE1)+0.137*as.numeric(lxdata1$COL1A2)-0.022*as.numeric(lxdata1$C4BPA)-0.031*as.numeric(lxdata1$JMJD7_PLA2G4B)-0.017*as.numeric(lxdata1$CR2)
lxdata1$Group <- ifelse(lxdata1$Risk>median(lxdata1$Risk),'High risk','Low risk')
colnames(lxdata1)[7] <- 'PD-L1'
colnames(lxdata1)[12] <- 'PD-1'
library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
library(tidyr)
library(rstatix)
##############################PD
# Fig5A
PDG <- c('PD-1','CTLA4','HAVCR2','LAG3','PD-L1')
rownames(lxdata1) <- lxdata1$sampleID
PDG_c1 <- lxdata1[,colnames(lxdata1)%in%PDG ]
PDG_all <- PDG_c1
PDG_all$sample <- rownames(lxdata1)
PDG_all$group <- lxdata1$Group
PDG_all$HAVCR2 <- as.numeric(PDG_all$HAVCR2)
PDG_all$CTLA4 <- as.numeric(PDG_all$CTLA4)
PDG_all$`PD-L1` <- as.numeric(PDG_all$`PD-L1`)
PDG_all$LAG3 <- as.numeric(PDG_all$LAG3)
PDG_all$`PD-1` <- as.numeric(PDG_all$`PD-1`)

pd <-PDG_all %>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(PDG_all)[-c(ncol(PDG_all)-1,ncol(PDG_all))]),names_to = 'Celltype',values_to = 'Composition')
pd$Celltype <- factor(pd $Celltype)
df_p_val1 <- pd  %>% 
  group_by(Celltype) %>% 
  wilcox_test(formula = Composition ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='Celltype')
write.xlsx(PDG_all,file = 'files/Fig5.xlsx',sheetName = 'Fig5A',row.names = T)
write.xlsx(df_p_val1,file = 'files/Fig5.xlsx',sheetName = 'Fig5A_pval',row.names = T,append = T)

PP<- ggboxplot(pd ,x = 'Celltype', y='Composition',fill='group')+theme_bw()+ylab('Expression')+xlab('')+
  scale_fill_manual(values=c("#DC143C", "#0000FF"))+
  stat_pvalue_manual(size=10,remove.bracket=T,y.position = 16,vjust=1,hjust=1,df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))

ggsave('figures/05A_PDG.tiff',PP,width = 8,height =6,dpi=200)

###########################IPS
# Fig5B
IPS <- read_tsv('files/TCIA-ClinicalData.tsv')
IPS$sample <- IPS$barcode
pd_stad <- lxdata1[,c(1,41)]
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
write.xlsx(IPS,file = 'files/Fig5.xlsx',sheetName = 'Fig5B',row.names = F,append = T)
write.xlsx(df_p_val1,file = 'files/Fig5.xlsx',sheetName = 'Fig5B_pval',row.names = F,append = T)
write_xlsx(df_p_val1,path = 'files/Fig5-1.xlsx')

PP<- ggboxplot(pd ,x = 'Celltype', y='Composition',fill='group')+theme_bw()+ylab('IPS')+xlab('')+
  scale_fill_manual(values=c("#DC143C", "#0000FF"))+
  stat_pvalue_manual(size=10,remove.bracket=T,y.position =12,vjust=1,hjust=1,df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))

ggsave('figures/05B_IPS.tiff',PP,width = 10,height =6,dpi=200)

