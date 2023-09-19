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
bsrgenes <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','ITGB1','C4BPA','JMJD7_PLA2G4B','CR2')
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
lxdata1$Age <- ifelse(lxdata1$Age=='young',1,2 )
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$Stage <- ifelse(lxdata1$Stage=='I',1,ifelse(lxdata1$Stage=='II',2 ,ifelse(lxdata1$Stage=='III',3,ifelse(lxdata1$Stage=='IV',4,0))))
lxdata1$Stage <- factor(lxdata1$Stage,levels=c(0,1,2,3,4),labels=c('UNKOWN','I','II','III','IV'))
colnames(lxdata1) <- gsub('-','_',colnames(lxdata1))
library(survival)
library(survminer)
library(rms)
library(survivalROC)
library(dplyr)

lxdata1$Risk=0.083*as.numeric(lxdata1$F2)-0.256*as.numeric(lxdata1$P2RX)+0.127*as.numeric(lxdata1$PLAUR)-0.163*as.numeric(lxdata1$PRKCZ)+0.098*as.numeric(lxdata1$F12)+0.083*as.numeric(lxdata1$SERPINE1)+
  0.076*as.numeric(lxdata1$ITGB1)-0.029*as.numeric(lxdata1$C4BPA)-0.021*as.numeric(lxdata1$CR2)-0.055*as.numeric(lxdata1$GUCY1A1)-0.035*as.numeric(lxdata1$JMJD7_PLA2G4B)
lxdata1$Group <- ifelse(lxdata1$Risk>median(lxdata1$Risk),'High risk','Low risk')

immunotherapy_score <- function(exp,lxdata1){
  genes <- c('PDCD1','CTLA4','CD274')
  c1 <- exp[rownames(exp)%in%genes,colnames(exp)%in%lxdata1$sampleID[lxdata1$Group=='High risk']]
  c2 <- exp[rownames(exp)%in%genes,colnames(exp)%in%lxdata1$sampleID[lxdata1$Group=='Low risk']]
  c1 <- c1 %>% t %>% as.data.frame()
  c1$sample <- rownames(c1)
  c1$group <- rep("High risk",nrow(c1))
  c2 <- c2 %>% t %>% as.data.frame()
  c2$sample <- rownames(c2)
  c2$group <- rep("Low risk",nrow(c2))
  all <- rbind(c1,c2)
  colnames(all)[3] <- 'PD-L1'
  colnames(all)[2] <- 'PD-1'
  return(all)
}


scores <- immunotherapy_score(exp,lxdata1)
write.xlsx(scores,file = 'figures/06_data.xlsx',sheetName = 'PDG',row.names = F,append = T)

immunotherapy_plot <- function(scores){
  library(dplyr)
  library(ggstatsplot)
  library(tidyverse)
  library(ggpubr)
  library(tidyr)
  library(rstatix)
  pd <-scores %>% 
    rownames_to_column('Sample')%>%
    pivot_longer(cols=c(colnames(scores)[-c(ncol(scores)-1,ncol(scores))]),names_to = 'Celltype',values_to = 'Composition')
  pd$Celltype <- factor(pd $Celltype)
  df_p_val <- pd  %>% 
    group_by(Celltype) %>% 
    wilcox_test(formula = Composition ~ group) %>% 
    add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
    add_xy_position(x='Celltype')
  p<- ggboxplot(pd ,x = 'Celltype', y='Composition',fill='group')+theme_bw()+ylab('Expression')+xlab('')+
    scale_fill_manual(values=c("#DC143C", "#0000FF"))+
    stat_pvalue_manual(size=10,remove.bracket=T,y.position = 16,vjust=1,hjust=1,df_p_val,label = '{p.signif}',tip.length = 0)+
    theme(axis.text.x = element_text(size=18,face='bold',color='black'),
          axis.text.y = element_text(size=18,face='bold',color='black'),
          text=element_text(size=25,face='bold'),
          legend.title = element_blank(),legend.position = "top",
          panel.grid = element_blank(),panel.border = element_rect(size=2))
  return(p)
}
p <- immunotherapy_plot(scores)
p
ggsave('figures/06A_PDG.png',p,width = 8,height =6,dpi=200)


###########################IPS
IPS_score <- function(){
  IPS <- read_tsv('files/TCIA-ClinicalData.tsv')
  IPS$sample <- IPS$barcode
  pd_stad <- lxdata1[,c('sampleID','Group')]
  pd_stad$sample <- substr(pd_stad$sampleID,1,12)
  IPS <- merge(IPS,pd_stad,by.x='sample')
  IPS <- IPS [,c(30:31,25:28)]
  colnames(IPS)[2:6] <- c('group','IPS','IPS_PD1','IPS_CTLA4','IPS_CTLA4_PD1')
  rownames(IPS) <- IPS$sampleID
  IPS <- na.omit(IPS)
  return(IPS)
}
scores <- IPS_score()
write.xlsx(scores,file = 'figures/06_data.xlsx',sheetName = 'IPS',row.names = F,append = T)
IPS_plot <- function(scores){
  pd <-scores %>% 
    rownames_to_column('Sample')%>%
    pivot_longer(cols=c(colnames(scores)[-c(1,2)]),names_to = 'Celltype',values_to = 'Composition')
  pd $Celltype <- factor(pd $Celltype)
  df_p_val <- pd  %>% 
    group_by(Celltype) %>% 
    wilcox_test(formula = Composition ~ group) %>% 
    add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
    add_xy_position(x='Celltype')
  p<- ggboxplot(pd ,x = 'Celltype', y='Composition',fill='group')+theme_bw()+ylab('IPS')+xlab('')+
    scale_fill_manual(values=c("#DC143C", "#0000FF"))+
    stat_pvalue_manual(size=10,remove.bracket=T,y.position =12,vjust=1,hjust=1,df_p_val,label = '{p.signif}',tip.length = 0)+
    theme(axis.text.x = element_text(size=18,face='bold',color='black'),
          axis.text.y = element_text(size=18,face='bold',color='black'),
          text=element_text(size=25,face='bold'),
          legend.title = element_blank(),legend.position = "top",
          panel.grid = element_blank(),panel.border = element_rect(size=2))
  
}
p <- IPS_plot(scores)
ggsave('figures/06_IPS.png',p,width = 10,height =6,dpi=200)


