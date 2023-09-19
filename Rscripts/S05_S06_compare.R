library(survival)
library(survminer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(readxl)
CRGs <- read_xlsx('files/01_CRGs_209.xlsx')
exp <- read.table('rawdata/TCGA_LUAD/TCGA-LUAD.htseq_counts.tsv.gz',header=T,row.names = 1)
colnames(exp) <- gsub('\\.','-',colnames(exp))
colnames(exp) <- substr(colnames(exp),1,15)
rownames(exp) <- str_sub(rownames(exp), start = 1, end = 15)
exp$ENSEMBL <- rownames(exp)
df <- bitr( rownames(exp), fromType = "ENSEMBL", toType = c( "ENTREZID" ), OrgDb = org.Hs.eg.db )
exp <- merge(exp, df, by='ENSEMBL')
gf=bitr(rownames(exp), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
exp <- merge(exp,gf, by='ENTREZID')
exp<-exp[!duplicated(exp$SYMBOL),]
rownames(exp)<-exp$SYMBOL[!duplicated(exp$SYMBOL)]
exp <- exp[,-which(colnames(exp)%in%c( "ENSEMBL", "ENTREZID","SYMBOL"))]
group <- ifelse(as.numeric(substr(colnames(exp),14,15))<10,'Tumor','Normal')
CRRS <- lasso_genes
Jin <- c('SERPINA1','CFHR3','PPP1CB','P2RX1','PLCB3','PLCB4','PIK3R1','GP1BA')
Chen <- c('SERPINE1','VWF','F2R','ANXA5','CD59','AXL','MMRN1')
Jia <- c('SERPINA1', 'HMGCS2', 'MMP7', 'PLAT')
Song <- c('ANG', 'C1QA', 'CFB', 'DUSP6', 'KLKB1', 'MMP7', 'RABIF')
hub <- unique(c(CRRS,Jin,Chen,Jia,Song))
exp_hub <- exp[rownames(exp)%in%hub,]
exp_CRGs <- exp[rownames(exp)%in%CRGs$gene,]
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

exp_hub_tumor <- exp_hub[,colnames(exp_hub)%in%clin_tcga$sampleID[clin_tcga$group=='Tumor']] #515
exp_hub_tumor <- exp_hub_tumor %>% t%>%as.data.frame()
exp_hub_tumor$sampleID <- rownames(exp_hub_tumor)

# exp_CRGs_tumor <- exp_CRGs[,colnames(exp_CRGs)%in%clin_tcga$sampleID[clin_tcga$group=='Tumor']] #515
# exp_CRGs_tumor <- exp_CRGs_tumor %>% t%>%as.data.frame()
# exp_CRGs_tumor$sampleID <- rownames(exp_CRGs_tumor)

# clin_tumor<- clin_tcga[clin_tcga$sampleID%in%exp_CRGs_tumor$sampleID,]
clin_tumor<- clin_tcga[clin_tcga$sampleID%in%exp_hub_tumor$sampleID,]

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
clin_tumor$ps1 <- ifelse(clin_tumor$Age=='old',1,0)
clin_tumor$ps2 <- ifelse(clin_tumor$Metastasis=='M0',0,3)
clin_tumor$ps3 <- ifelse(clin_tumor$radiation_therapy=='YES',3,0)
#clin_tumor$ps4 <- ifelse(clin_tumor$pharmaceutical_therapy=='YES',3,0)
#clin_tumor$ps5 <- ifelse(clin_tumor$targeted_molecular_therapy=='YES',3,0)
#clin_tumor$ps6 <- ifelse(clin_tumor$additional_surgery_locoregional_days!='NA',2,0)
clin_tumor$ps <- ifelse(clin_tumor$ps2==3|clin_tumor$ps3==3,5,2)
clin_tumor$ps <- clin_tumor$ps+clin_tumor$ps1
# lxdata1 <- merge(exp_CRGs_tumor,clin_tumor,by.x = 'sampleID')
lxdata1 <- merge(exp_hub_tumor,clin_tumor,by.x = 'sampleID')
colnames(lxdata1) <- gsub('-','_',colnames(lxdata1))
#Jin
mulcox1<- coxph(Surv(OS, OS.E) ~ SERPINA1+CFHR3+PPP1CB+P2RX1+PLCB3+PLCB4+PIK3R1+GP1BA, data = lxdata1 )
lxdata1$Jin = predict(mulcox1,lxdata1)
#Chen
mulcox2<- coxph(Surv(OS, OS.E) ~ SERPINE1+VWF+F2R+ANXA5+CD59+AXL+MMRN1, data = lxdata1 )
lxdata1$Chen=predict(mulcox2,lxdata1)
#Jia
mulcox3<- coxph(Surv(OS, OS.E) ~ SERPINA1+HMGCS2+MMP7+PLAT, data = lxdata1 )
lxdata1$Jia=predict(mulcox3,lxdata1)
#Song
mulcox4<- coxph(Surv(OS, OS.E) ~ ANG+C1QA+CFB+DUSP6+KLKB1+MMP7+RABIF, data = lxdata1 )
lxdata1$Song=predict(mulcox4,lxdata1)
#CRRS
mulcox5<- coxph(Surv(OS, OS.E) ~  COL1A2+F2+P2RX1+PLAUR+C4BPA, data = lxdata1 )
lxdata1$CRRS=predict(mulcox5,lxdata1)
library(xlsx)
write.xlsx(lxdata1[,c(1,53:54,57:61)],file = 'files/SupFig5.xlsx',row.names = F)
library(rms)
library(survivalROC)
library(dplyr)
troc.1= survivalROC(Stime=lxdata1$OS,  
                    status=lxdata1$OS.E,   
                    marker = lxdata1$Jin,     
                    predict.time = 12, method="KM")
troc.2= survivalROC(Stime=lxdata1$OS,  
                    status=lxdata1$OS.E,   
                    marker = lxdata1$Chen,     
                    predict.time = 12, method="KM")
troc.3= survivalROC(Stime=lxdata1$OS,  
                    status=lxdata1$OS.E,   
                    marker = lxdata1$Jia,     
                    predict.time = 12, method="KM")
troc.4= survivalROC(Stime=lxdata1$OS,  
                    status=lxdata1$OS.E,   
                    marker = lxdata1$Song,     
                    predict.time = 12, method="KM")
troc.5= survivalROC(Stime=lxdata1$OS,  
                    status=lxdata1$OS.E,   
                    marker = lxdata1$CRRS,     
                    predict.time = 12, method="KM")
png("figures/S05_roc.png",width = 1500,height = 1500,res=200)
plot(troc.1$FP, troc.1$TP, type="l",col="orange", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
title(main='TCGA',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.2$FP, troc.2$TP, type="l",lwd=4,col="blue", xlim=c(0,1), ylim=c(0,1))
lines(troc.3$FP, troc.3$TP, type="l",lwd=4,col="green", xlim=c(0,1), ylim=c(0,1))
lines(troc.4$FP, troc.4$TP, type="l",lwd=4,col="purple", xlim=c(0,1), ylim=c(0,1))
lines(troc.5$FP, troc.5$TP, type="l",lwd=4,col="red", xlim=c(0,1), ylim=c(0,1))

legend(0.55,0.58,c(paste('Jin=',round(troc.1$AUC,3)),
                   paste('Chen=',round(troc.2$AUC,3)),
                   paste('Jia=',round(troc.3$AUC,3)),
                   paste('Song=',round(troc.4$AUC,3)),
                   paste('CRRS=',round(troc.5$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('orange','blue',"green","purple",'red'),
       bty='n',seg.len = 1.5,cex = 1.5)

dev.off()

lxdata1$Group <- ifelse(lxdata1$Jin>median(lxdata1$Jin),'High risk','Low risk')
lxdata1$Group <- ifelse(lxdata1$Chen>median(lxdata1$Chen),'High risk','Low risk')
lxdata1$Group <- ifelse(lxdata1$Jia>median(lxdata1$Jia),'High risk','Low risk')
lxdata1$Group <- ifelse(lxdata1$Song>median(lxdata1$Song),'High risk','Low risk')
lxdata1$Group <- ifelse(lxdata1$CRRS>median(lxdata1$CRRS),'High risk','Low risk')
for(i in c('Jin','Chen','Jia','Song','CRRS')){

  lxdata1$Group <- ifelse(as.numeric(lxdata1[,i])>median(as.numeric(lxdata1[,i])),'High risk','Low risk')
  png(paste0("figures/S05_",i,"_OSKM.png"),width = 1500,height = 1500,res=200)
  ggsurvplot(survfit(Surv(OS,OS.E) ~ Group, data=lxdata1),
             surv.median.line = "hv",size=3,
             title = paste0(i," model"),
             legend.labs = c('High risk','Low risk'),legend='top',font.legend=c(20,'bold'),
             palette = c("#E7B800", "#2E9FDF",'grey'),
             pval = T,pval.size=18,legend.title = "",
             ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
             font.main = c(30, "bold"),
             ylab='Survival probability of OS',xlab='Time(month)',
             font.x = c(30, "bold"),
             font.y = c(30, "bold"))
  dev.off()
}


library(pec)
library(survival)

cox1 <- coxph(Surv(OS, OS.E) ~ CRRS, data = lxdata1,x=T,y=T)
cox2 <- coxph(Surv(OS, OS.E) ~ Jin, data = lxdata1,x=T,y=T)
cox3 <- coxph(Surv(OS, OS.E) ~ Chen, data = lxdata1,x=T,y=T)
cox4 <- coxph(Surv(OS, OS.E) ~ Jia, data = lxdata1,x=T,y=T)
cox5 <- coxph(Surv(OS, OS.E) ~ Song, data = lxdata1,x=T,y=T)

ApparrentCindex  <- pec::cindex(list("CRRS"=cox1,
                                     "Jin"=cox2,
                                     "Chen"=cox3,
                                     "Jia"=cox4,
                                     "Song"=cox5),
                                # formula=Surv(OS, OS.E)~Jin,
                                data=lxdata1,
                                eval.times=seq(1,242,20))

png("figures/S05_cindex.png",width = 1500,height = 1500,res=200)
plot(ApparrentCindex,lwd=4,cex.axis=2,legend.x=190,legend.y=1,legend.cex=1,col=c('red',"purple3",'blue',"green",'orange'))
dev.off()

########################PADUA_scores
lxdata <- lxdata1[,c(1,205:214)]
table(lxdata$ps)
library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
library(reshape2)
a <-data.frame(group=c('High risk','Low risk'),High_PADUAscores=c(41.2,44.1),Low_PADUAscores=c(58.8,55.9))
b=melt(a,id.vars='group')

P <- ggplot(data=b,aes(x=group,y=value,fill=variable)) + theme_bw()+ylab('Percent')+xlab('')+ggtitle('PADUA scores')+ylim(0,100)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_col(position = 'stack', width = 0.6)+
  geom_bar(position = "stack", stat = "identity", width = 0.6) +
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 18),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('figures/PADUA_scores.png',P,width = 6,height = 5,dpi=200)


################################SupFig 6
library(readxl)
oncogenes <- c('EGFR','ALK','TP53','KRAS','RET','MET','ERBB2','BRAF','ROS1')
Mut_tcga=read.table('rawdata/TCGA_LUAD/LUAD_mutation.txt.gz',header = T,fill = TRUE)
Mut_tcga$Hugo_Symbol <- Mut_tcga$gene
Mut_tcga$Chromosome <- Mut_tcga$chr
Mut_tcga$Start_Position <- Mut_tcga$start
Mut_tcga$End_Position <- Mut_tcga$end
Mut_tcga$Reference_Allele <- Mut_tcga$reference
Mut_tcga$Tumor_Seq_Allele2 <- Mut_tcga$alt
Mut_tcga$Variant_Classification<- Mut_tcga$effect
Mut_tcga$Variant_Type<- Mut_tcga$DNA_VAF
Mut_tcga$Tumor_Sample_Barcode<- Mut_tcga$sample
Mut_tcga_high <- Mut_tcga[Mut_tcga$sample%in%lxdata1$sampleID[lxdata1$Group=='High risk'],]  
Mut_tcga_low <- Mut_tcga[Mut_tcga$sample%in%lxdata1$sampleID[lxdata1$Group=='Low risk'],]  
library(maftools)
laml = read.maf(maf = Mut_tcga_high)
laml = read.maf(maf = Mut_tcga_low)
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
oncoplot(maf = laml,genes = oncogenes,
         sortByAnnotation = TRUE, 
         bgCol='white',
         colors = vc_cols)

#######################################################
library(dplyr)
library(writexl)
library(GEOquery)
library(stringr)
gset <- getGEO("GSE100797",destdir = '.', AnnotGPL=FALSE,getGPL =FALSE )
pd <- gset$GSE100797_series_matrix.txt.gz@phenoData@data
