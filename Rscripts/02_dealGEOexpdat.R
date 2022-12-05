library(dplyr)
library(writexl)
library(GEOquery)
library(stringr)
gene <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7-PLA2G4B','CR2')

gset <- getGEO('GSE72094', GSEMatrix =TRUE, AnnotGPL=FALSE)
exp <- exprs(gset[[1]])%>%as.data.frame()
exp[1:4,1:4]
pd <- pData(gset[[1]])
#判断pd的行名顺序与exp列名完全一致
p <- identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
# 提取芯片平台编号
gpl <- getGEO(gset[[1]]@annotation)
gpl <- gpl@dataTable@table
gpl <- gpl[,c(1,4)]
gene%in%x3$GeneSymbol
exp$ID <- rownames(exp)
x1 <- merge(gpl,exp,by.x='ID') 
x2 <- x1[x1$GeneSymbol!='',]
x3 <- x2[!duplicated(x2$GeneSymbol),]
rownames(x3) <- x3$GeneSymbol
#22115
x4 <- x3[,3:ncol(x3)]
exp442 <- x4
#无RFS
clindata <- pd[,c('geo_accession','gender:ch1','age_at_diagnosis:ch1','smoking_status:ch1','Stage:ch1','survival_time_in_days:ch1','vital_status:ch1')]
colnames(clindata) <- c('sample','gender','age','smoking','stage','time','status')
clindata398 <- clindata[clindata$status!='NA'&clindata$time!='NA',]
exp398 <- exp442[,colnames(exp442)%in%clindata398$sample]
save(exp398,clindata398,file = "LUADdata/GEO/GSE72094.Rdata")
load("LUADdata/GEO/GSE72094.Rdata")
#############################################################JAPAN 117
gset <- getGEO("GSE13213",filename = 'LUADdata/TCGA/GSE13213_series_matrix.txt.gz',destdir = '.', AnnotGPL=FALSE,getGPL =FALSE )
exp <-gset@assayData$exprs %>%as.data.frame()
pd <- gset@phenoData@data
gpl <-getGEO("GPL6480",filename = 'LUADdata/TCGA/GSE13213_family.soft',destdir = '.')
gpl <- gpl@gpls$GPL6480@dataTable@table
gpl <- gpl[,c(1,7)]
exp$ID <- rownames(exp)
x1 <- merge(gpl,exp,by.x='ID') 
x2 <- x1[x1$GENE_SYMBOL!='',]
x3 <- x2[!duplicated(x2$GENE_SYMBOL),]
rownames(x3) <- x3$GENE_SYMBOL
gene%in%x3$GENE_SYMBOL
#19595
x4 <- x3[,3:ncol(x3)]
exp117 <- x4
clindata117 <- pd[,c('geo_accession','Sex:ch1','Age:ch1','Smoking (BI):ch1','Stage (Pathological ):ch1','Survival (days):ch1','Status:ch1')]
colnames(clindata117) <- c('sample','gender','age','smoking','stage','time','status')
save(exp117,clindata117,file = "LUADdata/GEO/GSE13213.Rdata")
load("LUADdata/GEO/GSE13213.Rdata")
#############################################################USA 133
gset <- getGEO("GSE42127", filename = 'LUADdata/TCGA/GSE42127_series_matrix.txt.gz',destdir = '.', AnnotGPL=FALSE,getGPL =FALSE )
exp <-gset@assayData$exprs %>%as.data.frame()
pd <- gset@phenoData@data
gpl <-getGEO("GPL6884",filename = 'LUADdata/TCGA/GSE42127_family.soft',destdir = '.')
gpl <- gpl@gpls$GPL6884@dataTable@table
gpl <- gpl[,c(1,13)]
gene%in%gpl$Symbol
exp$ID <- rownames(exp)
x1 <- merge(gpl,exp,by.x='ID') 
x2 <- x1[x1$Symbol!='',]
x3 <- x2[!duplicated(x2$Symbol),]
rownames(x3) <- x3$Symbol
gene%in%x3$Symbol
#25440
x4 <- x3[,3:ncol(x3)]

clindata133<- pd[pd$`histology:ch1`=="Adenocarcinoma",c('geo_accession','gender:ch1','age at surgery:ch1','final.pat.stage:ch1','overall survival months:ch1','survival status:ch1')]
colnames(clindata133) <- c('sample','gender','age','stage','time','status')
exp133 <- x4[,colnames(x4)%in%clindata133$sample]
exp133 <- na.omit(exp133)
save(exp133,clindata133,file = "LUADdata/GEO/GSE42127.Rdata")
load("LUADdata/GEO/GSE42127.Rdata")

#############################################################USA 443 km 无意义
gset <- getGEO("GSE68465",filename = 'LUADdata/TCGA/GSE68465_series_matrix.txt.gz',destdir = '.', AnnotGPL=FALSE,getGPL =FALSE )
exp <-gset@assayData$exprs %>%as.data.frame()
pd <- gset@phenoData@data
gpl <-getGEO("GPL96",filename = 'LUADdata/TCGA/GSE68465_family.soft',destdir = '.')
gpl <- gpl@gpls$GPL96@dataTable@table
gpl <- gpl[,c(1,11)]
gene%in%gpl$`Gene Symbol`
exp$ID <- rownames(exp)
x1 <- merge(gpl,exp,by.x='ID') 
x2 <- x1[x1$`Gene Symbol`!='',]
x3 <- x2[!duplicated(x2$`Gene Symbol`),]
rownames(x3) <- x3$`Gene Symbol`
gene%in%x3$`Gene Symbol`
#13515
x4 <- x3[,3:ncol(x3)]
#table(pd$`disease_state:ch1`)
clindata443 <- pd[pd$`disease_state:ch1`!='Normal',c('geo_accession','Sex:ch1','age:ch1','smoking_history:ch1','disease_stage:ch1','months_to_last_contact_or_death:ch1','vital_status:ch1')]
colnames(clindata443) <- c('sample','gender','age','smoking','stage','time','status')
exp443 <- x4[,colnames(x4)%in%clindata443$sample]
exp443 <- log2(exp443+1)
gene%in%rownames(exp443)
save(exp443,clindata443,file = "LUADdata/GEO/GSE68465.Rdata")
load("LUADdata/GEO/GSE68465.Rdata")
#############################################################
gset <- getGEO("GSE31210", filename = 'LUADdata/TCGA/GSE31210_series_matrix.txt.gz',destdir = '.', AnnotGPL=FALSE,getGPL =FALSE )
exp <-gset@assayData$exprs %>%as.data.frame()
pd <- gset@phenoData@data
gpl <-getGEO("GPL570",filename = 'LUADdata/TCGA/GSE31210_family.soft',destdir = '.')
gpl <- gpl@gpls$GPL570@dataTable@table
gpl <- gpl[,c(1,11)]
gene%in%gpl$`Gene Symbol`
exp$ID <- rownames(exp)
x1 <- merge(gpl,exp,by.x='ID') 
x2 <- x1[x1$`Gene Symbol`!='',]
x3 <- x2[!duplicated(x2$`Gene Symbol`),]
rownames(x3) <- x3$`Gene Symbol`
gene%in%x3$`Gene Symbol`
#23520
x4 <- x3[,3:ncol(x3)]
table(pd$`tissue:ch1`)
clindata226<- pd[pd$`tissue:ch1`=='primary lung tumor',c('geo_accession','gender:ch1','age (years):ch1','smoking status:ch1','pstage iorii:ch1','months before relapse/censor:ch1','death:ch1')]
colnames(clindata226) <- c('sample','gender','age','smoking','stage','time','status')
exp226 <- x4[,colnames(x4)%in%clindata226$sample]
exp226 <- log2(exp226+1)
clindata226$time <- gsub(';.*','',clindata226$time )
clindata226$status <- ifelse(clindata226$status=='alive','Alive','Dead')
save(exp226,clindata226,file = "LUADdata/GEO/GSE31210.Rdata")
load("LUADdata/GEO/GSE31210.Rdata")
#############################################################
gset <- getGEO("GSE37745",  filename = 'LUADdata/TCGA/GSE37745_series_matrix.txt.gz',destdir = '.', AnnotGPL=FALSE,getGPL =FALSE )
exp <-gset@assayData$exprs %>%as.data.frame()
pd <- gset@phenoData@data
gpl <-getGEO("GPL570",filename = 'LUADdata/TCGA/GSE37745_family.soft',destdir = '.')
gpl <- gpl@gpls$GPL570@dataTable@table
gpl <- gpl[,c(1,11)]
gene%in%gpl$`Gene Symbol`
exp$ID <- rownames(exp)
x1 <- merge(gpl,exp,by.x='ID') 
x2 <- x1[x1$`Gene Symbol`!='',]
x3 <- x2[!duplicated(x2$`Gene Symbol`),]
rownames(x3) <- x3$`Gene Symbol`
gene%in%x3$`Gene Symbol`
#23520
x4 <- x3[,3:ncol(x3)]
table(pd$`histology:ch1`)
clindata106 <- pd[pd$`histology:ch1`=='adeno',c('geo_accession','gender:ch1','age:ch1','tumor stage:ch1','days to determined death status:ch1','dead:ch1')]
colnames(clindata106) <- c('sample','gender','age','stage','time','status')
table(clindata106$status)
clindata106$status <- ifelse(clindata106$status=='no','Alive','Dead')

exp106 <- x4[,colnames(x4)%in%clindata106$sample]

save(exp106,clindata106,file = "LUADdata/GEO/GSE37745.Rdata")
load("LUADdata/GEO/GSE37745.Rdata")
############################################################
gset <- getGEO("GSE19151", GSEMatrix =TRUE, AnnotGPL=FALSE)
gpl <- gpl[,c(1,11)]
x2 <- x1[x1$`Gene Symbol`!='',]
x3 <- x2[!duplicated(x2$`Gene Symbol`),]
rownames(x3) <- x3$`Gene Symbol`
#13515
x4 <- x3[,3:ncol(x3)]
clindata133 <- pd[,c('geo_accession','gender:ch1','age:ch1','source_name_ch1')]
colnames(clindata133) <- c('sample','gender','age','group')
clindata133$group[grep('VTE',clindata133$group)] <- 'VTE'
clindata133$group <- ifelse(clindata133$group =='VTE','VTE','Control')
exp133 <- x4[,clindata133$sample]
save(exp133,clindata133,file = "LUADdata/GEO/GSE19151.Rdata")
load("LUADdata/GEO/GSE19151.Rdata")
# Control     VTE 
# 63           70 （18 PE）
###########################################################
gset <- getGEO("GSE48000", GSEMatrix =TRUE, AnnotGPL=FALSE)
gpl <- gpl[,c(1,13)]
x2 <- x1[x1$Symbol!='',]
x3 <- x2[!duplicated(x2$Symbol),]
rownames(x3) <- x3$Symbol
#31412
x4 <- x3[,3:ncol(x3)]
clindata133 <- pd[,c('geo_accession','phenotype/group:ch1')]
colnames(clindata133) <- c('sample','group')
clindata133$group[grep('Risk',clindata133$group)] <- 'VTE'
clindata133$group <- ifelse(clindata133$group =='VTE','VTE','Control')
exp133 <- x4[,clindata133$sample]
save(exp133,clindata133,file = "LUADdata/GEO/GSE48000.Rdata")
load("LUADdata/GEO/GSE48000.Rdata")
# Control     VTE
# 25     107（58 PE）




#免疫治疗数据
############################################################黑色素瘤 PD1
library(readxl)
response <- read_xlsx('LUADdata/TCGA/melanoma.xlsx')
m_naive <- read.table('LUADdata/TCGA/Riaz2017_PD1_Melanoma_RNASeq_Ipi.Naive/ICB.Riaz2017_Nivolumab_Melanoma_Naive.self_subtract')
m_prog <- read.table('LUADdata/TCGA/Riaz2017_PD1_Melanoma_RNASeq_Ipi.Prog/ICB.Riaz2017_Nivolumab_Melanoma_Prog.self_subtract')
m_naive_clin <- read.table('LUADdata/TCGA/Riaz2017_PD1_Melanoma_RNASeq_Ipi.Naive/ICB.Riaz2017_Nivolumab_Melanoma_Naive.clinical',sep='\t',header=T,row.names=1)
m_prog_clin <- read.table('LUADdata/TCGA/Riaz2017_PD1_Melanoma_RNASeq_Ipi.Prog//ICB.Riaz2017_Nivolumab_Melanoma_Prog.clinical',sep='\t',header=T,row.names=1)
clin <- rbind(m_naive_clin ,m_prog_clin)
melanoma <- cbind(m_naive,m_prog)
response$Group <- ifelse(response$Response%in%c('PD','SD'),'SD/PD',ifelse(response$Response%in%c('CR','PR'),'CR/PR','NE'))
response <- response[response$Group!='NE',]
colnames(response)[1] <- 'sample'
library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

melanoma$ENTREZID <- rownames(melanoma)
df <- bitr( rownames(melanoma), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
head(df)
melanoma <- merge(melanoma ,df,by.x='ENTREZID')
rownames(melanoma) <- melanoma$SYMBOL
bsrgenes <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7-PLA2G4B','CR2')
bsrgenes%in%melanoma$SYMBOL
exp <- melanoma[melanoma$SYMBOL%in%bsrgenes,-c(1,53)]
rownames(exp) <- gsub('-','_',rownames(exp) )
exp1 <- exp %>% t %>% as.data.frame()
exp1$Risk=-0.299*exp1$P2RX1-0.106*exp1$GUCY1A1+0.107*exp1$PLAUR-0.124*exp1$PRKCZ+0.093*exp1$F12+0.056*exp1$SERPINE1+
  0.137*exp1$COL1A2-0.022*exp1$C4BPA-0.031*exp1$JMJD7_PLA2G4B-0.017*exp1$CR2
exp1$sample <- rownames(exp1)
clin$sample <- rownames(clin)
survdata <- merge(exp1,clin,by.x='sample')
survdata$time <- round(survdata$OS/30,2)
survdata$Time <- round(survdata$PFS/30,2)
survdata$group <- ifelse(survdata$Risk>median(survdata$Risk),'High risk','Low risk')
response1 <- response[,c(1,13)]
survdata1 <- merge(survdata,response1,by.x='sample')

png("CRG/result/figs/03_OSKM_risk_melanoma.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(time, OS.Event) ~ group, data=survdata1),size=3,
                surv.median.line = "hv",
                title = "Melanoma",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("pink", "purple",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("CRG/result/figs/03_RFSKM_risk_melanoma.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(Time, PFS.Event) ~ group, data=survdata1),size=3,
                surv.median.line = "hv",
                title = "Melanoma",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("lightblue", "yellowgreen",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of PFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
#风险评分box图
cb <- ggboxplot(survdata1,x = 'Group',y = 'Risk',fill='Group')+theme_bw()+ylab('Riskscore')+xlab('')+ggtitle('Melanoma')+ylim(-2,3)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_signif(size=1,textsize=12,comparisons = list(c("CR/PR", "SD/PD")),test = t.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('CRG/result/figs/03_riskbox.tiff',cb,width = 5,height = 5,dpi=200)
library(reshape2)
a <-data.frame(group=c('High risk','Low risk'),CR_PR=c(21.7,19.2),SD_PD=c(78.3,80.8))
b=melt(a,id.vars='group')
b$variable <- gsub('_','/',b$variable)
P <- ggplot(data=b,aes(x=group,y=value,fill=variable)) + geom_col(position = 'stack', width = 0.6)+
  geom_bar(position = "stack", stat = "identity", width = 0.6) +theme_classic()+theme(legend.title = element_blank())

###############################################################################黑色素瘤  PD1
response <- read_xlsx('LUADdata/TCGA/GSE78220.xlsx')
Hugo <- read.table('LUADdata/TCGA/Hugo2016_PD1_Melanoma_RNASeq/ICB.Hugo2016_Pembrolizumab_Melanoma.self_subtract',header = T,row.names = 1)
clin <- read.table('LUADdata/TCGA/Hugo2016_PD1_Melanoma_RNASeq/ICB.Hugo2016_Pembrolizumab_Melanoma.clinical',header = T,row.names = 1,sep = '\t')
response$Group <- ifelse(response$irRECIST%in%c('Progressive Disease'),'SD/PD',ifelse(response$irRECIST%in%c('Complete Response','Partial Response'),'CR/PR','NE'))
response <- response[response$Group!='NE',]
colnames(response)[1] <- 'sample'

Hugo$ENTREZID <- rownames(Hugo)
df <- bitr( rownames(Hugo), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
Hugo <- merge(Hugo ,df,by.x='ENTREZID')
rownames(Hugo) <- Hugo$SYMBOL
bsrgenes <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7-PLA2G4B','CR2')
bsrgenes%in%df$SYMBOL
exp <- Hugo[Hugo$SYMBOL%in%bsrgenes,-c(1,28)]
exp1 <- exp %>% t %>% as.data.frame()
exp1$Risk=-0.299*exp1$P2RX1-0.106*exp1$GUCY1A1+0.107*exp1$PLAUR-0.124*exp1$PRKCZ+0.093*exp1$F12+0.056*exp1$SERPINE1+
  0.137*exp1$COL1A2-0.022*exp1$C4BPA-0.017*exp1$CR2
exp1$sample <- rownames(exp1)
clin$sample <- rownames(clin)
survdata <- merge(exp1,clin,by.x='sample')
survdata$group <- ifelse(survdata$Risk>median(survdata$Risk),'High risk','Low risk')
response1 <- response[,c(1,11)]
survdata1 <- merge(survdata,response1,by.x='sample')
survdata1$time <- round(survdata1$OS/30,2)

png("CRG/result/figs/09_OSKM_risk_GSE78220.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(time, OS.Event) ~ group, data=survdata1),size=3,
                surv.median.line = "hv",
                title = "GSE78220",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()


library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
#风险评分box图
survdata1$Group <- factor(survdata1$Group,levels = c("CR/PR", "SD/PD"))
cb <- ggboxplot(survdata1,x = 'Group',y = 'Risk',fill='Group')+theme_bw()+ylab('Riskscore')+xlab('')+ggtitle('GSE78220')+ylim(-2,2)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_signif(size=1,textsize=12,comparisons = list(c("CR/PR", "SD/PD")),test = t.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('CRG/result/figs/09_riskbox_GSE78220.tiff',cb,width = 5,height = 5,dpi=200)
library(reshape2)
table(survdata1$Group[survdata1$group=='High risk'])
table(survdata1$Group[survdata1$group=='Low risk'])
a <-data.frame(group=c('High risk','Low risk'),CR_PR=c(38.5,69.2),SD_PD=c(61.5,30.8))
b=melt(a,id.vars='group')
b$variable <- gsub('_','/',b$variable)

P <- ggplot(data=b,aes(x=group,y=value,fill=variable)) + theme_bw()+ylab('Percent')+xlab('')+ggtitle('GSE78220')+ylim(0,100)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_col(position = 'stack', width = 0.6)+
  geom_bar(position = "stack", stat = "identity", width = 0.6) +
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('CRG/result/figs/09_barplot_GSE78220.tiff',P,width = 5,height = 5,dpi=200)

#######################################################################################肺癌+黑色素瘤 PD1
Prat <- read.table('LUADdata/TCGA/Prat2017_PD1_NSCLC-HNSC-Melanoma_Nanostring/ICB.Prat2017_Pembrolizumab-Nivolumab_NSCLC-HNSC-Melanoma.self_subtract',header = T,row.names = 1)
clin <- read.table('LUADdata/TCGA/Prat2017_PD1_NSCLC-HNSC-Melanoma_Nanostring/ICB.Prat2017_Pembrolizumab-Nivolumab_NSCLC-HNSC-Melanoma.clinical',header = T)
Prat$ENTREZID <- rownames(Prat)
df <- bitr( rownames(Prat), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
Prat<- merge(Prat ,df,by.x='ENTREZID')
rownames(Prat) <- Prat$SYMBOL
bsrgenes <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7-PLA2G4B','CR2')
bsrgenes%in%df$SYMBOL
exp <- Prat[Prat$SYMBOL%in%bsrgenes,-c(1,35)]
exp1 <- exp %>% t %>% as.data.frame()
exp1$Risk=0.107*exp1$PLAUR+0.093*exp1$F12-0.022*exp1$C4BPA-0.017*exp1$CR2
exp1$sample <- rownames(exp1)
clin$sample <- rownames(clin)
survdata <- merge(exp1,clin,by.x='sample')
survdata$Group <- ifelse(survdata$RECIST%in%c('PD','SD'),'SD/PD',ifelse(survdata$RECIST%in%c('CR','PR'),'CR/PR','NE'))
survdata <- survdata[survdata$Group!='NE',]
survdata$group <- ifelse(survdata$Risk>median(survdata$Risk),'High risk','Low risk')
# survdata1 <- survdata[survdata$Type%in%c('ADENO','SQUAMOUS'),]
# survdata1$group <- ifelse(survdata1$Risk>median(survdata1$Risk),'High risk','Low risk')

png("CRG/result/figs/09_PFSKM_risk_GSE93157.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(PFS, PFS.Event) ~ group, data=survdata),size=3,
                surv.median.line = "hv",
                title = "GSE93157",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of PFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
#风险评分box图
survdata$Group <- factor(survdata$Group,levels = c("CR/PR", "SD/PD"))
cb <- ggboxplot(survdata,x = 'Group',y = 'Risk',fill='Group')+theme_bw()+ylab('Riskscore')+xlab('')+ggtitle('GSE93157')+ylim(-0.15,0.15)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_signif(size=1,textsize=12,comparisons = list(c("CR/PR", "SD/PD")),test = t.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('CRG/result/figs/09_riskbox_GSE93157.tiff',cb,width = 5,height = 5,dpi=200)
library(reshape2)
table(survdata$Group[survdata$group=='High risk'])
table(survdata$Group[survdata$group=='Low risk'])
a <-data.frame(group=c('High risk','Low risk'),CR_PR=c(43.8,35.3),SD_PD=c(56.2,64.7))
b=melt(a,id.vars='group')
b$variable <- gsub('_','/',b$variable)
P <- ggplot(data=b,aes(x=group,y=value,fill=variable)) + theme_bw()+ylab('Percent')+xlab('')+ggtitle('GSE93157')+ylim(0,100)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_col(position = 'stack', width = 0.6)+
  geom_bar(position = "stack", stat = "identity", width = 0.6) +
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('CRG/result/figs/09_barplot_GSE93157.tiff',P,width = 5,height = 5,dpi=200)

