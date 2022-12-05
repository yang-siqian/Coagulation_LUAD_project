library(readxl)
CRGs <- read_xlsx('CRG/result/files/CRGs_209.xlsx')
MNC <- read.csv('CRG/result/files/MNC_top30.csv',sep = ',')
Degree <- read.csv('CRG/result/files/Degree_top30.csv',sep = ',')
MD <- intersect(MNC$name,Degree$name)
#1
load("LUADdata/GEO/GSE19151.Rdata")
# group_list <- clindata133[colnames(exp133),c(1,4)]
# count <- round(2**exp133-1,0)
# write.table(count,'~/Desktop/test1.txt',sep = '\t',quote = F)
# write_xlsx(group_list ,'~/Desktop/test1.xlsx')
B <- read_xlsx('~/Downloads/GSE19151_degs_1.3_0.01_limma.xlsx')
#2
load("LUADdata/GEO/GSE48000.Rdata")
# group_list <- clindata133[colnames(exp133),]
# count <- round(2**exp133-1,0)
# write.table(count,'~/Desktop/test.txt',sep = '\t',quote = F)
# write_xlsx(group_list ,'~/Desktop/test.xlsx')
A <- read_xlsx('~/Downloads/GSE48000_degs_1.2_0.05_limma.xlsx')
AB <- intersect(A$Tag,B$Tag)
hubgenes <- intersect(AB,MD)
#########################hub基因正常和肿瘤boxplot
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
pd <-hubdat%>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(hubdat)[-c(1,ncol(hubdat))]),names_to = 'Gene',values_to = 'Composition')
pd$Gene <- factor(pd$Gene)
df_p_val1 <- pd %>% 
  group_by(Gene) %>% 
  wilcox_test(formula = Composition ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='Gene')

pp <- ggboxplot(pd,x = 'Gene', y='Composition',fill='group')+theme_bw()+ylab('Expression')+xlab('TCGA')+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  stat_pvalue_manual(df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(axis.text.x = element_text(size=10,face='bold',color='black'),
        axis.text.y = element_text(size=13,face='bold',color='black'),
        text=element_text(size=14,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/02_hub_TCGA.tiff',pp,width = 8,height = 6,dpi=1000)

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
pd <-hub_meth1 %>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(hub_meth1)[-c(1,ncol(hub_meth1))]),names_to = 'Gene',values_to = 'Composition')
pd$Gene <- factor(pd$Gene)
df_p_val1 <- pd %>% 
  group_by(Gene) %>% 
  wilcox_test(formula = Composition ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='Gene')

pp <- ggboxplot(pd,x = 'Gene', y='Composition',fill='group')+theme_bw()+ylab('Methylation')+xlab('TCGA')+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  stat_pvalue_manual(df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(axis.text.x = element_text(size=10,face='bold',color='black'),
        axis.text.y = element_text(size=13,face='bold',color='black'),
        text=element_text(size=14,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/02_hubMETH_TCGA.tiff',pp,width = 8,height = 6,dpi=1000)

#############TCGA cox
load( "CRG/extdata/featureselect_data.Rdata")
load('LUADdata/TCGA/clinic.Rdata')
clin_tcga <- CLINIC[,c("sampleID","gender","age_at_initial_pathologic_diagnosis","pathologic_stage","pathologic_T","pathologic_N","pathologic_M","vital_status","tobacco_smoking_history_indicator")]
colnames(clin_tcga) <- c("sampleID",'Gender','Age','Stage',"Tumor","Lymph_Node","Metastasis",'Status','Smoking')
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
clin_tumor <- clin_tumor[,-c(10:12)]
clin_tumor <- clin_tumor[clin_tumor$time!=0,]
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
covariates <- c(hubgenes[3],'Age','Stage','Tumor','Lymph_Node','Metastasis')
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


mulcox<- coxph(Surv(time, status) ~ PIK3R1+Age+Stage+Tumor+Lymph_Node, data = lxdata1);mulcox

forest <- ggforest(mulcox,  #coxph得到的Cox回归结果
                   data = lxdata1,  #数据集
                   main = 'Hazard ratio of all data',  #标题
                   cpositions = c(0.05, 0.15, 0.35),  #前三列数值的距离
                   fontsize = 1.4, #字体大小
                   refLabel = 'Reference', #相对变量的标签
                   noDigits = 3 #HR、95%CI的小数位
)
ggsave(filename="CRG/result/figs/02_PIK3R1muti_forest.tiff",forest, width=17, height=10)
modeldat <- unidata
modeldat$group <- ifelse(modeldat$PIK3R1 >median(modeldat$PIK3R1),'High PIK3R_expr','Low PIK3R_expr')
png("CRG/result/figs/02_KM_PIK3R_expr_TCGA.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(time, status) ~ group, data=modeldat),surv.median.line = "hv",palette = c('red','blue','grey'),pval = T,pval.size=10);p
dev.off()

###########################################
load("LUADdata/GEO/GSE19151.Rdata")
load("LUADdata/GEO/GSE48000.Rdata")
group <- clindata133[colnames(exp133),4]
group <- clindata133[colnames(exp133),2]
exp <- as.data.frame(t(exp133[hubgenes,]))
exp1 <- do.call(data.frame,exp)
exp1$sample <- rownames(exp)
exp1 <- cbind(group,exp1 )
rownames(exp1) <- exp1$sample
pd <-exp1 %>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(exp1)[-c(1,ncol(exp1))]),names_to = 'Gene',values_to = 'Composition')
pd$Gene <- factor(pd$Gene)
df_p_val1 <- pd %>% 
  group_by(Gene) %>% 
  wilcox_test(formula = Composition ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='Gene')

pp <- ggboxplot(pd,x = 'Gene', y='Composition',fill='group')+theme_bw()+ylab('Expression')+xlab('GSE19151')+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  stat_pvalue_manual(df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(axis.text.x = element_text(size=10,face='bold',color='black'),
        axis.text.y = element_text(size=13,face='bold',color='black'),
        text=element_text(size=14,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))
pp <- ggboxplot(pd,x = 'Gene', y='Composition',fill='group')+theme_bw()+ylab('Expression')+xlab('GSE48000')+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  stat_pvalue_manual(df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(axis.text.x = element_text(size=10,face='bold',color='black'),
        axis.text.y = element_text(size=13,face='bold',color='black'),
        text=element_text(size=14,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/02_VTEexp_GSE19151.tiff',pp,width = 8,height = 6,dpi=1000)
ggsave('CRG/result/figs/02_VTEexp_GSE48000.tiff',pp,width = 8,height = 6,dpi=1000)

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
pd <-hubdat%>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(hubdat)[-c(1,ncol(hubdat))]),names_to = 'Gene',values_to = 'Composition')
pd$Gene <- factor(pd$Gene)
df_p_val1 <- pd %>% 
  group_by(Gene) %>% 
  wilcox_test(formula = Composition ~ group_list) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='Gene')

pp <- ggboxplot(pd,x = 'Gene', y='Composition',fill='group_list')+theme_bw()+ylab('Expression')+xlab('FUSCC')+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  stat_pvalue_manual(df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(axis.text.x = element_text(size=10,face='bold',color='black'),
        axis.text.y = element_text(size=13,face='bold',color='black'),
        text=element_text(size=14,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/02_hub_FUSCC.tiff',pp,width = 8,height = 6,dpi=1000)

