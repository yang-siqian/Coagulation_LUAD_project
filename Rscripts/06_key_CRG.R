library(readxl)
library(writexl)

# Fig6A
CRRS <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7_PLA2G4B','CR2')
A <- read_xlsx('files/GSE48000_degs_1.2_0.05_limma.xlsx')
B <- read_xlsx('files/GSE19151_degs_1.2_0.05_limma.xlsx')
load('Rdata/DEG_TN_TCGA.Rdata')
C <- rownames(DEG_limma_voom)[DEG_limma_voom$adj.P.Val<0.05&abs(DEG_limma_voom$logFC)>log2(1.2)]
#C <- DEG_limma_voom[DEG_limma_voom$adj.P.Val<0.05&abs(DEG_limma_voom$logFC)>log2(1.2),]
hubgenes <-intersect(intersect(intersect(A$Tag,B$Tag),C),CRRS)
length(intersect(C,CRRS))
length(intersect(B$Tag,CRGs$gene))
length(intersect(C,CRGs$gene))
length(intersect(intersect(A$Tag,B$Tag),C))
length(intersect(intersect(B$Tag,C),CRRS))
length(intersect(intersect(A$Tag,C),CRGs$gene))
length(intersect(intersect(B$Tag,C),CRGs$gene))
length(intersect(intersect(intersect(A$Tag,B$Tag),C),CRRS))

write.xlsx(A,file = 'files/Fig6.xlsx',sheetName = "Fig6A_GSE48000_DEGs",row.names = F)
write.xlsx(B,file = 'files/Fig6.xlsx',sheetName = "Fig6A_GSE19151_DEGs",row.names = F,append=T)
write.xlsx(C,file = 'files/Fig6.xlsx',sheetName = "Fig6A_TCGA_DEGs",row.names = T,append=T)

########################################
library(grid)
library(futile.logger)
library(VennDiagram)
venn.plot <- draw.quad.venn(
  area1 = 2972,
  area2 = 9002,
  area3 = 2093,
  area4 = 11,
  n12 = 810,
  n13 = 106,
  n14 = 4,
  n23 = 534,
  n24 = 2,
  n34 = 4,
  n123 = 32,
  n124 = 2,
  n134 = 2,
  n234 = 1,
  n1234 = 1,
  category = c("GSE48000 DEGs", "GSE19151 DEGs", "TCGA DEGs", "CRRS"),
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3"),
  lwd = rep(2, 4), lty = rep("solid",4),
  cex = 2,
  fontface = 'bold',
  cat.cex = 2,
  cat.fontface = 'bold'
)
png(filename = "figures/06A_Venn_diagram_new.png",width = 3000,height = 3000,res=200)
grid.draw(venn.plot)
dev.off()

#########################hub基因正常和肿瘤boxplot
# Fig6D
library(tidyverse)
library(rstatix)
library(ggtext)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
library(tidyr)
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
exp_CRGs <- exp[rownames(exp)%in%CRGs$gene,]

group <- ifelse(as.numeric(substr(colnames(exp_CRGs),14,15))<10,'Tumor','Normal')
hubdat <- as.data.frame(t(exp_CRGs[rownames(exp_CRGs)%in%hubgenes,]))
hubdat <- do.call(data.frame,hubdat )
hubdat <- cbind(group,hubdat )
rownames(hubdat) <- colnames(exp_CRGs)
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
write.xlsx(hubdat,file = 'files/Fig6.xlsx',sheetName = 'Fig6D_exp',row.names = F,append = T)
write.xlsx(df_p_val1,file = 'files/Fig6.xlsx',sheetName = 'Fig6D_pval',row.names = F,append = T)

pp <- ggboxplot(pd,x = 'Gene', y='Composition',fill='group')+theme_bw()+ylab('Expression')+xlab('TCGA')+ylim(min(pd$Composition)-3,max(pd$Composition)+5)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  stat_pvalue_manual(size=14,df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('figures/06D_exp_TCGA.tiff',pp,width = 5,height = 5,dpi=1000)

PP <- ggboxplot(hubdat,x = 'group',y = 'CR2',fill='group')+theme_bw()+ylab('Expression')+xlab('')+ggtitle('TCGA')+ylim(min(hubdat$CR2)-2,max(hubdat$CR2)+2)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(size=1,textsize=12,comparisons = list(c('Tumor','Normal')),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))

ggsave('figures/06D_CR2_TCGA.tiff',PP,width = 5,height = 5,dpi=1000)

#######################CNV
cnv_tcga <- read.table("rawdata/TCGA_LUAD/TCGA_CNA.txt",header = T,row.names = 1,sep='\t')
cnv_tcga <- cnv_tcga[rownames(cnv_tcga)=='CR2',]
colnames(cnv_tcga) <- gsub('\\.','-',colnames(cnv_tcga))
cnv_tcga <- cnv_tcga %>% t %>% as.data.frame()
cnv_tcga$sample <- rownames(cnv_tcga)
cnv_tcga$group <- ifelse(cnv_tcga$CR2>0,'Gain',ifelse(cnv_tcga$CR2==0,'Diploid','Loss'))
cnv_tcga <- cnv_tcga[cnv_tcga$sample%in%clin_tumor$sampleID,]
table(cnv_tcga$group)
# Diploid    Gain    Loss 
# 125     352      20 
write.xlsx(cnv_tcga,file = 'files/Fig6.xlsx',sheetName = 'Fig6B_cnv',row.names = F,append = T)

############TCGA 甲基化
# Fig6C
meth <- read.table('rawdata/TCGA_LUAD/TCGA.LUAD.sampleMap_HumanMethylation450.txt',header = T,row.names = 1)
colnames(meth) <- gsub('\\.','-',colnames(meth) )
probe <- read.table('rawdata/TCGA_LUAD/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy',header = T,sep = '\t')
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
write.xlsx(hub_meth1,file = 'files/Fig6.xlsx',sheetName = 'Fig6C_methy',row.names = F,append = T)
write.xlsx(df_p_val1,file = 'files/Fig6.xlsx',sheetName = 'Fig6C_pval',row.names = F,append = T)

pp <- ggboxplot(pd,x = 'Gene', y='Composition',fill='group')+theme_bw()+ylab('Methylation')+xlab('TCGA')+ylim(min(pd$Composition)-1,max(pd$Composition)+2)+
  scale_fill_manual(values=c("orange", "darkgreen"))+
  stat_pvalue_manual(size=14,df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('figures/06C_meth_TCGA.tiff',pp,width = 5,height = 5,dpi=1000)

PP <- ggboxplot(hub_meth1,x = 'group',y = 'CR2',fill='group')+theme_bw()+ylab('Methylation')+xlab('')+ggtitle('TCGA')+ylim(min(hub_meth1$CR2)-1,max(hub_meth1$CR2)+1)+
  scale_fill_manual(values=c("orange", "darkgreen"))+
  geom_signif(size=1,textsize=12,comparisons = list(c('Tumor','Normal')),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('figures/06C_CR2meth_TCGA.tiff',PP,width = 5,height = 5,dpi=1000)

#############TCGA cox
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
exp_CRGs_tumor <- exp_CRGs[,colnames(exp_CRGs)%in%clin_tcga$sampleID[clin_tcga$group=='Tumor']] #515
exp_CRGs_tumor <- exp_CRGs_tumor %>% t%>%as.data.frame()
exp_CRGs_tumor$sampleID <- rownames(exp_CRGs_tumor)

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
clin_tumor$RFS <- clin_tumor$additional_surgery_locoregional_days
clin_tumor$RFS[which(is.na(clin_tumor$RFS))] <- clin_tumor$days_to_death[which(is.na(clin_tumor$RFS))]
clin_tumor$RFS[which(is.na(clin_tumor$RFS))] <- clin_tumor$days_to_last_followup[which(is.na(clin_tumor$RFS))]
clin_tumor$RFS <- round(clin_tumor$RFS/30,3)
clin_tumor$RFS.E <- clin_tumor$OS.E
clin_tumor$RFS.E[grep('Recurrence',clin_tumor$New_tumor_type)] <- 1
clin_tumor <- clin_tumor[which(clin_tumor$OS!=0),] #502

hubdat <-  exp_CRGs_tumor[,c('sampleID',hubgenes)]
modeldat<- merge(hubdat,clin_tumor,by.x = 'sampleID')
lxdata <- modeldat
lxdata1 <- modeldat
lxdata1$Age <- ifelse(lxdata$Age=='young',1,2 )
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$Stage <- ifelse(lxdata$Stage=='I',1,ifelse(lxdata$Stage=='II',2 ,ifelse(lxdata$Stage=='III',3,ifelse(lxdata$Stage=='IV',4,0))))
lxdata1$Stage <- factor(lxdata1$Stage,levels=c(0,1,2,3,4),labels=c('UNKOWN','I','II','III','IV'))

# Fig6H
library(survival)
library(survminer)
library(rms)
library(survivalROC)
library(dplyr)
covariates <- c(hubgenes,'Age','Stage','Tumor','Lymph_Node','Metastasis')
covariates <-fd$gene
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
mulcox<- coxph(Surv(time, status) ~ CR2+Age+Tumor+Lymph_Node+Metastasis, data = lxdata1);mulcox

forest <- ggforest(mulcox,  #coxph得到的Cox回归结果
                   data = lxdata1,  #数据集
                   main = 'Hazard ratio of all data',  #标题
                   cpositions = c(0.05, 0.15, 0.35),  #前三列数值的距离
                   fontsize = 1.4, #字体大小
                   refLabel = 'Reference', #相对变量的标签
                   noDigits = 3 #HR、95%CI的小数位
)
ggsave(filename="figures/06H_CR2muti_forest.tiff",forest, width=17, height=10)
##########################################################################
# Fig6G
lxdata1$Group <- ifelse(modeldat$CR2>median(modeldat$CR2),'High CR2_expr','Low CR2_expr')
png("figures/06G_OSKM_CR2_expr_TCGA.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(time, status) ~ Group, data=lxdata1),
                surv.median.line = "hv",size=3,
                title = "TCGA CR2",
                legend.labs = c('High CR2','Low CR2'),legend=c(0.8,0.87),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("figures/06G_RFSKM_CR2_expr_TCGA.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(RFS, RFS.E) ~ Group, data=lxdata1),size=3,
                surv.median.line = "hv",
                title = "TCGA CR2",
                legend.labs = c('High CR2','Low CR2'),legend=c(0.8,0.87),font.legend=c(20,'bold'),
                palette = c("red", "blue",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of RFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

write.xlsx(lxdata1[,c(1,21:25)],file = 'files/Fig6.xlsx',sheetName = 'Fig6G',row.names = F,append = T)

#####################################################
# Fig6E
# Female   Male 
#  270      232
# 53.78%   46.22%
pp <- ggboxplot(lxdata1,x = 'Gender',y = 'CR2',fill='Gender')+theme_bw()+ylab('Expression')+xlab('')+ggtitle('TCGA CR2')+ylim(min(lxdata1$CR2)-1,max(lxdata1$CR2)+2)+
  scale_fill_manual(values=c("#354e97", "#df5b3f"))+
  geom_signif(size=1,textsize=12,comparisons = list(c("Male", "Female")),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('figures/06E_gender_TCGA.tiff',pp,width = 5,height = 5,dpi=200)

#total 502
# UNKOWN      I     II    III     IV 
#      8    270    119     80     25 
#1.59%    53.78% 23.71% 15.94%  4.98%
pp <- ggboxplot(lxdata1[lxdata1$Stage!='UNKOWN',],x = 'Stage',y = 'CR2',fill='Stage')+theme_bw()+ylab('Expression')+xlab('')+ggtitle('TCGA CR2')+ylim(min(lxdata1$CR2)-1,max(lxdata1$CR2)+13)+
  scale_fill_manual(values=c("#354e97","70a3c4",'#f5b46f', "#df5b3f"))+
  geom_signif(size=1,textsize=10,comparisons = list(c("I", "II"),c("II", "III"),c("III",'IV'),c("II", "IV"),c("I", "III"),c("I", "IV")),y_position = c(13,14,15,18,21,24,27),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('figures/06E_Stage_TCGA.tiff',pp,width = 5,height = 5,dpi=200)
# young   old 
#   330   162
#67.07% 32.93%
pp <- ggboxplot(lxdata1[lxdata1$Age%in%c("young", "old"),],x = 'Age',y = 'CR2',fill='Age')+theme_bw()+ylab('Expression')+xlab('')+ggtitle('TCGA CR2')+ylim(min(lxdata1$CR2)-1,max(lxdata1$CR2)+3)+
  scale_fill_manual(values=c("#354e97","70a3c4"))+
  geom_signif(size=1,textsize=12,comparisons = list(c("young", "old")),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('figures/06E_Age_TCGA.tiff',pp,width = 5,height = 5,dpi=200)

write.xlsx(lxdata1[,1:5],file = 'files/Fig6.xlsx',sheetName = 'Fig6E',row.names = F,append = T)

###########################################
# Fig6F
TMB <- read_xlsx('files/02_tmb_msi.xlsx',sheet = 'TMB')
ESTIMATE <- read_xlsx('files/02_tmb_msi.xlsx',sheet = 'ESTIMATE')
TMB$group <- ifelse(TMB$Expression>median(TMB$Expression),'High CR2','Low CR2')
TMB <- TMB[TMB$SampleNameme%in%clin_tumor$sampleID,]
write.xlsx(TMB,file = 'files/Fig6.xlsx',sheetName = 'Fig6F_TMB',row.names = F,append = T)

pp <- ggboxplot(TMB,x = 'group',y = 'TMB',fill='group')+theme_bw()+ylab('TMB')+xlab('')+ggtitle('TCGA')+ylim(min(TMB$TMB)-3,max(TMB$TMB)+5)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(size=1,textsize=12,comparisons = list(c('High CR2','Low CR2')),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('Cfigures/06F_TMB_TCGA.tiff',pp,width = 5,height = 5,dpi=200)

ESTIMATE$group <- ifelse(ESTIMATE$exp>median(ESTIMATE$exp),'High CR2','Low CR2')
ESTIMATE <- ESTIMATE[ESTIMATE$SampleName%in%clin_tumor$sampleID,]
write.xlsx(ESTIMATE,file = 'files/Fig6.xlsx',sheetName = 'Fig6F_Immune',row.names = F,append = T)

imm <- ggboxplot(ESTIMATE,x = 'group',y = 'ImmuneScore',fill='group')+theme_bw()+ylab('ImmuneScore')+xlab('')+ggtitle('TCGA')+ylim(min(ESTIMATE$ImmuneScore)-1000,max(ESTIMATE$ImmuneScore)+1000)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(size=1,textsize=12,comparisons = list(c('High CR2','Low CR2')),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('figures/06F_imm_TCGA.tiff',imm,width = 5,height = 5,dpi=200)

strom <- ggboxplot(ESTIMATE,x = 'group',y = 'StromalScore',fill='group')+theme_bw()+ylab('StromalScore')+xlab('')+ggtitle('TCGA')+ylim(min(ESTIMATE$StromalScore)-1000,max(ESTIMATE$StromalScore)+1000)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(size=1,textsize=12,comparisons = list(c('High CR2','Low CR2')),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('figures/06F_strom_TCGA.tiff',strom,width = 5,height = 5,dpi=200)

est <- ggboxplot(ESTIMATE,x = 'group',y = 'ESTIMATEScore',fill='group')+theme_bw()+ylab('ESTIMATEScore')+xlab('')+ggtitle('TCGA')+ylim(min(ESTIMATE$ESTIMATEScore)-1000,max(ESTIMATE$ESTIMATEScore)+1000)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(size=1,textsize=12,comparisons = list(c('High CR2','Low CR2')),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1),
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('figures/06F_est_TCGA.tiff',est,width = 5,height = 5,dpi=200)
