# library(dplyr)
# library(writexl)
# library(KEGGREST)
# #download 209 CRGs
# ccc <- keggGet('hsa04610') #Complement and coagulation cascades
# pa <- keggGet('hsa04611')  #Platelet activation
# ccc_gs <- unlist(lapply(ccc[[1]]$GENE,function(x) strsplit(x,';')))
# ccc_gl <- ccc_gs[1:length(ccc_gs)%%3==2]
# pa_gs <- unlist(lapply(pa[[1]]$GENE,function(x) strsplit(x,';')))
# pa_gl <- pa_gs[1:length(pa_gs)%%3==2]
# genelist <- c(ccc_gl,pa_gl) %>% as.data.frame()
# write_xlsx(genelist,path = 'files/01_CRGs_209.xlsx')

#Fig1A
library(readxl)
CRGs <- read_xlsx('files/01_CRGs_209.xlsx')
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

# table(as.numeric(substr(Mut_tcga$sample,14,15))<10)       #检查全是肿瘤
# table(unique(Mut_tcga$gene)%in%CRGs$gene)                 #195/209 CRGs   195个CRG都发生mutations
# length(unique(Mut_tcga[Mut_tcga$gene%in%CRGs$gene,1]))    #450 samples
# length(unique(Mut_tcga$sample))                           #450/513 samples  87.7%

Mut_tcga_CRGs <- Mut_tcga[Mut_tcga$gene%in%CRGs$gene,]  
Mut_tcga_CRGs <- Mut_tcga_CRGs[Mut_tcga_CRGs$Tumor_Sample_Barcode%in%mut_clin$Tumor_Sample_Barcode,]
write.xlsx(Mut_tcga_CRGs,file = 'files/Fig1.xlsx',sheetName = 'Fig1A',row.names = F)
library(maftools)
mut_clin <- clin
mut_clin <- na.omit(mut_clin)
g <- read.xlsx('figures/02_data.xlsx',sheetName = 'CRRS_surv')
colnames(g)[1] <- 'Tumor_Sample_Barcode' 
colnames(mut_clin)[1] <-'Tumor_Sample_Barcode' 
mut_clin <- merge(g,mut_clin,by.x='Tumor_Sample_Barcode' )
laml = read.maf(maf = Mut_tcga_CRGs,clinicalData = mut_clin)
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
oncoplot(maf = laml,clinicalFeatures =c('group','Age','Gender','Stage') ,
         top=30,
         sortByAnnotation = TRUE, 
         bgCol='white',
         colors = vc_cols)


  
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
exp_CRGs_tumor <- exp_CRGs[,colnames(exp_CRGs)%in%clin_tcga$sampleID[clin_tcga$group=='Tumor']] #515
exp_CRGs_tumor <- exp_CRGs_tumor %>% t%>%as.data.frame()
exp_CRGs_tumor$sampleID <- rownames(exp_CRGs_tumor)

clin_tumor<- clin_tcga[clin_tcga$sampleID%in%exp_CRGs_tumor$sampleID,]
# clin_tumor$Gender <- ifelse(clin_tumor$Gender=='MALE','Male','Female')
# clin_tumor$Stage <- ifelse(clin_tumor$Stage%in%c('Stage I','Stage IA','Stage IB'),'I',
#                            ifelse(clin_tumor$Stage%in%c('Stage II','Stage IIA','Stage IIB'),"II",
#                                   ifelse(clin_tumor$Stage%in%c('Stage IIIA','Stage IIIB'),"III",
#                                          ifelse(clin_tumor$Stage%in%c('Stage IV'),'IV','UNKOWN'))))
# clin_tumor$Tumor <- ifelse(clin_tumor$Tumor%in%c('T1', 'T1a', 'T1b'),'T1',ifelse(clin_tumor$Tumor%in%c('T2', 'T2a', 'T2b'),'T2',"T3-4"))
# clin_tumor$`Lymph_Node` <- ifelse(clin_tumor$`Lymph_Node`%in%c('N0'),'N0',ifelse(clin_tumor$`Lymph_Node`%in%c('N1'),'N1',"N2-3"))
# clin_tumor$Metastasis <- ifelse(clin_tumor$Metastasis%in%c('M0'),'M0','M1')
# clin_tumor$Smoking <- ifelse(clin_tumor$Smoking%in%c(''),'UNKOWN',ifelse(clin_tumor$Smoking%in%c('Lifelong Non-smoker'),'Never smoker','Current smoker'))
# clin_tumor$Status <- ifelse(clin_tumor$Status%in%c('LIVING'),'Alive','Dead')
# clin_tumor$Age <- ifelse(clin_tumor$Age>70,'old','young')
clin_tumor$RFS <- clin_tumor$additional_surgery_locoregional_days
clin_tumor$RFS[which(is.na(clin_tumor$RFS))] <- clin_tumor$days_to_death[which(is.na(clin_tumor$RFS))]
clin_tumor$RFS[which(is.na(clin_tumor$RFS))] <- clin_tumor$days_to_last_followup[which(is.na(clin_tumor$RFS))]
clin_tumor$RFS <- round(clin_tumor$RFS/30,3)
clin_tumor$RFS.E <- clin_tumor$OS.E
clin_tumor$RFS.E[grep('Recurrence',clin_tumor$New_tumor_type)] <- 1
clin_tumor <- clin_tumor[which(clin_tumor$OS!=0),] #502

# SupFig1C
#mut and non-mut OS AND RFS  405|97
survdata<- clin_tumor[,c(1,20,21,22,23)]
survdata$group <- ifelse(survdata$sampleID%in%unique(Mut_tcga_CRGs$sample[Mut_tcga_CRGs$effect%in%c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation')]),'Mutation','Non-mutation')
library(xlsx)
write.xlsx(survdata,file='files/SupFig1.xlsx',sheetName = 'SupFig1C',row.names = F)
png("figures/S01C_OSKM_mut_TCGA.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(OS, OS.E) ~ group, data=survdata),
                surv.median.line = "hv",size=3,
                title = "TCGA",
                legend.labs = c('Mutation','Non-mutation'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("figures/S01C_RFSKM_mut_TCGA.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(RFS, RFS.E) ~ group, data=survdata),size=3,
                surv.median.line = "hv",
                title = "TCGA",
                legend.labs = c('Mutation','Non-mutation'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("red", "blue",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of RFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

######################CNV
cnv_tcga <- read.table("rawdata/TCGA_LUAD/TCGA_CNA.txt",header = T,row.names = 1,sep='\t')
#table(as.numeric(substr(colnames(cnv_tcga),14,15))<10)   # 516 检查全是肿瘤
#table(rownames(cnv_tcga)%in%CRGs$gene)                   # 199/209 CRGs
cnv_tcga <- cnv_tcga[rownames(cnv_tcga)%in%CRGs$gene,]
colnames(cnv_tcga) <- gsub('\\.','-',colnames(cnv_tcga))

# SupFig1B
#CNA(2和-2) AND NON-CNA OS AND RFS  358|139  
x <- c()
for(i in 1:ncol(cnv_tcga)){
  if(length(rownames(cnv_tcga)[abs(cnv_tcga[,i])==2])==0){x <- c(x,i)}
}
survdata1<- clin_tumor[,c(1,20,21,22,23)]
survdata1$group <- ifelse(survdata$sampleID%in%colnames(cnv_tcga)[x],'Non-CNA',ifelse(survdata$sampleID%in%colnames(cnv_tcga)[-x],'CNA','NA'))
survdata1 <- survdata1[survdata1$group!='NA',]
library(xlsx)
write.xlsx(survdata1,file='files/SupFig1.xlsx',sheetName = 'SupFig1B',row.names = F,append = T)
png("figures/S01B_OSKM_CNA_TCGA.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(OS, OS.E) ~ group, data=survdata1),
                surv.median.line = "hv",size=3,
                title = "TCGA",
                legend.labs = c('CNA','Non-CNA'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("figures/S01B_RFSKM_CNA_TCGA.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(RFS, RFS.E) ~ group, data=survdata1),size=3,
                surv.median.line = "hv",
                title = "TCGA",
                legend.labs = c('CNA','Non-CNA'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("red", "blue",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of RFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()
############################AMP ADN DEL
# Fig1D
#AMP AND DEL OS AND RFS  271|77  
a <- c()
d <- c()
for(i in 1:ncol(cnv_tcga)){
  if(sum(cnv_tcga[,i][abs(cnv_tcga[,i])==2])<0){d <- c(d,i)}
  if(sum(cnv_tcga[,i][abs(cnv_tcga[,i])==2])>0){a <- c(a,i)}
}
survdata2<- clin_tumor[,c(1,20,21,22,23)]
survdata2$group <- ifelse(survdata2$sampleID%in%colnames(cnv_tcga)[a],'Amplification',ifelse(survdata2$sampleID%in%colnames(cnv_tcga)[d],'Deep deletion','NA'))
survdata2 <- survdata2[survdata2$group!='NA',]
library(xlsx)
write.xlsx(survdata2,file='files/Fig1.xlsx',sheetName = 'Fig1D',row.names = F,append = T)

png("figures/01D_OSKM_ampdel_TCGA.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(OS, OS.E) ~ group, data=survdata2),
                surv.median.line = "hv",size=3,
                title = "TCGA",
                legend.labs = c('Amplification','Deep deletion'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("figures/01_RFSKM_ampdel_TCGA.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(RFS, RFS.E) ~ group, data=survdata2),size=3,
                surv.median.line = "hv",
                title = "TCGA",
                legend.labs = c('Amplification','Deep deletion'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("red", "blue",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of RFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

# Fig1F
freq <- data.frame(matrix(ncol = 3))
colnames(freq) <- c('gene','Amplification','Deletion')
for(i in rownames(cnv_tcga)){
  a <- c(i,length(cnv_tcga[i,][,cnv_tcga[i,]==2])/ncol(cnv_tcga),length(cnv_tcga[i,][,cnv_tcga[i,]==-2])/ncol(cnv_tcga))
  freq <- rbind(freq,a)
}
freq <- na.omit(freq)
freq <- freq[order(freq$Amplification,decreasing = T),]
fd <- freq[1:30,]
fd$gene <- factor(fd$gene,levels=freq$gene)
fd$Amplification <- round(as.numeric(fd$Amplification)*100)
fd$Deletion <- round(as.numeric(fd$Deletion)*100)
library(xlsx)
write.xlsx(fd,file='files/Fig1.xlsx',sheetName = 'Fig1F',row.names = F,append = T)
library(ggalt)
cp <- ggplot(fd,aes(y = gene, x=Amplification,xend=Deletion))+
  geom_segment(aes(x=Amplification,xend=Deletion,y=gene,yend=gene),color='grey',size=2)+
  geom_dumbbell(size_x=5, size_xend = 5,color="grey",colour_x = "red", colour_xend = "blue")+
  geom_text(size=4,aes(x=Amplification,label=Amplification,vjust=-1.1))+
  geom_text(size=4,aes(x=Deletion,label=Deletion,vjust=1.8))+
  labs(x='CNA frequency%', y=NULL) + scale_x_continuous(breaks=seq(0,100,10))+
  theme_classic() + theme(axis.text.x = element_text(size=16,face='bold',color='black',angle = 90,hjust = 0.8,vjust = 0.9),
                          axis.text.y = element_text(size=10,face='bold',color='black'),
                          text=element_text(size=16,face='bold'),axis.line = element_line(size=1))+coord_flip()

ggsave('figures/01_dot_CNA.tiff',cp,width = 10,height = 8,dpi=1000)

# SupFig1D
expdata<- exp_CRGs_tumor[,colnames(exp_CRGs_tumor)%in%fd$gene[1:6]]
a <- c()
d <- c()
for(i in 1:ncol(cnv_tcga)){
  if(sum(cnv_tcga[,i][abs(cnv_tcga[,i])==2])<0){d <- c(d,i)}
  if(sum(cnv_tcga[,i][abs(cnv_tcga[,i])==2])>0){a <- c(a,i)}
}

for(gene in fd$gene[1:6]){
  group <- paste0('group_',gene)
  expdata$group <- ifelse(rownames(expdata)%in%colnames(cnv_tcga)[cnv_tcga[gene,]==2],'Amplification',ifelse(rownames(expdata)%in%colnames(cnv_tcga)[cnv_tcga[gene,]==0],'Diploid',''))
  colnames(expdata)[ncol(expdata)] <- group
  expdata1 <- expdata[,which(colnames(expdata)%in%c(gene,group))]
  expdata1 <- expdata1[expdata1[,2]!='',]
  cb <- ggboxplot(expdata1,x = group,y = gene,fill=group)+theme_bw()+ylab('Expression')+xlab('')+ggtitle(gene)+ylim(min(expdata1[,gene])-3,max(expdata1[,gene])+3)+
    scale_fill_manual(values=c("#70a3c4","#ed6f59"))+
    geom_signif(size=1,textsize=12,comparisons = list(c('Amplification','Diploid')),test = wilcox.test,map_signif_level = T)+
    theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
          axis.text.x = element_text(size=18,face='bold',color='black'),
          axis.text.y = element_text(size=18,face='bold',color='black'),
          text=element_text(size=25,face='bold'),
          legend.title = element_blank(),legend.position = "none",
          panel.grid = element_blank(),panel.border = element_rect(size=1.5))
  ggsave(paste0('figures/S01D_cnvbox_',gene,'.png'),cb,width = 5,height = 5,dpi=200)
}
write.xlsx(expdata,file='files/SupFig1.xlsx',sheetName = 'SupFig1D',append = T)

# SupFig1A
M <- Mut_tcga[Mut_tcga$effect%in%c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation'),]
for(gene in c('COL3A1','F8','ITGAX','ADCY2','PLCB1','ADCY8')){
  p=length(unique(M$sample[M$gene==gene]))/length(unique(M$sample))
  print(paste0(gene,':',round(p*100,1),'%'))
}
# "COL3A1:10.7%"
# "F8:8.4%"
# "ITGAX:8%"
# "ADCY2:7.8%"
# "PLCB1:7.8%"
# "ADCY8:6.8%"


# SupFig1E
expdata2<- exp_CRGs_tumor[,colnames(exp_CRGs_tumor)%in%c('COL3A1','F8','ITGAX','ADCY2','PLCB1','ADCY8')]
for(gene in c('COL3A1','F8','ITGAX','ADCY2','PLCB1','ADCY8')){
  group <- paste0('group_',gene)
  expdata2$group <- ifelse(rownames(expdata2)%in%unique(M$sample[M$gene==gene]),'Mutation','Non-mutation')
  colnames(expdata2)[ncol(expdata2)] <- group
  expdata3 <- expdata2[,which(colnames(expdata2)%in%c(gene,group))]
  expdata3 <- expdata3[expdata3[,2]!='',]
  cb <- ggboxplot(expdata3,x = group,y = gene,fill=group)+theme_bw()+ylab('Expression')+xlab('')+ggtitle(gene)+ylim(min(expdata3[,gene])-3,max(expdata3[,gene])+3)+
    scale_fill_manual(values=c("#ebd657", "#334c81"))+
    geom_signif(size=1,textsize=12,comparisons = list(c('Mutation','Non-mutation')),test = wilcox.test,map_signif_level = T)+
    theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
          axis.text.x = element_text(size=18,face='bold',color='black'),
          axis.text.y = element_text(size=18,face='bold',color='black'),
          text=element_text(size=25,face='bold'),
          legend.title = element_blank(),legend.position = "none",
          panel.grid = element_blank(),panel.border = element_rect(size=1.5))
  ggsave(paste0('figures/S01E_mutbox_',gene,'.png'),cb,width = 5,height = 5,dpi=200)

}
write.xlsx(expdata2,file='files/SupFig1.xlsx',sheetName = 'SupFig1E',append = T)

