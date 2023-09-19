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

bsrgenes <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','ITGB1','C4BPA','JMJD7_PLA2G4B','CR2')
bsrgenes[10] <- gsub('_','-',bsrgenes[10])
hubdat <-  exp_CRGs_tumor[,c('sampleID',bsrgenes)]
modeldat<- merge(hubdat,clin_tumor,by.x = 'sampleID')
lxdata <- modeldat
lxdata1 <- merge(exp_CRGs_tumor,clin_tumor,by.x = 'sampleID')
lxdata1$OS <- round(lxdata1$OS,3)
lxdata1$Age <- ifelse(lxdata$Age=='young',1,2 )
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$Stage <- ifelse(lxdata$Stage=='I',1,ifelse(lxdata$Stage=='II',2 ,ifelse(lxdata$Stage=='III',3,ifelse(lxdata$Stage=='IV',4,0))))
lxdata1$Stage <- factor(lxdata1$Stage,levels=c(0,1,2,3,4),labels=c('UNKOWN','I','II','III','IV'))
colnames(lxdata1) <- gsub('-','_',colnames(lxdata1))
#mulcox<- coxph(Surv(OS, OS.E) ~ CR2+F2+F12+GUCY1A1+ITGB1+P2RX1+SERPINE1+PLAUR+PRKCZ+C4BPA+JMJD7_PLA2G4B, data = lxdata1 );mulcox
lxdata1$Risk=predict(mulcox,lxdata1)
lxdata1$Group <- ifelse(lxdata1$Risk>0.192116,'High risk','Low risk')

#############gsva
library(readxl)
enrich_score <- read.table('files/enrichScore.txt',header = T,row.names = 1)
group <- read_xlsx('files/Fig3.xlsx',sheet = 'Fig3B')
DEG <- read_xlsx('files/Fig4.xlsx',sheet = 'Fig4A')
colnames(DEG)[1] <- 'Names'
DEG_kegg <- DEG[DEG$P.Value<0.05,]

############gsea
groupFile <- as.data.frame(lxdata1[,c(1,210)])
write.xlsx(groupFile,file = 'files/GSEA_groupfile.xlsx')
load('Rdata/DEG_TN_TCGA.Rdata')
gseaexp <- exp[rownames(exp)%in%rownames(DEG_limma_voom)[DEG_limma_voom$adj.P.Val<0.05&abs(DEG_limma_voom$logFC)>log2(1)],]
write.table(gesaexp ,'files/GSEA_exp.txt',sep = '\t',quote = F)
##############TMB
# Fig4C_TMB
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
library(xlsx)
write.xlsx(tmb,file = 'files/Fig4.xlsx',sheetName = 'Fig4C_TMB',row.names = F,append = T)
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
ggsave('figures/04C_TMB.png',tb,width = 5,height =5,dpi=200)

##############Estimate


ssGSEA <- function(exp,lxdata1){
  library(estimate)
  library(GSVA)
  library(GSEABase)
  library(limma)
  library(readxl)
  gene_set<-read_xlsx('files/immueGenelist.xlsx')
  gene_set <- gene_set[,c(1,2)]
  colnames(gene_set) <- c('gene','type') 
  list<- split(as.matrix(gene_set)[,1], gene_set[,2])
  
  normalization<-function(x){
    return((x-min(x))/(max(x)-min(x)))}
  
  for(i in c('High risk','Low risk')){
    immueExpr <- exp[rownames(exp)%in%gene_set$gene,colnames(exp)%in%lxdata1$sampleID[lxdata1$Group==i]]
    ssGSEA_Score <- gsva(as.matrix(immueExpr),list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
    norm_ssGSEA_Score <- normalization(ssGSEA_Score)%>%t%>%as.data.frame()
    norm_ssGSEA_Score$group <- rep(i,nrow(norm_ssGSEA_Score))
    if(i=='High risk'){
      norm_ssGSEA_Score1 <- norm_ssGSEA_Score
    }else if(i=='Low risk'){
      norm_ssGSEA_Score2 <- norm_ssGSEA_Score
    }
  }
  norm_ssGSEA <- rbind(norm_ssGSEA_Score1,norm_ssGSEA_Score2)
  colnames(norm_ssGSEA ) <- c('Act_B','Act_CD4','Act_CD8','Act_DC','CD56bright',
                                    'CD56dim','Tcm_CD4','Tcm_CD8','Tem_CD4','Tem_CD8',
                                    "Eosinophil",'Tgd','Imm_B','iDC',"Macrophage" ,
                                    'Mast',"MDSC",'Mem_B',"Monocyte",'NK',
                                    'NKT',"Neutrophil" ,'pDC','Treg','Tfh',
                                    'Th1','Th17','Th2',"group" )
  
  return(norm_ssGSEA)
}

norm_ssGSEA_Score <- ssGSEA(exp,lxdata1)

box_plot <- function(score){
  library(tidyverse)
  library(rstatix)
  library(ggtext)
  library(ggstatsplot)
  library(ggpubr)
  library(tidyr)
  fd <-score %>% 
    rownames_to_column('sample')%>%
    pivot_longer(cols=c(colnames(score)[-ncol(score)]),names_to = 'Celltype',values_to = 'Composition')
  fd$Celltype <- factor(fd$Celltype)
  df_p_val1 <- fd %>% 
    group_by(Celltype) %>% 
    wilcox_test(formula = Composition ~ group) %>% 
    add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
    add_xy_position(x='Celltype')
  
  # library(writexl)
  # write_xlsx(df_p_val1,path = 'figures/05_data.xlsx',sheetName = 'ssGSEA_pval',row.names = F,append = T)
  # 
  # Fig4D
  p<- ggboxplot(fd,x = 'Celltype', y='Composition',fill='group')+theme_bw()+ylab('')+xlab('')+
    scale_fill_manual(values=c("#DC143C", "#0000FF"))+
    stat_pvalue_manual(size=8,remove.bracket=T,
                       y.position = 1.2,
                       vjust=1,df_p_val1,label = '{p.signif}',tip.length = 0)+
    theme(axis.text.x = element_text(size=18,face='bold',color='black',angle = 90,hjust = 0.8,vjust = 0.9),
          axis.text.y = element_text(size=18,face='bold',color='black'),
          text=element_text(size=18,face='bold'),
          legend.title = element_blank(),legend.position = "top",
          panel.grid = element_blank(),panel.border = element_rect(size=2))
  
  return(p)
  
}

im <- box_plot(norm_ssGSEA_Score)
im
write.xlsx(norm_ssGSEA_Score,file = 'figures/05_data.xlsx',sheetName = 'norm_ssGSEA_Score',row.names = T,append = T)
ggsave('figures/05D_immue.png',im,width = 15,height = 10,dpi=200)


#############################################################################################################
immue_scores <- function(exp,lxdata1){
  gene_set<-read_xlsx('files/immueGenelist.xlsx')
  gene_set <- gene_set[,c(1,2)]
  colnames(gene_set) <- c('gene','type')
  for(i in 1:2){
    if(i==1){
      exp_cluster <- exp[rownames(exp)%in%gene_set$gene,colnames(exp)%in%lxdata1$sampleID[lxdata1$Group=='High risk']]
    }else{
      exp_cluster <- exp[rownames(exp)%in%gene_set$gene,colnames(exp)%in%lxdata1$sampleID[lxdata1$Group=='Low risk']]
    }
    write.table(exp_cluster ,file=paste0("files/exp_cluster",i,".txt"),sep = '\t',quote = F)
    filterCommonGenes(input.f=paste0("files/exp_cluster",i,".txt"),
                      output.f=paste0("files/exp_cluster",i,".gct"),
                      id="GeneSymbol")
    
    estimateScore(input.ds = paste0("files/exp_cluster",i,".gct"), 
                  output.ds=paste0("files/estimateScore_cluster",i,".gct"),
                  platform="illumina")
    scores <- read.table(paste0("files/estimateScore_cluster",i,".gct"),skip = 2,header = T)
    scores <- scores %>% t %>% as.data.frame()
    colnames(scores) <- scores[1,]
    scores <- scores[-c(1,2),]
    scores$sample <- gsub('\\.','-',rownames(scores))
    if(i==1){
      scores$group <- rep('High risk',nrow(scores))
      scores1 <- scores
    }else{
      scores$group <- rep('Low risk',nrow(scores))
      scores2 <- scores
    }
  }
  scores <- rbind(scores1,scores2)
  scores$StromalScore <- as.numeric(scores$StromalScore)
  scores$ImmuneScore <- as.numeric(scores$ImmuneScore)
  scores$ESTIMATEScore <- as.numeric(scores$ESTIMATEScore)
  return(scores)
}

scores <- immue_scores(exp,lxdata1)
rownames(scores) <- scores$sample
scores <- merge(TMB,scores,by.x='sample')
write.xlsx(scores,file='figures/05_data.xlsx',sheetName = 'estimate_tmb',row.names = T,append = T)

estimate_plot <- function(scores){
  par(mfrow=c(2,2))
  tb <- ggboxplot(scores,x="group",y="TMB",fill = "group")+theme_classic()+ylab('TMB')+xlab('')+ggtitle('')+ylim(0,40)+
    geom_signif(size=1,textsize=12,comparisons = list(c('High risk','Low risk')),test = wilcox.test,map_signif_level = T)+
    theme( plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
           axis.text.x = element_text(size=18,face='bold',color='black'),
           axis.text.y = element_text(size=18,face='bold',color='black'),
           text=element_text(size=25,face='bold'),
           legend.title = element_blank(),legend.position = "none",
           axis.line = element_line(size=1))
  # Fig4C_StromalScore
  ss <- ggboxplot(scores,x="group",y="StromalScore",fill = "group")+theme_classic()+ylab('StromalScore')+xlab('')+ggtitle('')+ylim(-100,130)+
    geom_signif(size=1,textsize=12,comparisons = list(c('High risk','Low risk')),test = wilcox.test,map_signif_level = T)+
    theme( plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
           axis.text.x = element_text(size=18,face='bold',color='black'),
           axis.text.y = element_text(size=18,face='bold',color='black'),
           text=element_text(size=25,face='bold'),
           legend.title = element_blank(),legend.position = "none",
           axis.line = element_line(size=1))
  # Fig4C_ImmuneScore
  is <- ggboxplot(scores,x='group',y='ImmuneScore',fill = "group")+theme_classic()+ylab('ImmuneScore')+xlab('')+ggtitle('')+ylim(-60,140)+
    geom_signif(size=1,textsize=12,comparisons = list(c('High risk','Low risk')),test = wilcox.test,map_signif_level = T)+
    theme( plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
           axis.text.x = element_text(size=18,face='bold',color='black'),
           axis.text.y = element_text(size=18,face='bold',color='black'),
           text=element_text(size=25,face='bold'),
           legend.title = element_blank(),legend.position = "none",
           axis.line = element_line(size=1))
  # Fig4C_ESTIMATEScore
  es <- ggboxplot(scores,x='group',y='ESTIMATEScore',fill = "group")+theme_classic()+ylab('ESTIMATEScore')+xlab('')+ggtitle('')+ylim(-150,220)+
    geom_signif(size=1,textsize=12,comparisons = list(c('High risk','Low risk')),test = wilcox.test,map_signif_level = T)+
    theme( plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
           axis.text.x = element_text(size=18,face='bold',color='black'),
           axis.text.y = element_text(size=18,face='bold',color='black'),
           text=element_text(size=25,face='bold'),
           legend.title = element_blank(),legend.position = "none",
           axis.line = element_line(size=1))
  
  p=tb+ss+is+es
  return(p)
}
p <- estimate_plot(scores)
p
ggsave('figures/05_estimate.png',p,width = 10,height = 10,dpi=200)

##############################MHC
# Fig4E
YZ_score <- function(exp,lxdata1,X){
  if(X=='MHC'){
    genes <- c("HLA-F","HLA-E","HLA-G","HLA-DRB5","HLA-DRB1",
             "HLA-DRA","HLA-DQB2","HLA-DQB1","HLA-DQA2","HLA-DQA1",
             "HLA-DPB1","HLA-DPA1","HLA-DOB","HLA-DOA","HLA-DMB","HLA-DMA","HLA-C","HLA-B","HLA-A")
  }else if(X=='TCAF'){
    genes <- c('TNFRSF8','TNFRSF25','TNFRSF17','TNFRSF14',
              'ICOS','CD40LG','CD28','CD27','CD226','CD2')
  }
  
  c1 <- exp[rownames(exp)%in%genes,colnames(exp)%in%lxdata1$sampleID[lxdata1$Group=='High risk']]
  c2 <- exp[rownames(exp)%in%genes,colnames(exp)%in%lxdata1$sampleID[lxdata1$Group=='Low risk']]
  c1 <- c1 %>% t %>% as.data.frame()
  c1$sample <- rownames(c1)
  c1$group <- rep("High risk",nrow(c1))
  c2 <- c2 %>% t %>% as.data.frame()
  c2$sample <- rownames(c2)
  c2$group <- rep("Low risk",nrow(c2))
  all <- rbind(c1,c2)
  return(all)
}

MHC_score <- YZ_score(exp,lxdata1,'MHC')
TCAF_score <- YZ_score(exp,lxdata1,'TCAF')
write.xlsx(MHC_score,file = 'figures/05_data.xlsx',sheetName = 'MHC',row.names = T,append = T)
#write_xlsx(df_p_val1,path = 'figures/05_data.xlsx',sheetName = 'MHC_pval',row.names = F,append = T)
write.xlsx(TCAF_score,file = 'figures/05_data.xlsx',sheetName = 'TCAF',row.names = T,append = T)
#write_xlsx(df_p_val1,path = 'figures/05_data.xlsx',sheetName = 'TCAF_pval',row.names = F,append = T)
YZ_plot <- function(score1,score2){
  library(tidyr)
  par(mfrow=c(2,1))
  pd <-score1 %>% 
    rownames_to_column('Sample')%>%
    pivot_longer(cols=c(colnames(score1)[-c(ncol(score1)-1,ncol(score1))]),names_to = 'Celltype',values_to = 'Composition')
  pd $Celltype <- factor(pd$Celltype)
  df_p_val <- pd  %>% 
    group_by(Celltype) %>% 
    wilcox_test(formula = Composition ~ group) %>% 
    add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
    add_xy_position(x='Celltype')
  
  p1<- ggboxplot(pd ,x = 'Celltype', y='Composition',fill='group')+theme_bw()+ylab('')+xlab('')+
    scale_fill_manual(values=c("#DC143C", "#0000FF"))+
    stat_pvalue_manual(size=14,remove.bracket=T,y.position = 25,vjust=1,df_p_val,label = '{p.signif}',tip.length = 0)+
    theme(axis.text.x = element_text(size=18,face='bold',color='black'),
          axis.text.y = element_text(size=18,face='bold',color='black'),
          text=element_text(size=18,face='bold'),
          legend.title = element_blank(),legend.position = "top",
          panel.grid = element_blank(),panel.border = element_rect(size=2))+coord_flip()
  td <-score2 %>% 
    rownames_to_column('Sample')%>%
    pivot_longer(cols=c(colnames(score2)[-c(ncol(score2)-1,ncol(score2))]),names_to = 'Celltype',values_to = 'Composition')
  td$Celltype <- factor(td$Celltype)
  df_p_val1 <- td  %>% 
    group_by(Celltype) %>% 
    wilcox_test(formula = Composition ~ group) %>% 
    add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
    add_xy_position(x='Celltype')
  
  p2<- ggboxplot(td ,x = 'Celltype', y='Composition',fill='group')+theme_bw()+ylab('')+xlab('')+
    scale_fill_manual(values=c("#DC143C", "#0000FF"))+
    stat_pvalue_manual(size=14,remove.bracket=T,y.position = 15,vjust=1,df_p_val1,label = '{p.signif}',tip.length = 0)+
    theme(axis.text.x = element_text(size=18,face='bold',color='black'),
          axis.text.y = element_text(size=18,face='bold',color='black'),
          text=element_text(size=18,face='bold'),
          legend.title = element_blank(),legend.position = "top",
          panel.grid = element_blank(),panel.border = element_rect(size=2))+coord_flip()
  
  p=p1+p2
  return(p)
  
}
p <- YZ_plot(MHC_score,TCAF_score)
p
ggsave('figures/05_YZ.png',p,width = 18,height = 12,dpi=200)




