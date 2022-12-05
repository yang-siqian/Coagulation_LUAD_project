library(readxl)
library(maftools)
CRGs <- read_xlsx('CRG/result/files/CRGs_209.xlsx')
MNC <- read.csv('CRG/result/files/MNC20.csv',sep = ',')
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
table(unique(Mut_tcga$gene)%in%CRGs$gene) #195/209 CRGs 450/513 samples  87.7%
Mut_tcga <- Mut_tcga[Mut_tcga$gene%in%CRGs$gene,]
table(as.numeric(substr(Mut_tcga$sample,14,15))<10)   # 检查全是肿瘤
length(unique(Mut_tcga$sample))  # mutation 513 samples
laml = read.maf(maf = Mut_tcga)
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
oncoplot(maf = laml,
         top=30,
         sortByAnnotation = TRUE, 
         bgCol='white',
         colors = vc_cols)

M <- Mut_tcga[Mut_tcga$effect%in%c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation'),]
length(unique(M$sample[M$gene=='COL3A1']))/511   #10.76%
length(unique(M$sample[M$gene=='F8']))/511   #8.41%
length(unique(M$sample[M$gene=='ITGAX']))/511   #8.02%
length(unique(M$sample[M$gene=='ADCY2']))/511   #7.83%
length(unique(M$sample[M$gene=='PLCB1']))/511   #7.83%
length(unique(M$sample[M$gene=='ADCY8']))/511   #6.85%

expdata<- exp_CRGs_tumor[,colnames(exp_CRGs_tumor)%in%c('COL3A1','F8','ITGAX','ADCY2','PLCB1','ADCY8')]

for(gene in c('COL3A1','F8','ITGAX','ADCY2','PLCB1','ADCY8')){
  gene='ADCY8'
  expdata$group <- ifelse(rownames(expdata)%in%unique(M$sample[M$gene==gene]),'Mutation','Non-mutation')
  
  cb <- ggboxplot(expdata,x = 'group',y = gene,fill='group')+theme_bw()+ylab('Expression')+xlab('')+ggtitle(gene)+ylim(min(expdata[,gene])-3,max(expdata[,gene])+3)+
    scale_fill_manual(values=c("#ebd657", "#334c81"))+
    geom_signif(size=1,textsize=12,comparisons = list(c('Mutation','Non-mutation')),test = wilcox.test,map_signif_level = T)+
    theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
          axis.text.x = element_text(size=18,face='bold',color='black'),
          axis.text.y = element_text(size=18,face='bold',color='black'),
          text=element_text(size=25,face='bold'),
          legend.title = element_blank(),legend.position = "none",
          panel.grid = element_blank(),panel.border = element_rect(size=1.5))
  ggsave(paste0('CRG/result/figs/01_mutbox_',gene,'.png'),cb,width = 5,height = 5,dpi=200)
  
}

###################################################################
#mut and non-mut OS AND RFS   106/405
library(survival)
library(survminer)
length(unique(Mut_tcga$sample))
survdata<- clin_tumor[,c(1,20,21,23,24)]
survdata$group <- ifelse(survdata$sampleID%in%unique(Mut_tcga$sample[Mut_tcga$effect%in%c('Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation')]),'Mutation','Non-mutation')
png("CRG/result/figs/01_OSKM_mut_TCGA.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(time, status) ~ group, data=survdata),
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

png("CRG/result/figs/01_RFSKM_mut_TCGA.png",width = 1500,height = 1500,res=200)
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


###################################################################
cnv_tcga <- read.table("LUADdata/TCGA/TCGA_CNA.txt",header = T,row.names = 1,sep='\t')
table(as.numeric(substr(colnames(cnv_tcga),14,15))<10)   # 516 检查全是肿瘤
table(rownames(cnv_tcga)%in%CRGs$gene) #199/209 CRGs
cnv_tcga <- cnv_tcga[rownames(cnv_tcga)%in%CRGs$gene,]
colnames(cnv_tcga) <- gsub('\\.','-',colnames(cnv_tcga))

table(colSums(abs(cnv_tcga))==0) #检查有6人没有CRGs的SCNA

##############################################################
#CNA(2和-2) AND NON-CNA OS AND  139/358   72%
x <- c()
for(i in 1:ncol(cnv_tcga)){
  if(length(rownames(cnv_tcga)[abs(cnv_tcga[,i])==2])==0){x <- c(x,i)}
}

survdata$Group <- ifelse(survdata$sampleID%in%colnames(cnv_tcga)[x],'Non-CNA',ifelse(survdata$sampleID%in%colnames(cnv_tcga)[-x],'CNA','NA'))

png("CRG/result/figs/01_OSKM_CNA_TCGA.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(time, status) ~ Group, data=survdata[survdata$Group!='NA',]),
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

png("CRG/result/figs/01_RFSKM_CNA_TCGA.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(RFS, RFS.E) ~ Group, data=survdata[survdata$Group!='NA',]),size=3,
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

#############################################################
a <- c()
d <- c()
for(i in 1:ncol(cnv_tcga)){
  if(sum(cnv_tcga[,i][abs(cnv_tcga[,i])==2])<0){d <- c(d,i)}
  if(sum(cnv_tcga[,i][abs(cnv_tcga[,i])==2])>0){a <- c(a,i)}
}

survdata$Group1 <- ifelse(survdata$sampleID%in%colnames(cnv_tcga)[a],'Amplification',ifelse(survdata$sampleID%in%colnames(cnv_tcga)[d],'Deep deletion','NA'))
#a:d=281:81
#56.5%:27.3%
#amp:del=271:77

png("CRG/result/figs/01_OSKM_ampdel_TCGA.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(time, status) ~ Group1, data=survdata[survdata$Group1!='NA',]),
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

png("CRG/result/figs/01_RFSKM_ampdel_TCGA.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(RFS, RFS.E) ~ Group1, data=survdata[survdata$Group1!='NA',]),size=3,
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


##############################################################

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

ggsave('CRG/result/figs/01_dot_CNA.tiff',cp,width = 10,height = 8,dpi=1000)

############################################################
expdata<- exp_CRGs_tumor[,colnames(exp_CRGs_tumor)%in%fd$gene]
a <- c()
d <- c()
for(i in 1:ncol(cnv_tcga)){
  if(sum(cnv_tcga[,i][abs(cnv_tcga[,i])==2])<0){d <- c(d,i)}
  if(sum(cnv_tcga[,i][abs(cnv_tcga[,i])==2])>0){a <- c(a,i)}
}


for(gene in fd$gene){
  # gene='C6'
  expdata$group <- ifelse(rownames(expdata)%in%colnames(cnv_tcga)[cnv_tcga[gene,]==2],'Amplification',ifelse(rownames(expdata)%in%colnames(cnv_tcga)[cnv_tcga[gene,]==0],'Diploid',''))
  expdata1 <- expdata[expdata$group!='',]
  cb <- ggboxplot(expdata1,x = 'group',y = gene,fill='group')+theme_bw()+ylab('Expression')+xlab('')+ggtitle(gene)+ylim(min(expdata[,gene])-3,max(expdata[,gene])+3)+
    scale_fill_manual(values=c("#70a3c4","#ed6f59"))+
    geom_signif(size=1,textsize=12,comparisons = list(c('Amplification','Diploid')),test = wilcox.test,map_signif_level = T)+
    theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
          axis.text.x = element_text(size=18,face='bold',color='black'),
          axis.text.y = element_text(size=18,face='bold',color='black'),
          text=element_text(size=25,face='bold'),
          legend.title = element_blank(),legend.position = "none",
          panel.grid = element_blank(),panel.border = element_rect(size=1.5))
  ggsave(paste0('CRG/result/figs/01_cnvbox_',gene,'.png'),cb,width = 5,height = 5,dpi=200)
  
}

expdata<- exp_CRGs_tumor[,colnames(exp_CRGs_tumor)%in%c('RASGRP2','CR2','THBD')]
expdata$group <- ifelse(rownames(expdata)%in%colnames(cnv_tcga)[cnv_tcga['RASGRP2',]>0],'Gain',ifelse(rownames(expdata)%in%colnames(cnv_tcga)[cnv_tcga[gene,]==0],'Diploid','Loss'))
table(expdata$group)
# Diploid    Gain    Loss 
#  151     133     231 
expdata$group <- ifelse(rownames(expdata)%in%colnames(cnv_tcga)[cnv_tcga['CR2',]>0],'Gain',ifelse(rownames(expdata)%in%colnames(cnv_tcga)[cnv_tcga[gene,]==0],'Diploid','Loss'))
table(expdata$group)
# Diploid    Gain    Loss 
#    75     359      81 
expdata$group <- ifelse(rownames(expdata)%in%colnames(cnv_tcga)[cnv_tcga['THBD',]>0],'Gain',ifelse(rownames(expdata)%in%colnames(cnv_tcga)[cnv_tcga[gene,]==0],'Diploid','Loss'))
table(expdata$group)
# Diploid    Gain    Loss 
#    219     195     101