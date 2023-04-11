library(survival)
library(survminer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(readxl)

exp <- read.table('rawdata/Melanoma/TCGA/TCGA-SKCM.htseq_counts.tsv.gz',header=T,row.names = 1)
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
exp_tumor <- exp[,which(group=='Tumor')] %>% t %>% as.data.frame()
exp_tumor$sample <- rownames(exp_tumor)
clin <- read.table('rawdata/Melanoma/TCGA/TCGA-SKCM.survival.tsv',header=T,row.names = 1)
clin$sample <- substr(rownames(clin),1,15)
survdata <- merge(exp_tumor,clin,by.x='sample')
CRRS <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','ITGB1','C4BPA','JMJD7-PLA2G4B','CR2')
survdata1 <- survdata[,c('sample',CRRS,'OS','OS.time')]
survdata1$Risk=0.083*survdata1$F2-0.256*survdata1$P2RX+0.127*survdata1$PLAUR-0.163*survdata1$PRKCZ+0.098*survdata1$F12+0.083*survdata1$SERPINE1+
  0.076*survdata1$ITGB1-0.029*survdata1$C4BPA-0.021*survdata1$CR2-0.055*survdata1$GUCY1A1
survdata1$group <- ifelse(survdata1$Risk>median(survdata1$Risk),'High risk','Low risk')
write.xlsx(survdata1[,c(1,12:14,19)],file = 'files/Fig5.xlsx',sheetName = 'Fig5D_TCGA_melanoma',row.names = F,append = T)
png("figures/05D_OSKM_risk_TCGA.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(OS.time, OS) ~ group, data=survdata1),size=3,
                surv.median.line = "hv",
                title = "TCGA Melanoma",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()
