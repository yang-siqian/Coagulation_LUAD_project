library(dplyr)
library(writexl)
library(GEOquery)
library(stringr)
gene <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7-PLA2G4B','CR2')
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
save(exp117,clindata117,file = "Rdata/GSE13213.Rdata")

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
save(exp226,clindata226,file = "Rdata/GSE31210.Rdata")

#############################################################
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
  save(exp398,clindata398,file = "Rdata/GSE72094.Rdata")

##########
gset <- getGEO('GSE78220', GSEMatrix =TRUE, AnnotGPL=FALSE)
exp <-read_xlsx('GSE78220_PatientFPKM.xlsx') %>%as.data.frame()
colnames(exp) <- gsub('.baseline','',colnames(exp))
colnames(exp) <- gsub('.OnTx','',colnames(exp))
pd <- pData(gset[[1]])
rownames(exp) <- exp$Gene
exp <- exp[,-1]
exp <- t(exp[rownames(exp)=='COL11A1',])%>%as.data.frame()
exp$title <- rownames(exp)
pd <- pd[,c(1,8,63,70)]
pd$response <- ifelse(pd$source_name_ch1%in%c('tumor biopsy_Complete Response','tumor biopsy_Partial Response'),'Complete or Partial Response','Progressive Disease')
pd$status <- ifelse(pd$`vital status:ch1`=='Dead',1,0)
pd$time <- pd$`overall survival (days):ch1`
pd <- pd[,-c(2:4)]
x=merge(pd,exp,by.x='title')
library(xlsx)
write.xlsx(x,file='immunotherapy_yzj.xlsx',sheetName = 'GSE78220')

#############
gset <- getGEO('GSE63557', GSEMatrix =TRUE, AnnotGPL=FALSE)

if (length(gset) > 1) idx <- grep("GPL19103", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
pd <- pData(gset)
ex <- exprs(gset)
pd$group <- gsub('group: ','',pd$characteristics_ch1.2)
colnames(pd)[2] <- 'sample'
pd <- pd[pd$group!='Untreated',c('sample','group')]
ex <- ex[,pd$sample]%>%as.data.frame()
group_list <- pd$group
source("~/project/R/my_luad/scripts/DEG_analysis.R")
filename<-paste0("DEG_r-n.Rdata")
runDEG("Responder","Nonresponder",ex,group_list,filename)
load(filename)
nrDEG<-subset(nrDEG,nrDEG$pvalue<0.05)   #gene features减少到15963
logFC_cutoff <- with(nrDEG, mean( abs(log2FoldChange) ) + 2 * sd(abs(log2FoldChange)))
logFC_cutoff <- round(logFC_cutoff,2) #确定两位小数的logFC_cutoff
nrDEG$change = as.factor( ifelse( nrDEG$pvalue < 0.05 & abs(nrDEG$log2FoldChange) > logFC_cutoff,
                                  ifelse( nrDEG$log2FoldChange > logFC_cutoff , 'UP', 'DOWN' ), 'STABLE' ) )
#12814_at
nrDEG[rownames(nrDEG)=='12814_at',]
DEG_limma_voom[rownames(DEG_limma_voom)=='12814_at',]

ex1 <- t(ex['12814_at',])%>%as.data.frame()
colnames(ex1) <- 'COL11A1'
ex1$sample <- rownames(ex1)
x <-merge(ex1,pd,by.x='sample') 
write.xlsx(x,file='immunotherapy_yzj.xlsx',sheetName = 'GSE63557',append = T)
write.xlsx(ex,file='immunotherapy_yzj.xlsx',sheetName = 'GSE63557_exp',append = T)
write.xlsx(as.data.frame(pd),file='immunotherapy_yzj.xlsx',sheetName = 'GSE63557_group',append = T)


##############
gset <- getGEO('GSE161537', GSEMatrix =TRUE, AnnotGPL=FALSE)
pd <- pData(gset[[1]])
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE161537", "file=GSE161537_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
# pre-filter low count genes
# keep genes with at least 2 counts > 10
keep <- rowSums( tbl >= 10 ) >= 2
tbl <- tbl[keep, ]
library(clusterProfiler)
library(org.Hs.eg.db)

# log transform raw counts
# instead of raw counts can display vst(as.matrix(tbl)) i.e. variance stabilized counts
dat <- log10(tbl + 1)%>%as.data.frame()
df <-bitr(rownames(dat), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
dat$ENTREZID <- rownames(dat)
exp <- merge(dat,df,by.x='ENTREZID')
rownames(exp) <- exp$SYMBOL
exp <- exp[,!colnames(exp)%in%c('ENTREZID','SYMBOL')]
exp1 <- exp[rownames(exp)%in%CRRS,]
exp1 <- t(exp1)%>%as.data.frame()
exp1$sample <- rownames(exp1)
# exp <- t(dat[rownames(dat)=='1301',])%>%as.data.frame()
# colnames(exp) <- 'COL11A1'
pd <- pd[,c(2,16,17,50)]
colnames(pd) <- c('sample','time(month)','status','response')
pd$`time(month)` <- substr(pd$`time(month)`,13,30)
pd$status <- substr(pd$status,13,14)
exp$sample <- rownames(exp)
x <- merge(exp1,pd,by.x='sample')
x$time <- x$`time(month)`
mulcox<- coxph(Surv(time, status) ~ F2+GUCY1A1+ITGB1+SERPINE1+PLAUR+PRKCZ, data = x );mulcox


x$Risk <- 0.083*x$F2-0.055*exp1$GUCY1A1+0.076*x$ITGB1+0.083*x$SERPINE1+0.127*x$PLAUR-0.163*x$PRKCZ

x$group <- ifelse(x$response%in%c('PD','SD'),'SD/PD','CR/PR')
p <- ggboxplot(x,x="group",y="Risk",fill = "group")+theme_classic()+ylab('Risk')+xlab('')+ggtitle('')+ylim(0,1)+
  geom_signif(size=1,textsize=9,comparisons = list(c('SD/PD','CR/PR')),map_signif_level = T)+
  theme( plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
         axis.text.x = element_text(size=18,face='bold',color='black'),
         axis.text.y = element_text(size=18,face='bold',color='black'),
         text=element_text(size=25,face='bold'),
         legend.title = element_blank(),legend.position = "none",
         axis.line = element_line(size=1))

write.xlsx(x,file='immunotherapy_yzj.xlsx',sheetName = 'GSE161537',append = T)


