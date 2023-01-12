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

