###################################################TCGA Data Process###################################################
################z-score标准化和log转换###################

RNA_TCGA <- function(filename,group){

  exp <- read.table(filename,header=T,row.names = 1)
  colnames(exp) <- gsub('\\.','-',colnames(exp))
  colnames(exp) <- substr(colnames(exp),1,15)
  library(stringr)
  rownames(exp) <- str_sub(rownames(exp), start = 1, end = 15)
  exp$ENSEMBL <- rownames(exp)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  df <- bitr( rownames(exp), fromType = "ENSEMBL", toType = c( "ENTREZID" ), OrgDb = org.Hs.eg.db )
  exp <- merge(exp, df, by='ENSEMBL')
  gf=bitr(rownames(exp), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
  exp <- merge(exp,gf, by='ENTREZID')
  exp<-exp[!duplicated(exp$SYMBOL),]
  rownames(exp)<-exp$SYMBOL[!duplicated(exp$SYMBOL)]
  exp <- exp[,-which(colnames(exp)%in%c( "ENSEMBL", "ENTREZID","SYMBOL"))]
  library(dplyr)
  if(group=='tumor'){
    exp <- exp[,as.numeric(substr(colnames(exp),14,15))<10] 
  }else if(group=='normal'){
    exp <- exp[,as.numeric(substr(colnames(exp),14,15))>=10] 
  }else if(group=='all'){
    exp <- exp
  }
  exp_scale <- scale(exp) %>% as.data.frame()
  return(exp_scale)
  
}


exp <- RNA_TCGA('rawdata/TCGA_LUAD/TCGA-LUAD.htseq_counts.tsv.gz','tumor')

CLINIC_TCGA <- function(filename,group){
  load(filename)
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
  if(group=='tumor'){
    clin_tcga <- clin_tcga[clin_tcga$group=='Tumor',] 
  }else if(group=='normal'){
    clin_tcga <- clin_tcga[clin_tcga$group=='Normal',] 
  }else if(group=='all'){
    clin_tcga <- clin_tcga
  }
  clin_tcga <- clin_tcga[!is.na(clin_tcga$OS),]
  clin_tcga <- clin_tcga[!is.na(clin_tcga$OS.E),]
  clin_tcga <- clin_tcga[clin_tcga$OS!=0,]
  return(clin_tcga)
  
}

clin <- CLINIC_TCGA('Rdata/clinic.Rdata','tumor')

sampleID <- intersect(clin$sampleID,colnames(exp)) #506
library(readxl)
CRGs <- read_xlsx('files/01_CRGs_209.xlsx')
exp_CRGs <- exp[CRGs$gene,sampleID] %>% t %>% as.data.frame()
exp_CRGs$sampleID <- rownames(exp_CRGs)
CRGs_surv_dat <- clin[clin$sampleID%in%sampleID,c('sampleID','OS','OS.E')]
CRGs_surv_dat <- merge(CRGs_surv_dat,exp_CRGs,by.x='sampleID')
CRGs_surv_dat <- na.omit(t(CRGs_surv_dat))%>% t %>% as.data.frame()
rownames(CRGs_surv_dat) <- CRGs_surv_dat$sampleID

###################################################FUSCC Data Process###################################################



###################################################GEO Data Process###################################################
#######################GSE30219########################
library(dplyr)
library(writexl)
library(GEOquery)
library(stringr)
gene <- c('F2','P2RX1','PLAUR','COL1A2','C4BPA')
gset <- getGEO('GSE30219', GSEMatrix =TRUE, AnnotGPL=FALSE)
exp <- exprs(gset[[1]])%>%as.data.frame()
boxplot(exp[1:100,1:20])
exp[1:4,1:4]
pd <- pData(gset[[1]])
#判断pd的行名顺序与exp列名完全一致
p <- identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
# 提取芯片平台编号
gpl <- getGEO(gset[[1]]@annotation)
gpl <- gpl@dataTable@table
gpl <- gpl[,c(1,11)]
gene%in%x3$`Gene Symbol`
exp$ID <- rownames(exp)
x1 <- merge(gpl,exp,by.x='ID') 
x2 <- x1[x1$`Gene Symbol`!='',]
x3 <- x2[!duplicated(x2$`Gene Symbol`),]
rownames(x3) <- x3$`Gene Symbol` 
#23520
x4 <- x3[,3:ncol(x3)]
exp308 <- x4
#无RFS
clindata <- pd[,c('geo_accession','gender:ch1',"age at surgery:ch1" ,"pm stage:ch1","pn stage:ch1","pt stage:ch1" ,"disease free survival in months:ch1" ,"follow-up time (months):ch1" ,"relapse (event=1; no event=0):ch1" , "status:ch1" )]
colnames(clindata) <- c('sample','gender','age','pM','pN','pT','RFS','OS','RFS.E','OS.E')
clindata293 <- clindata[clindata$OS.E!='NTL',]
clindata293$Status <- ifelse(clindata293$OS.E=='DEAD',1,0)
exp293 <- exp308[,colnames(exp308)%in%clindata$sample]
save(exp293,clindata293,file = "Rdata/GSE30219.Rdata")
#######################GSE41271########################
gene <- c('F2','P2RX1','PLAUR','COL1A2','C4BPA')
gset <- getGEO('GSE41271', GSEMatrix =TRUE, AnnotGPL=FALSE)
exp <- exprs(gset[[1]])%>%as.data.frame()
boxplot(exp[1:100,1:20])
exp[1:4,1:4]
pd <- pData(gset[[1]])
#判断pd的行名顺序与exp列名完全一致
p <- identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
# 提取芯片平台编号
gpl <- getGEO(gset[[1]]@annotation)
gpl <- gpl@dataTable@table
gpl <- gpl[,c(1,13)]
gene%in%x3$Symbol
exp$ID <- rownames(exp)
x1 <- merge(gpl,exp,by.x='ID') 
x2 <- x1[x1$Symbol!='',]
x3 <- x2[!duplicated(x2$Symbol),]
rownames(x3) <- x3$Symbol 
#25440
x4 <- x3[,3:ncol(x3)]
exp276 <- x4
#无RFS
pd <- pd[pd$`histology:ch1`=='Adenocarcinoma',]
pd <- pd[-182,]
pd$`date of birth:ch1`[173] <- '1942-03-20'
pd$`date of birth:ch1`[175] <- '1946-10-10'
pd$`date of birth:ch1`[176] <- '1936-04-21'

age <- c()
for(i in 1:182){
  print(i)
  date1 <- as.Date(pd$`date of birth:ch1`[i])
  date2 <- as.Date(pd$`date of surgery:ch1`[i])
  a <- seq(from=date1,to=date2,by='day')
  age[i] <- (length(a)-1)/365
}

pd$age <- age

RFS <- c()
for(i in 1:182){
  print(i)
  date1 <-as.Date(pd$`date of surgery:ch1`[i])
  date2 <- as.Date(pd$`last follow-up recurrence:ch1`[i])
  a <- seq(from=date1,to=date2,by='day')
  RFS[i] <- (length(a)-1)/30
}
pd$RFS <- RFS


OS <- c()
for(i in 1:182){
  print(i)
  date1 <-as.Date(pd$`date of surgery:ch1`[i])
  date2 <- as.Date(pd$`last follow-up survival:ch1`[i])
  a <- seq(from=date1,to=date2,by='day')
  OS[i] <- (length(a)-1)/30
}
pd$OS <- OS
pd$RFS.E <- ifelse(pd$`recurrence:ch1`=='Y',1,0)
pd$OS.E <- ifelse(pd$`vital statistics:ch1`=='D',1,0)
clindata <- pd[,c('geo_accession','gender:ch1',"age" ,"final patient stage:ch1",'RFS','OS','RFS.E','OS.E' )]
colnames(clindata) <- c('sample','gender','age','stage','RFS','OS','RFS.E','OS.E')
clindata182 <- clindata
exp182 <- exp276[,colnames(exp276)%in%clindata$sample]
save(exp182,clindata182,file = "Rdata/GSE41271.Rdata")
#######################GSE68465########################
library(dplyr)
library(writexl)
library(GEOquery)
library(stringr)
gene <- c('F2','P2RX1','PLAUR','COL1A2','C4BPA')
gset <- getGEO('GSE68465', GSEMatrix =TRUE, AnnotGPL=FALSE)
exp <- exprs(gset[[1]])%>%as.data.frame()
exp <- log10(exp+1)
boxplot(exp[1:100,1:20])
exp[1:4,1:4]
pd <- pData(gset[[1]])
pd <- pd[pd$`disease_state:ch1`=="Lung Adenocarcinoma",]
pd <- pd[pd$`mths_to_last_clinical_assessment:ch1`!='--',]
#判断pd的行名顺序与exp列名完全一致
p <- identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
# 提取芯片平台编号
gpl <- getGEO(gset[[1]]@annotation)
gpl <- gpl@dataTable@table
gpl <- gpl[,c(1,11)]
gene%in%x3$`Gene Symbol`
exp$ID <- rownames(exp)
x1 <- merge(gpl,exp,by.x='ID') 
x2 <- x1[x1$`Gene Symbol`!='',]
x3 <- x2[!duplicated(x2$`Gene Symbol`),]
rownames(x3) <- x3$`Gene Symbol` 
#13515
x4 <- x3[,3:ncol(x3)]
exp462 <- x4
#无RFS
clindata <- pd[,c('geo_accession',"Sex:ch1","age:ch1" ,"disease_stage:ch1" ,'mths_to_last_clinical_assessment:ch1', "vital_status:ch1" )]
colnames(clindata) <- c('sample','gender','age','stage','OS','OS.E')
clindata433 <- clindata
clindata433$Status <- ifelse(clindata433$OS.E=='Dead',1,0)
exp433 <- exp462[,colnames(exp462)%in%clindata433$sample]
save(exp433,clindata433,file = "Rdata/GSE68465.Rdata")
#######################GSE37745########################
library(dplyr)
library(writexl)
library(GEOquery)
library(stringr)
gene <- c('F2','P2RX1','PLAUR','COL1A2','C4BPA')
gset <- getGEO('GSE37745', GSEMatrix =TRUE, AnnotGPL=FALSE)
exp <- exprs(gset[[1]])%>%as.data.frame()

boxplot(exp[1:100,1:20])
exp[1:4,1:4]
pd <- pData(gset[[1]])
pd <- pd[pd$`histology:ch1`=="adeno",]
pd <- pd[pd$`mths_to_last_clinical_assessment:ch1`!='--',]
#判断pd的行名顺序与exp列名完全一致
p <- identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
# 提取芯片平台编号
gpl <- getGEO(gset[[1]]@annotation)
gpl <- gpl@dataTable@table
gpl <- gpl[,c(1,11)]
gene%in%x3$`Gene Symbol`
exp$ID <- rownames(exp)
x1 <- merge(gpl,exp,by.x='ID') 
x2 <- x1[x1$`Gene Symbol`!='',]
x3 <- x2[!duplicated(x2$`Gene Symbol`),]
rownames(x3) <- x3$`Gene Symbol` 
#23520
x4 <- x3[,3:ncol(x3)]
exp196 <- x4

clindata <- pd[,c('geo_accession',"gender:ch1","age:ch1" ,"tumor stage:ch1" ,'days to determined death status:ch1', "dead:ch1","days to recurrence / to last visit:ch1" ,"recurrence:ch1"  )]
colnames(clindata) <- c('sample','gender','age','stage','OS','OS.E','RFS','RFS.E')
clindata106 <- clindata
clindata106$Time <- clindata106$OS/30
clindata106$Status <- ifelse(clindata106$OS.E=='yes',1,0)
exp106 <- exp196[,colnames(exp196)%in%clindata106$sample]
save(exp106,clindata106,file = "Rdata/GSE37745.Rdata")

#######################GSE3141########################
library(dplyr)
library(writexl)
library(GEOquery)
library(stringr)
gene <- c('F2','P2RX1','PLAUR','COL1A2','C4BPA')
gset <- getGEO('GSE3141', GSEMatrix =TRUE, AnnotGPL=FALSE)
exp <- exprs(gset[[1]])%>%as.data.frame()
exp <- log2(exp+1)
boxplot(exp[1:100,1:20])
exp[1:4,1:4]
pd <- pData(gset[[1]])

#判断pd的行名顺序与exp列名完全一致
p <- identical(rownames(pd),colnames(exp));p
if(!p) exp = exp[,match(rownames(pd),colnames(exp))]
# 提取芯片平台编号
gpl <- getGEO(gset[[1]]@annotation)
gpl <- gpl@dataTable@table
gpl <- gpl[,c(1,11)]
gene%in%x3$`Gene Symbol`
exp$ID <- rownames(exp)
x1 <- merge(gpl,exp,by.x='ID') 
x2 <- x1[x1$`Gene Symbol`!='',]
x3 <- x2[!duplicated(x2$`Gene Symbol`),]
rownames(x3) <- x3$`Gene Symbol` 
#23520
x4 <- x3[,3:ncol(x3)]
exp112 <- x4

clindata <- pd[,c('geo_accession',"characteristics_ch1"  )]
a <- as.data.frame(t(matrix(unlist(strsplit(pd$characteristics_ch1,';')),nrow = 3,ncol = 111)))
clindata <- cbind(clindata,a)
clindata <- clindata[ grep('*A',clindata$V1),]
clindata$OS <- substr(clindata$V2,15,30)
clindata$OS.E <-substr(clindata$V3,27,27)

clindata58 <- clindata[,c(1,6,7)]
colnames(clindata58)[1] <- 'sample'

exp58<- exp112[,colnames(exp112)%in%clindata58$sample]
save(exp58,clindata58,file = "Rdata/GSE3141.Rdata")

