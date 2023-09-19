#1.CRGs+DEG+unicox+multicox-->TCGA multi-omics training+ validate 分三等分，验证集中AUC 0.65
library(ChAMP)
methy.beta <- read.table('rawdata/TCGA_LUAD/TCGA.LUAD.sampleMap_HumanMethylation450.txt',header = T,row.names = 1)
xx <- c()
for(i in 1:ncol(methy.beta)) {
  if(table(is.na(methy.beta[,i]))[2]/nrow(methy.beta)<0.2){
    xx <- c(xx,i)
  }
}
methy.beta <- methy.beta[,xx]

table(substr(colnames(methy.beta),14,15))
group <- ifelse(substr(colnames(methy.beta),14,15)>10,'Normal','Tumor')
pd <- cbind(colnames(methy.beta),group)
colnames(pd) <- c('sample','sampleType')
pd <- as.data.frame(pd)

myDMP <- champ.DMP(beta=methy.beta,
                   pheno=pd$sampleType,
                   compare.group = c('Tumor','Normal'),
                   adjPVal = 0.05,
                   adjust.method = 'BH',
                   arraytype = '450K')

dim(myDMP[[1]][myDMP[[1]]$deltaBeta>0.15&myDMP[[1]]$adj.P.Val<0.05,])
DEDMP <- myDMP[[1]][myDMP[[1]]$deltaBeta>0.15&myDMP[[1]]$adj.P.Val<0.05,]

library(readxl)
CRGs <- read_xlsx('files/01_CRGs_209.xlsx')
CRDMP <- DEDMP[DEDMP$gene%in%CRGs$gene,]
methylation <- methy.beta[rownames(CRDMP),substr(colnames(methy.beta),14,15)<10]
colnames(methylation) <- gsub('\\.','-',colnames(methylation))


library()
MUT <- read.table('rawdata/TCGA_LUAD/mc3_gene_level_LUAD_mc3_gene_level.txt',header = T)
mutation <- MUT[MUT$sample%in%CRGs$gene,]
rownames(mutation) <- mutation$sample
mutation <- mutation[,-1]
mutation <- mutation[rowSums(mutation)!=0,substr(colnames(mutation),14,15)<10]
colnames(mutation) <- gsub('\\.','-',colnames(mutation))

library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
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
load('Rdata/DEG_TN_TCGA.Rdata')
DEGs <- rownames(DEG_limma_voom)[DEG_limma_voom$adj.P.Val<0.05&abs(DEG_limma_voom$logFC)>log2(1)]
CRGs <- CRGs$gene[CRGs$gene%in%DEGs]
exp_CRGs <- exp[rownames(exp)%in%CRGs,substr(colnames(exp),14,15)<10]

samples <- intersect(intersect(colnames(methylation),colnames(mutation)),colnames(exp_CRGs))

methylation <- methylation[,colnames(methylation)%in%samples]
mutation <- mutation[,colnames(mutation)%in%samples]
exp_CRGs <- exp_CRGs[,colnames(exp_CRGs)%in%samples]
rownames(mutation) <- paste0(rownames(mutation),'_mut')
rownames(exp_CRGs) <- paste0(rownames(exp_CRGs),'_exp')
library(dplyr)

methylation <- t(methylation)%>% as.data.frame()
methylation$sample <- rownames(methylation)
mutation <- t(mutation)%>% as.data.frame()
mutation$sample <- rownames(mutation)
exp_CRGs <- t(exp_CRGs)%>% as.data.frame()
exp_CRGs$sample <- rownames(exp_CRGs)
data <- merge(methylation,mutation,by='sample')
df <-merge(data,exp_CRGs,by='sample') 

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
clin <- clin_tcga[,c('sampleID','OS','OS.E')]
colnames(clin)[1] <- 'sample'
clin <- na.omit(clin)
clin <- clin[clin$OS!=0,]
CRGs_surv_dat <- merge(df,clin,by.x='sample')
CRGs_surv_dat[is.na(CRGs_surv_dat)] <- 0
saveRDS(CRGs_surv_dat,'Rdata/CRGs_surv_dat.RDS')
fold <- sample(c(rep(1,146),rep(2,146),rep(3,146)),438,replace = F)
library(glmnet)
set.seed(10)
X.train <- as.matrix(CRGs_surv_dat[fold!=3,2:(ncol(CRGs_surv_dat)-2)])
Y.train <- as.matrix(cbind(time=CRGs_surv_dat$OS[fold!=3],status=CRGs_surv_dat$OS.E[fold!=3]))
fit.cv <- cv.glmnet(X.train,Y.train, family = "cox", nfolds = 5,alpha=1,keep=TRUE)
fit.cv$lambda.min
#0.07656243
png("figures/02B_lassocoef.png",width = 1500,height = 1500,res=200)
plot(fit.cv, lwd=3,cex.axis=1.5,ann = F)
title(xlab="log(λ)", ylab="Partial likelihood Deviance",cex.lab=1.5,font.lab=2)
dev.off()

idmin = match(fit.cv$lambda.min, fit.cv$lambda)
lsocoef<-as.numeric(coef(fit.cv ,s="lambda.min"))
lasso_genes <- colnames(X.train)[lsocoef!=0]
model <- glmnet(X.train,Y.train, family = "cox", alpha=1)

# Fig2A
png("figures/02A_coefficient.png",width = 1500,height = 1500,res=200)
plot(model, lwd=3,cex.axis=1.5,ann = F)
title(xlab="log(λ)", ylab="Coefficients",cex.lab=1.5,font.lab=2)
dev.off()

library(survival)
library(survminer)
library(rms)
library(survivalROC)
library(dplyr)
lxdata1 <- CRGs_surv_dat[fold!=3,]
# covariates <- c(bsrgenes)
covariates <- c(lasso_genes)
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS, OS.E)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = lxdata1)})
#提取HR，95%置信区间和p值
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         coef <- round(x$coef[1],3)
                         #获取HR
                         HR <-signif(x$coef[2], digits=2);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR,coef)
                         names(res)<-c("p.value","HR (95% CI for HR)",'coef')
                         return(res)
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <-as.data.frame(res,stringsAsFactors=F)
HR=gsub("[\\(\\)]","",res$`HR (95% CI for HR)`)
HR=gsub("-","",HR)
HR=as.data.frame(do.call(cbind,strsplit(HR," ")),stringsAsFactors=F)
names(HR)=rownames(res)
HR <- t(HR)

signatures <- rownames(res)[res$p.value<0.05]

mulcox<- coxph(Surv(OS, OS.E) ~ cg19707677+cg24815934+cg04914221+PRKG2_mut+C7_mut+CR2_exp+F2_exp+F10_exp+GP9_exp+P2RX1_exp, data = lxdata1 );mulcox

forest <- ggforest(mulcox,  #coxph得到的Cox回归结果
                   data = lxdata1,  #数据集
                   main = 'Hazard ratio of all data',  #标题
                   cpositions = c(0.05, 0.15, 0.35),  #前三列数值的距离
                   fontsize = 2, #字体大小
                   refLabel = 'Reference', #相对变量的标签
                   noDigits = 3 #HR、95%CI的小数位
)
ggsave(filename="02D_muti_forest.png",forest, width=28, height=14)

lxdata2 <- CRGs_surv_dat[fold==3,]
lxdata2$Risk <- predict(mulcox,lxdata2)
troc.1= survivalROC(Stime=lxdata2$OS,  
                    status=lxdata2$OS.E,   
                    marker = lxdata2$Risk,     
                    predict.time = 12, method="KM")
troc.1$AUC
#0.65
