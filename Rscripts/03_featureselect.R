library(readxl)
CRGs <- read_xlsx('CRG/result/files/CRGs_209.xlsx')
MNC <- read.csv('CRG/result/files/MNC20.csv',sep = ',')
exp <- read.table('LUADdata/TCGA/TCGA-LUAD.htseq_counts.tsv.gz',header=T,row.names = 1)
colnames(exp) <- gsub('\\.','-',colnames(exp))
table(as.numeric(substr(colnames(exp),14,15))<10) #59 normal 526 tumor
load('LUADdata/TCGA/clinic.Rdata')
table(substr(colnames(exp),1,15)%in%CLINIC$sampleID) #临床数据都对应
colnames(exp) <- substr(colnames(exp),1,15)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
rownames(exp) <- str_sub(rownames(exp), start = 1, end = 15)
exp$ENSEMBL <- rownames(exp)
# bitr in clusterProfiler
df <- bitr( rownames(exp), fromType = "ENSEMBL", toType = c( "ENTREZID" ), OrgDb = org.Hs.eg.db )
#### 第三步，将GI号合并到DEG数据框内
exp <- merge(exp, df, by='ENSEMBL')
gf=bitr(rownames(exp), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
exp <- merge(exp,gf, by='ENTREZID')
exp<-exp[!duplicated(exp$SYMBOL),]
rownames(exp)<-exp$SYMBOL[!duplicated(exp$SYMBOL)]
which(colnames(exp)%in%c( "ENSEMBL", "ENTREZID","SYMBOL"))
exp <- exp[,-c(1,2,588)]
group <- ifelse(as.numeric(substr(colnames(exp),14,15))<10,'Tumor','Normal')
exp_mnc <- exp[rownames(exp)%in%MNC$Name,]
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
#########################F2基因正常和肿瘤boxplot
F2_dat <- t(exp_mnc)
F2_dat <- cbind(group,F2_dat )
F2_dat <- as.data.frame(F2_dat)
F2_dat$F2 <- as.numeric(F2_dat$F2)
cb <- ggboxplot(F2_dat,x = 'group',y = 'F2',fill='group')+theme_classic()+ylab('Expression')+xlab('')+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(textsize=5,comparisons = list(c("Tumor", "Normal")),test = wilcox.test,map_signif_level = F)+
  theme(axis.text.x = element_text(size=16,face='bold',color='black'),
        axis.text.y = element_text(size=16,face='bold',color='black'),
        text=element_text(size=16,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1))
ggsave('CRG/result/figs/08_box_F2exp_TCGA.tiff',cb,width = 5,height = 5,dpi=1000)
###########甲基化cg13527922  F2 boxplot
meth <- read.table('LUADdata/TCGA/TCGA.LUAD.sampleMap_HumanMethylation450.txt',header = T,row.names = 1)
colnames(meth) <- gsub('\\.','-',colnames(meth) )
F2_meth <- meth[rownames(meth)=='cg13527922',]
F2_meth <- as.data.frame(t(F2_meth))
F2_meth$group <- ifelse(as.numeric(substr(rownames(F2_meth),14,15))<10,'Tumor','Normal')
mb <- ggboxplot(F2_meth,x = 'group',y = 'cg13527922',fill='group')+theme_classic()+ylab('Methylation')+xlab('')+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(textsize=5,comparisons = list(c("Tumor", "Normal")),test = wilcox.test,map_signif_level = F)+
  theme(axis.text.x = element_text(size=16,face='bold',color='black'),
        axis.text.y = element_text(size=16,face='bold',color='black'),
        text=element_text(size=16,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1))
ggsave('CRG/result/figs/08_box_F2meth_TCGA.tiff',mb,width = 5,height = 5,dpi=1000)

source("~/project/R/my_luad/scripts/DEG_analysis.R")
filename<-paste0("CRG/extdata/DEG_TN_TCGA.Rdata")
runDEG("Tumor","Normal",exp,group,filename)
load(filename)
nrDEG<-DEG_limma_voom
nrDEG<-subset(nrDEG,nrDEG$P.Value<0.05)   #gene features减少到7276
logFC_cutoff <- with(nrDEG, mean( abs(logFC) ) + 2 * sd(abs(logFC)))
logFC_cutoff <- round(logFC_cutoff,2) #确定两位小数的logFC_cutoff=0.88
nrDEG$change = as.factor( ifelse( nrDEG$P.Value< 0.05 & abs(nrDEG$logFC) > logFC_cutoff,
                                  ifelse( nrDEG$logFC > logFC_cutoff , 'UP', 'DOWN' ), 'STABLE' ) )
DEGgenes=rownames(subset(nrDEG,change!='STABLE') ) #差异很大的gene features只剩下400个
intersect(DEGgenes,MNC$Name)  #只有F2是差异表达的基因
colnames(nrDEG)[1] <- 'log2FoldChange'
colnames(nrDEG)[4] <- 'pvalue'
source("~/project/R/my_luad/scripts/volcano.R")
volname <- "CRG/result/figs/07_TNtcga_volcano.png"
volcanoplot(nrDEG,volname)

clin_rna <- CLINIC[CLINIC$sampleID%in%colnames(exp),c("sampleID", "days_to_death" ,"days_to_last_followup","vital_status")] #574/585,11 个人没有临床数据
table(clin_rna$vital_status)
clin_rna$time <- ifelse(clin_rna$vital_status=='LIVING',clin_rna$days_to_last_followup/30,clin_rna$days_to_death/30)
clin_rna$status <- ifelse(clin_rna$vital_status=='LIVING',0,1)
table(as.numeric(substr(clin_rna$sampleID,14,15))<10) # 59 normal 515 tumor
clin_rna$group <-  ifelse(as.numeric(substr(clin_rna$sampleID,14,15))<10,'Tumor','Normal')
library(survival)
library(survminer)
ggsurvplot(survfit(Surv(time, status) ~ group, data=clin_rna),surv.median.line = "hv",palette = c('red','blue','grey'),pval = T,pval.size=10)

exp_tumor <- exp_mnc[,colnames(exp_mnc)%in%clin_rna$sampleID[clin_rna$group=='Tumor']] #515
library(dplyr)
exp_tumor <- exp_tumor %>% t%>%as.data.frame()
exp_tumor$sampleID <- rownames(exp_tumor)
tumor_clin <- as.data.frame(clin_rna[clin_rna$group=='Tumor'])
tumor_surv <- merge(tumor_clin,exp_tumor,by.y='sampleID')

tumor_surv$Group <- ifelse(as.numeric(tumor_surv$F2)>=median(as.numeric(tumor_surv$F2)),"High_F2","Low_F2")
ggsurvplot(survfit(Surv(time, status) ~ Group, data=tumor_surv),surv.median.line = "hv",palette = c('red','blue'),pval = T,pval.size=10)


#######################################lasso######################

exp_CRGs <- exp[rownames(exp)%in%CRGs$gene,]
exp_CRGs_tumor <- exp_CRGs[,colnames(exp_CRGs)%in%clin_rna$sampleID[clin_rna$group=='Tumor']] #515
exp_CRGs_tumor <- exp_CRGs_tumor %>% t%>%as.data.frame()
exp_CRGs_tumor$sampleID <- rownames(exp_CRGs_tumor)
CRGs_clin <- as.data.frame(clin_rna[clin_rna$group=='Tumor',])
CRGs_surv <- merge(CRGs_clin,exp_CRGs_tumor,by.y='sampleID')
CRGs_surv_dat <- CRGs_surv[,c(1,5,6,8:ncol(CRGs_surv))]
CRGs_surv_dat <- na.omit(CRGs_surv_dat )
CRGs_surv_dat$time <- round(CRGs_surv_dat$time,3)
CRGs_surv_dat <- CRGs_surv_dat[CRGs_surv_dat$time!=0,]
library(glmnet)
set.seed(2)
X.train <- as.matrix(CRGs_surv_dat[,8:ncol(CRGs_surv_dat)])
Y.train <- as.matrix(cbind(time=CRGs_surv_dat$time,status=CRGs_surv_dat$status))
fit.cv <- cv.glmnet(X.train,Y.train, family = "cox", nfolds = 10,alpha=1,keep=TRUE)
plot(fit.cv)
idmin = match(fit.cv$lambda.min, fit.cv$lambda)
plot(roc.glmnet(fit.cv$fit.preval, newy = Y.train)[[idmin]])
lsocoef<-as.numeric(coef(fit.cv ,s="lambda.min"))
lasso_genes <- colnames(X.train)[lsocoef!=0]
model <- glmnet(X.train,Y.train, family = "cox", alpha=1)
plot(model)
library(ggpubr)
combinde_coef <- c(0.08691,-0.31064,-0.11088,0.10607,-0.13589,0.09316,0.05993,0.13770,-0.02263,-0.01604)
gene <- c('F2',    'P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','CR2')
dmf <-cbind(gene,combinde_coef )
colnames(dmf) <- c('gene','coef')
dmf <- as.data.frame(dmf)
dmf$coef <- as.numeric(dmf$coef)
dmf <- dmf[order(dmf$coef,decreasing = T),]
dmf$gene <- factor(dmf$gene,levels=dmf$gene) 
library(ggalt)
dp <- ggplot(dmf,aes(y = gene, x=coef,xend=0))+
  geom_segment(aes(x=coef,xend=0,y=gene,yend=gene),color='grey',size=2)+
  geom_dumbbell(size_x=6, size_xend = 0,color="grey",colour_x = "black")+
  labs(x='Coefficient', y=NULL) +geom_vline(aes(xintercept = 0),lty=2,lwd=1)+
  theme_bw() + theme(axis.text.y = element_text(size=16,face='bold',color='black'),
                          axis.text.x = element_text(size=15,face='bold',color='black'),
                          text=element_text(size=16,face='bold'),axis.line = element_line(size=1),panel.grid = element_blank(),panel.border = element_rect(size=1.5))

ggsave('CRG/result/figs/09_dot_coef_multicox.tiff',dp,width = 6,height = 8,dpi=1000)


write_xlsx(dmf,'CRG/result/files/10_lassocoef.xlsx') 

#对lasso筛选的参考基因和重叠基因进行最优子集回归
library(My.stepwise)
bsrdata <- CRGs_surv_dat[,c('time','status',lasso_genes[1:12])]
rownames(bsrdata) <- CRGs_surv_dat$sampleID
genelist <- c(lasso_genes[1:12])
bsr <- My.stepwise.coxph(Time = 'time',Status = 'status',variable.list = genelist,data = CRGs_surv_dat,sle = 0.25,sls = 0.25)
#只选了F2

library(StepReg)
formula = Surv(time, status) ~ . - status

bsr1 <- stepwiseCox(formula,
                    bsrdata,
                    include=NULL,
                    selection=c("bidirection"),
                    select="HQ",
                    method=c("efron"),
                    sle=0.15,
                    sls=0.15,
                    weights=NULL,
                    best=NULL)


#                 coef exp(coef)   se(coef)          z     Pr(>|z|)
# F2        0.08690764 1.0907959 0.02798698  3.1052885 1.900935e-03
# P2RX1    -0.31063514 0.7329813 0.07095748 -4.3777648 1.199026e-05
# GUCY1A1  -0.11087861 0.8950474 0.08025990 -1.3814945 1.671270e-01
# PLAUR     0.10607296 1.1119030 0.07418423  1.4298586 1.527576e-01
# PRKCZ    -0.13589472 0.8729345 0.07042829 -1.9295473 5.366295e-02
# F12       0.09315952 1.0976368 0.05444703  1.7110119 8.707892e-02
# SERPINE1  0.05993044 1.0617627 0.06367458  0.9411988 3.466030e-01
# COL1A2    0.13770395 1.1476357 0.07126920  1.9321664 5.333898e-02
# C4BPA    -0.02262912 0.9776250 0.02937875 -0.7702546 4.411489e-01
# CR2      -0.01603878 0.9840892 0.03594258 -0.4462335 6.554286e-01

mulcox<- coxph(Surv(time, status) ~ F2+P2RX1+GUCY1A1+PLAUR+PRKCZ+F12+SERPINE1+COL1A2+C4BPA+CR2, data = bsrdata);mulcox

#               coef exp(coef) se(coef)      z Pr(>|z|)    
# F2        0.08691   1.09080  0.02799  3.105   0.0019 **
# P2RX1    -0.31064   0.73298  0.07096 -4.378  1.2e-05 ***
# GUCY1A1  -0.11088   0.89505  0.08026 -1.381   0.1671
# PLAUR     0.10607   1.11190  0.07418  1.430   0.1528
# PRKCZ    -0.13589   0.87293  0.07043 -1.930   0.0537 .
# F12       0.09316   1.09764  0.05445  1.711   0.0871 .
# SERPINE1  0.05993   1.06176  0.06367  0.941   0.3466
# COL1A2    0.13770   1.14764  0.07127  1.932   0.0533 .
# C4BPA    -0.02263   0.97762  0.02938 -0.770   0.4411
# CR2      -0.01604   0.98409  0.03594 -0.446   0.6554
forest <- ggforest(mulcox,  #coxph得到的Cox回归结果
                   data = bsrdata,  #数据集
                   main = 'Hazard ratio of all data',  #标题
                   cpositions = c(0.05, 0.15, 0.35),  #前三列数值的距离
                   fontsize = 1.4, #字体大小
                   refLabel = 'Reference', #相对变量的标签
                   noDigits = 3 #HR、95%CI的小数位
)

ggsave(filename="CRG/result/figs/09_muti_forest.png",forest, width=17, height=5)
covariates <- c('F2',    'P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','CR2')
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data =  bsrdata)})
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
#转换成数据框，并转置
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <-as.data.frame(res,stringsAsFactors=F)
HR=gsub("[\\(\\)]","",res$`HR (95% CI for HR)`)
HR=gsub("-"," ",HR)
HR=as.data.frame(do.call(cbind,strsplit(HR," ")),stringsAsFactors=F)
names(HR)=rownames(res)
HR <- t(HR)



unicox<- summary(coxph(Surv(time, status) ~ F2, data = bsrdata));unicox
# coef exp(coef) se(coef)     z Pr(>|z|)    
# F2 0.10843   1.11453  0.02621 4.137 3.51e-05 ***
unicox<- summary(coxph(Surv(time, status) ~ P2RX1, data = bsrdata));unicox
# coef exp(coef) se(coef)      z Pr(>|z|)    
# P2RX1 -0.21305   0.80811  0.05522 -3.858 0.000114 ***
unicox<- summary(coxph(Surv(time, status) ~ GUCY1A1, data = bsrdata));unicox
# coef exp(coef) se(coef)      z Pr(>|z|)   
# GUCY1A1 -0.19247   0.82491  0.06342 -3.035   0.0024 **
unicox<- summary(coxph(Surv(time, status) ~ PLAUR, data = bsrdata));unicox
# coef exp(coef) se(coef)     z Pr(>|z|)  
# PLAUR 0.1453    1.1564   0.0588 2.471   0.0135 *
unicox<- summary(coxph(Surv(time, status) ~ PRKCZ, data = bsrdata));unicox
# coef exp(coef) se(coef)      z Pr(>|z|)   
# PRKCZ -0.19426   0.82344  0.06782 -2.864  0.00418 **
unicox<- summary(coxph(Surv(time, status) ~ SERPINE1, data = bsrdata));unicox
# coef exp(coef) se(coef)     z Pr(>|z|)  
# SERPINE1 0.10528   1.11102  0.04841 2.175   0.0296 *
unicox<- summary(coxph(Surv(time, status) ~ COL1A2, data = bsrdata));unicox
# coef exp(coef) se(coef)     z Pr(>|z|)  
# COL1A2 0.08340   1.08698  0.04827 1.728    0.084 .
unicox<- summary(coxph(Surv(time, status) ~ C4BPA, data = bsrdata));unicox
# coef exp(coef) se(coef)      z Pr(>|z|)    
# C4BPA -0.08427   0.91918  0.02502 -3.368 0.000758 ***
unicox<- summary(coxph(Surv(time, status) ~ CR2 , data = bsrdata));unicox
# exp(coef) exp(-coef) lower .95 upper .95
# CR2    0.8969      1.115    0.8419    0.9556
unicox<- summary(coxph(Surv(time, status) ~ F12 , data = bsrdata));unicox
# coef exp(coef) se(coef)     z Pr(>|z|)  
# F12 0.13263   1.14182  0.05211 2.545   0.0109 *
save(exp,exp_CRGs,exp_CRGs_tumor,exp_mnc,clin_rna,lasso_genes,mulcox,fit.cv,model,bsr1,bsrdata,CRGs_surv_dat,DEGgenes,file = "CRG/extdata/featureselect_data.Rdata")
load( "CRG/extdata/featureselect_data.Rdata")

clin <- bsrdata
clin$riskscore <- 0.08691*as.numeric(clin$F2)+
  0.09316*as.numeric(clin$F12)+
  0.13770*as.numeric(clin$COL1A2)-
  0.01604*as.numeric(clin$CR2)+
  0.10607*as.numeric(clin$PLAUR)-
  0.31064*as.numeric(clin$P2RX1) -
  0.13589*as.numeric(clin$PRKCZ) -
  0.02263*as.numeric(clin$C4BPA ) +
  0.05993*as.numeric(clin$SERPINE1)-
  0.11088*as.numeric(clin$GUCY1A1)

clin$Group <- ifelse(clin$riskscore>median(clin$riskscore),'High risk','Low risk')
library(survival)
library(survminer)
ggsurvplot(survfit(Surv(time, status) ~ Group, data=clin),surv.median.line = "hv",palette = c('red','blue','grey'),pval = T,pval.size=10)

clin$riskscore1 <- 0.10843*as.numeric(clin$F2)
clin$riskscore2 <- -0.21305*as.numeric(clin$P2RX1)
clin$riskscore3 <- -0.19247*as.numeric(clin$GUCY1A1)
clin$riskscore4 <- 0.1453*as.numeric(clin$PLAUR)
clin$riskscore5 <- -0.19426*as.numeric(clin$PRKCZ)
clin$riskscore6 <- 0.10528*as.numeric(clin$SERPINE1)
clin$riskscore7 <- 0.08340*as.numeric(clin$COL1A2)
clin$riskscore8 <- -0.08427*as.numeric(clin$C4BPA)
clin$riskscore9 <- 0.8969*as.numeric(clin$CR2)
clin$riskscore10 <- 0.13263 *as.numeric(clin$F12)


pd <- clin[order(clin$riskscore,decreasing = F),]
pd$Sample <- 1:nrow(pd)
pd$Status <- ifelse(pd$status=='0','Alive','Dead')
library(survivalROC)
troc= survivalROC(Stime=pd$time,  
                  status=pd$status,   
                  marker = pd$riskscore,     
                  predict.time = 12, method="KM")
troc.1= survivalROC(Stime=pd$time,  
                    status=pd$status,   
                    marker = pd$riskscore1,     
                    predict.time = 12, method="KM")
troc.2= survivalROC(Stime=pd$time,  
                    status=pd$status,   
                    marker = pd$riskscore2,     
                    predict.time = 12, method="KM")
troc.3= survivalROC(Stime=pd$time,  
                    status=pd$status,   
                    marker = pd$riskscore3,     
                    predict.time = 12, method="KM")
troc.4= survivalROC(Stime=pd$time,  
                    status=pd$status,   
                    marker = pd$riskscore4,     
                    predict.time = 12, method="KM")
troc.5= survivalROC(Stime=pd$time,  
                    status=pd$status,   
                    marker = pd$riskscore5,     
                    predict.time = 12, method="KM")
troc.6= survivalROC(Stime=pd$time,  
                    status=pd$status,   
                    marker = pd$riskscore6,     
                    predict.time = 12, method="KM")
troc.7= survivalROC(Stime=pd$time,  
                    status=pd$status,   
                    marker = pd$riskscore7,     
                    predict.time = 12, method="KM")
troc.8= survivalROC(Stime=pd$time,  
                    status=pd$status,   
                    marker = pd$riskscore8,     
                    predict.time = 12, method="KM")
troc.9= survivalROC(Stime=pd$time,  
                    status=pd$status,   
                    marker = pd$riskscore9,     
                    predict.time = 12, method="KM")
troc.10= survivalROC(Stime=pd$time,  
                    status=pd$status,   
                    marker = pd$riskscore10,     
                    predict.time = 12, method="KM")
png("CRG/result/figs/09_rocall.tiff",width = 1500,height = 1500,res=200)
plot(troc$FP, troc$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=3,cex.axis=1.5,ann = F)
title(xlab="1-Specificity", ylab="Sensitivity",cex.lab=1.5,font.lab=4)
abline(0,1,col='gray',lty=2,lwd=2)
lines(troc.1$FP, troc.1$TP, type="l",lwd=2,col="blue", xlim=c(0,1), ylim=c(0,1))
lines(troc.2$FP, troc.2$TP, type="l",lwd=2,col="green", xlim=c(0,1), ylim=c(0,1))
lines(troc.3$FP, troc.3$TP, type="l",lwd=2,col="yellow", xlim=c(0,1), ylim=c(0,1))
lines(troc.4$FP, troc.4$TP, type="l",lwd=2,col="greenyellow", xlim=c(0,1), ylim=c(0,1))
lines(troc.5$FP, troc.5$TP, type="l",lwd=2,col="pink", xlim=c(0,1), ylim=c(0,1))
lines(troc.6$FP, troc.6$TP, type="l",lwd=2,col="black", xlim=c(0,1), ylim=c(0,1))
lines(troc.7$FP, troc.7$TP, type="l",lwd=2,col="navy", xlim=c(0,1), ylim=c(0,1))
lines(troc.8$FP, troc.8$TP, type="l",lwd=2,col="maroon", xlim=c(0,1), ylim=c(0,1))
lines(troc.9$FP, troc.9$TP, type="l",lwd=2,col="orange", xlim=c(0,1), ylim=c(0,1))
lines(troc.10$FP, troc.10$TP,lwd=2, type="l",col="purple", xlim=c(0,1), ylim=c(0,1))
legend(0.45,0.3,c(paste('Risk=',round(troc$AUC,3)),
                  paste('F2=',round(troc.1$AUC,3)),
                  paste('P2RX1=',round(troc.2$AUC,3)),
                  paste('GUCY1A1=',round(troc.3$AUC,3)),
                  paste('PLAUR=',round(troc.4$AUC,3)),
                  paste('PRKCZ=',round(troc.5$AUC,3)),
                  paste('SERPINE1=',round(troc.6$AUC,3)),
                  paste('COL1A2=',round(troc.7$AUC,3)),
                  paste('C4BPA=',round(troc.8$AUC,3)),
                  paste('CR2=',round(troc.9$AUC,3)),
                  paste('F12=',round(troc.10$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty=1,lwd=2,col = c('red','blue',"green","yellow","greenyellow","pink","black","navy","maroon","orange","purple"),
       bty='n',seg.len = 2,cex = 0.78)
dev.off()

troc= survivalROC(Stime=pd$time,  
                  status=pd$status,     
                  marker = pd$riskscore,     
                  predict.time = 12, method="KM")
troc_3= survivalROC(Stime=pd$time,  
                    status=pd$status,   
                    marker = pd$riskscore,     
                    predict.time = 36, method="KM")
troc_5= survivalROC(Stime=pd$time,  
                    status=pd$status,  
                    marker = pd$riskscore,     
                    predict.time = 60, method="KM")
png("CRG/result/figs/09_roc_year.tiff",width = 1500,height = 1500,res=200)
plot(troc$FP, troc$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=3,cex.axis=1.5,ann = F)
title(xlab="1-Specificity", ylab="Sensitivity",cex.lab=1.5,font.lab=4)
abline(0,1,col='gray',lty=2,lwd=2)
lines(troc_3$FP, troc_3$TP, type="l",col="blue", xlim=c(0,1), ylim=c(0,1),lwd=2)
lines(troc_5$FP, troc_5$TP, type="l",col="green", xlim=c(0,1), ylim=c(0,1),lwd=2)
legend(0.45,0.3,c(paste('AUC of 1-year =',round(troc$AUC,3)),paste('AUC of 3-year =',round(troc_3$AUC,3)),paste('AUC of 5-year =',round(troc_5$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd = 2,col = c('red','blue',"green"),bty='n',seg.len = 2,cex = 1.3)
dev.off()
library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
#风险评分box图
cb <- ggboxplot(pd,x = 'Status',y = 'riskscore',fill='Status')+theme_classic()+ylab('Riskscore')+xlab('')+scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_signif(textsize=5,comparisons = list(c("Dead", "Alive")),test = wilcox.test,map_signif_level = F)+
  theme(axis.text.x = element_text(size=16,face='bold',color='black'),
        axis.text.y = element_text(size=16,face='bold',color='black'),
        text=element_text(size=16,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1))
ggsave('CRG/result/figs/09_riskbox.tiff',cb,width = 5,height = 5,dpi=1000)

#点图
pd <- clin[order(clin$riskscore,decreasing = F),]
pd$Sample <- 1:nrow(pd)
sp <- ggplot(pd,aes(Sample,time,color=Status))+geom_point()+geom_vline(aes(xintercept = 251),lty=2,lwd=1)+
  theme_bw()+scale_colour_manual(values = c('blue','red'))+ylab('Survival month')+xlab('')+
  theme(axis.text.x = element_text(size=16,face='bold',color='black'),
        axis.text.y = element_text(size=16,face='bold',color='black'),
        text=element_text(size=16,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
rp <- ggplot(pd,aes(Sample,riskscore,color=Group))+geom_point()+geom_vline(aes(xintercept = 251),lty=2,lwd=1)+
  theme_bw()+scale_colour_manual(values = c('red','blue'))+ylab('Risk score')+xlab('')+
  theme(axis.text.x = element_text(size=16,face='bold',color='black'),
        axis.text.y = element_text(size=16,face='bold',color='black'),
        text=element_text(size=16,face='bold'),
        legend.title = element_blank(),legend.justification=c(0,1),legend.position=c(0.01,0.99),
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('CRG/result/figs/09_survpoint.tiff',sp,width = 5,height =3,dpi=1000)
ggsave('CRG/result/figs/09_riskpoint.tiff',rp,width = 5,height = 3,dpi=1000)

library(pheatmap)
gene <- c('F2', 'P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','CR2')
dat <- pd[,gene]
library(pheatmap)
group=pd$Group%>%as.data.frame()
colnames(group) <- 'group'
rownames(group) <- rownames(pd)
group$group <- factor(group$group,levels=c('Low risk','High risk'))
data <- t(dat)
p <- pheatmap(data ,show_colnames = F,cluster_cols =F,annotation_col = group,cellheight = 15,
              color = colorRampPalette(c("navy", "white", "firebrick3"))(50)) 
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(p, "CRG/result/figs/09_pheat.pdf")
save(pd,gene,file = 'CRG/extdata/tcga_plot.Rdata')
load('CRG/extdata/tcga_plot.Rdata')


