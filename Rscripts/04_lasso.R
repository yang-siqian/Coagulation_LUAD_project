load( "CRG/extdata/featureselect_data.Rdata")
load('LUADdata/TCGA/clinic.Rdata')
clin_tcga <- CLINIC[,c("sampleID","gender","age_at_initial_pathologic_diagnosis","pathologic_stage","pathologic_T",
                       "pathologic_N","pathologic_M","vital_status","tobacco_smoking_history_indicator","new_neoplasm_event_type","radiation_therapy",
                       "additional_pharmaceutical_therapy",'history_of_neoadjuvant_treatment','targeted_molecular_therapy','days_to_additional_surgery_locoregional_procedure',
                       'days_to_additional_surgery_metastatic_procedure')]
colnames(clin_tcga) <- c("sampleID",'Gender','Age','Stage',"Tumor","Lymph_Node","Metastasis",'Status','Smoking','New_tumor_type',"radiation_therapy",
                         'pharmaceutical_therapy', 'neoadjuvant_treatment','targeted_molecular_therapy','additional_surgery_locoregional_days','additional_surgery_metastatic_days')
clin_tcga <- as.data.frame(clin_tcga)
rownames(clin_tcga) <- clin_tcga$sampleID
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
clin_tumor <- merge(clin_tumor,clin_rna,by.x='sampleID')
clin_tumor <- clin_tumor[clin_tumor$time!=0,]
# clin_tumor <- clin_tumor[,-c(10:12)]
clin_tumor$RFS <- clin_tumor$additional_surgery_locoregional_days
clin_tumor$RFS[which(is.na(clin_tumor$RFS))] <- clin_tumor$days_to_death[which(is.na(clin_tumor$RFS))]
clin_tumor$RFS[which(is.na(clin_tumor$RFS))] <- clin_tumor$days_to_last_followup[which(is.na(clin_tumor$RFS))]
clin_tumor$RFS <- round(clin_tumor$RFS/30,3)
clin_tumor$RFS.E <- clin_tumor$status
clin_tumor$RFS.E[grep('Recurrence',clin_tumor$New_tumor_type)] <- 1

library(glmnet)
set.seed(10)
X.train <- as.matrix(CRGs_surv_dat[,8:ncol(CRGs_surv_dat)])
Y.train <- as.matrix(cbind(time=CRGs_surv_dat$time,status=CRGs_surv_dat$status))
fit.cv <- cv.glmnet(X.train,Y.train, family = "cox", nfolds = 10,alpha=1,keep=TRUE)
fit.cv$lambda.min
#0.04818326
png("CRG/result/figs/03_lassocoef.tiff",width = 1500,height = 1500,res=200)
plot(fit.cv, lwd=3,cex.axis=1.5,ann = F)
title(xlab="log(λ)", ylab="Partial likelihood Deviance",cex.lab=1.5,font.lab=2)
dev.off()
idmin = match(fit.cv$lambda.min, fit.cv$lambda)
lsocoef<-as.numeric(coef(fit.cv ,s="lambda.min"))
lasso_genes <- colnames(X.train)[lsocoef!=0]
model <- glmnet(X.train,Y.train, family = "cox", alpha=1)
png("CRG/result/figs/03_coefficient.tiff",width = 1500,height = 1500,res=200)
plot(model, lwd=3,cex.axis=1.5,ann = F)
title(xlab="log(λ)", ylab="Coefficients",cex.lab=1.5,font.lab=2)
dev.off()

#对lasso筛选的参考基因和重叠基因进行最优子集回归
library(My.stepwise)
bsrdata <- CRGs_surv_dat[,c('time','status',lasso_genes)]
colnames(bsrdata)[15] <- gsub('-','_',colnames(bsrdata)[15])
lasso_genes[13] <- gsub('-','_',lasso_genes[13])
rownames(bsrdata) <- CRGs_surv_dat$sampleID
genelist <- c(lasso_genes)
bsr <- My.stepwise.coxph(Time = 'time',Status = 'status',variable.list = genelist,data = bsrdata,sle = 0.25,sls = 0.25)
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
####################################dot
library(ggpubr)
coef <- round(bsr1$Coefficients[,1],3)
gene <- bsr1$Variables
gene[10] <- gsub('_','-',gene[10] )
dmf <-cbind(gene,coef )
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

ggsave('CRG/result/figs/03_dot_coef_bsr.tiff',dp,width = 6,height = 6,dpi=1000)

####################################################
bsrgenes <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7_PLA2G4B','CR2')
# bsrgenes <- bsr1$Variables
bsrgenes[10] <- gsub('_','-',bsrgenes[10])
hubdat <-  exp_CRGs_tumor[,c('sampleID',bsrgenes)]
modeldat<- merge(hubdat,clin_tumor,by.x = 'sampleID')
lxdata <- modeldat
lxdata1 <- merge(exp_CRGs_tumor,clin_tumor,by.x = 'sampleID')
lxdata1$time <- round(lxdata1$time,3)
lxdata1$Age <- ifelse(lxdata$Age=='young',1,2 )
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$Stage <- ifelse(lxdata$Stage=='I',1,ifelse(lxdata$Stage=='II',2 ,ifelse(lxdata$Stage=='III',3,ifelse(lxdata$Stage=='IV',4,0))))
lxdata1$Stage <- factor(lxdata1$Stage,levels=c(0,1,2,3,4),labels=c('UNKOWN','I','II','III','IV'))
bsrgenes[10] <- gsub('-','_',bsrgenes[10])
colnames(lxdata1)[11] <- gsub('-','_',colnames(lxdata1)[11])
library(survival)
library(survminer)
library(rms)
library(survivalROC)
library(dplyr)
covariates <- c(bsrgenes)
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))
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
HR=gsub("-"," ",HR)
HR=as.data.frame(do.call(cbind,strsplit(HR," ")),stringsAsFactors=F)
names(HR)=rownames(res)
HR <- t(HR)
####################################
library(rms)
library(survivalROC)
library(dplyr)
JIN_genes <- c('SERPINA1','CFHR3','PPP1CB','P2RX1','PLCB3','PLCB4','PIK3R6','GP6','PIK3R1','GP1BA','PLA2G4F')
LI_genes <- c('CASP6','NLRP7','NOD1','NLRP1','NLRP2')
GENG_genes <- c('NLRP7','CASP1','PLCG1','NLRP1')
exp_tumor <- exp[,colnames(exp)%in%clin_rna$sampleID[clin_rna$group=='Tumor']]
exp_BJ_tumor <- exp_tumor[rownames(exp_tumor)%in%c(bsrgenes[-10],gsub('_','-',bsrgenes[10]),JIN_genes,LI_genes,GENG_genes),]%>%t%>%as.data.frame()
exp_BJ_tumor$sampleID <- rownames(exp_BJ_tumor)
exp_BJ_tumor <- do.call(data.frame,exp_BJ_tumor)
colnames(exp_BJ_tumor)[23] <- 'JMJD7_PLA2G4B'
lxdata1 <- merge(exp_BJ_tumor,clin_tumor,by.x = 'sampleID')
lxdata1$time <- round(lxdata1$time,3)
lxdata1$Age <- ifelse(lxdata$Age=='young',1,2 )
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$Stage <- ifelse(lxdata$Stage=='I',1,ifelse(lxdata$Stage=='II',2 ,ifelse(lxdata$Stage=='III',3,ifelse(lxdata$Stage=='IV',4,0))))
lxdata1$Stage <- factor(lxdata1$Stage,levels=c(0,1,2,3,4),labels=c('UNKOWN','I','II','III','IV'))
#3种别的模型
# JIN = (0.188*lxdata1$SERPINA1) + (-0.153*lxdata1$CFHR3) + 
#   (0.571*lxdata1$PPP1CB) + (-0.831*lxdata1$P2RX1) + (0.422*lxdata1$PLCB3) + 
#  (-0.359*lxdata1$PLCB4) + (0.537*lxdata1$PIK3R6) + (-0.262*lxdata1$GP6) + (-1.170*lxdata1$PIK3R1) + (-0.257*lxdata1$GP1BA) + (-0.598*lxdata1$PLA2G4F)
lxdata1$JIN = (0.188*lxdata1$SERPINA1) + (-0.153*lxdata1$CFHR3) + (0.571*lxdata1$PPP1CB) + (-0.831*lxdata1$P2RX1) + (0.422*lxdata1$PLCB3) + 
     (-0.359*lxdata1$PLCB4) + (-1.170*lxdata1$PIK3R1) + (-0.257*lxdata1$GP1BA)
# LI = (0.0946) *lxdata1$CASP6 + (0.1573) *lxdata1$NLRP7 + (-0.124) *NOD1 + (-0.1627) *lxdata1$NLRP1 + (-0.0262) *lxdata1$NLRP2
lxdata1$LI = (0.0946) *lxdata1$CASP6 + (-0.124) *lxdata1$NOD1+ (-0.1627) *lxdata1$NLRP1 
# GENG = -0.08983*lxdata1$NLRP7-0.0120*lxdata1$CASP1-0.0645*lxdata1$PLCG1 -0.0598*lxdata1$NLRP1
lxdata1$GENG = -0.0120*lxdata1$CASP1-0.0645*lxdata1$PLCG1 -0.0598*lxdata1$NLRP1


lxdata1$Risk=0.09*lxdata1$F2-0.299*lxdata1$P2RX1-0.106*lxdata1$GUCY1A1+0.107*lxdata1$PLAUR-0.124*lxdata1$PRKCZ+0.093*lxdata1$F12+0.056*lxdata1$SERPINE1+
  0.137*lxdata1$COL1A2-0.022*lxdata1$C4BPA-0.031*lxdata1$JMJD7_PLA2G4B-0.017*lxdata1$CR2

lxdata1$F2_model=0.108 *lxdata1$F2
lxdata1$P2RX1_model=-0.213 *lxdata1$P2RX1
lxdata1$GUCY1A1_model=-0.192 *lxdata1$GUCY1A1
lxdata1$PLAUR_model=0.145 *lxdata1$PLAUR
lxdata1$PRKCZ_model=-0.194 *lxdata1$PRKCZ 
lxdata1$F12_model=0.133 *lxdata1$F12
lxdata1$SERPINE1_model=0.105 *lxdata1$SERPINE1
lxdata1$COL1A2_model=0.083 *lxdata1$COL1A2
lxdata1$C4BPA_model=-0.084 *lxdata1$C4BPA
lxdata1$JMJD7_PLA2G4B_model=-0.144 *lxdata1$JMJD7_PLA2G4B
lxdata1$CR2_model=-0.109 *lxdata1$CR2

troc.0= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$JIN,     
                    predict.time = 12, method="KM")
troc.01= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$LI,     
                    predict.time = 12, method="KM")
troc.02= survivalROC(Stime=lxdata1$time,  
                     status=lxdata1$status,   
                     marker = lxdata1$GENG,     
                     predict.time = 12, method="KM")
troc.1= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$Risk,     
                    predict.time = 12, method="KM")
troc.13= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$Risk,     
                    predict.time = 36, method="KM")
troc.15= survivalROC(Stime=lxdata1$time,  
                     status=lxdata1$status,   
                     marker = lxdata1$Risk,     
                     predict.time = 60, method="KM")
troc.2= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$F2_model,     
                    predict.time =12, method="KM")
troc.3= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$P2RX1_model,     
                    predict.time = 12, method="KM")
troc.4= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$GUCY1A1_model,     
                    predict.time = 12, method="KM")

troc.5= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$PLAUR_model,     
                    predict.time = 12, method="KM")
troc.6= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$PRKCZ_model,     
                    predict.time = 12, method="KM")
troc.7= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$F12_model,     
                    predict.time = 12, method="KM")
troc.8= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$SERPINE1_model,     
                    predict.time = 12, method="KM")
troc.9= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$COL1A2_model,     
                    predict.time = 12, method="KM")
troc.10= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$C4BPA_model,     
                    predict.time = 12, method="KM")
troc.11= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$JMJD7_PLA2G4B_model,     
                    predict.time = 12, method="KM")
troc.12= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$CR2_model,     
                    predict.time = 12, method="KM")



png("CRG/result/figs/03_rocmethods_bsr_TCGA.tiff",width = 1500,height = 1500,res=200)
plot(troc.1$FP, troc.1$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
title(main='TCGA',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.2$FP, troc.2$TP, type="l",lwd=4,col="blue", xlim=c(0,1), ylim=c(0,1))
lines(troc.3$FP, troc.3$TP, type="l",lwd=4,col="green", xlim=c(0,1), ylim=c(0,1))
lines(troc.4$FP, troc.4$TP, type="l",lwd=4,col="purple", xlim=c(0,1), ylim=c(0,1))
lines(troc.5$FP, troc.5$TP, type="l",lwd=4,col="orange", xlim=c(0,1), ylim=c(0,1))
lines(troc.6$FP, troc.6$TP, type="l",lwd=4,col="black", xlim=c(0,1), ylim=c(0,1))
lines(troc.7$FP, troc.7$TP, type="l",lwd=4,col="yellow", xlim=c(0,1), ylim=c(0,1))
lines(troc.8$FP, troc.8$TP, type="l",lwd=4,col="pink", xlim=c(0,1), ylim=c(0,1))
lines(troc.9$FP, troc.9$TP, type="l",lwd=4,col="cyan", xlim=c(0,1), ylim=c(0,1))
lines(troc.10$FP, troc.10$TP, type="l",lwd=4,col="magenta", xlim=c(0,1), ylim=c(0,1))
lines(troc.11$FP, troc.11$TP, type="l",lwd=4,col="seagreen4", xlim=c(0,1), ylim=c(0,1))
lines(troc.12$FP, troc.12$TP, type="l",lwd=4,col="red4", xlim=c(0,1), ylim=c(0,1))
legend(0.55,0.58,c(paste('Risk=',round(troc.1$AUC,3)),
                  paste('F2=',round(troc.2$AUC,3)),
                  paste('P2RX1=',round(troc.3$AUC,3)),
                  paste('GUCY1A1=',round(troc.4$AUC,3)),
                  paste('PLAUR=',round(troc.5$AUC,3)),
                  paste('PRKCZ=',round(troc.6$AUC,3)),
                  paste('LF12=',round(troc.7$AUC,3)),
                  paste('SERPINE1=',round(troc.8$AUC,3)),
                  paste('COL1A2=',round(troc.9$AUC,3)),
                  paste('C4BPA=',round(troc.10$AUC,3)),
                  paste('JMJD7-PLA2G4B=',round(troc.11$AUC,3)),
                  paste('CR2=',round(troc.12$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','blue',"green","purple",'orange','black',"yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()
###############################################################################
library(pec)
library(survival)
dat1 <- lxdata1

cox1 <- coxph(Surv(time, status) ~ Risk, data = dat1,x=T,y=T)
cox2 <- coxph(Surv(time, status) ~ F2, data = dat1,x=T,y=T)
cox3 <- coxph(Surv(time, status) ~ P2RX1, data = dat1,x=T,y=T)
cox4 <- coxph(Surv(time, status) ~ GUCY1A1, data = dat1,x=T,y=T)
cox5 <- coxph(Surv(time, status) ~ PLAUR, data = dat1,x=T,y=T)
cox6 <- coxph(Surv(time, status) ~ PRKCZ, data = dat1,x=T,y=T)
cox7 <- coxph(Surv(time, status) ~ F12, data = dat1,x=T,y=T)
cox8 <- coxph(Surv(time, status) ~ SERPINE1, data = dat1,x=T,y=T)
cox9 <- coxph(Surv(time, status) ~ COL1A2, data = dat1,x=T,y=T)
cox10 <- coxph(Surv(time, status) ~ C4BPA, data = dat1,x=T,y=T)
cox11 <- coxph(Surv(time, status) ~ JMJD7_PLA2G4B, data = dat1,x=T,y=T)
cox12 <- coxph(Surv(time, status) ~ CR2, data = dat1,x=T,y=T)


ApparrentCindex  <- pec::cindex(list("Risk"=cox1,
                                     "F2"=cox2,
                                     "P2RX1"=cox3,
                                     "GUCY1A1"=cox4,
                                     "PLAUR"=cox5,
                                     "PRKCZ"=cox6,
                                     "LF12"=cox7,
                                     "SERPINE1"=cox8,
                                     "COL1A2"=cox9,
                                     "C4BPA"=cox10,
                                     "JMJD7-PLA2G4B"=cox11,
                                     "CR2"=cox12),
                                formula=Surv(time,status)~Risk,
                                data=dat1,
                                eval.times=seq(1,242,20))

png("CRG/result/figs/03_cindex.tiff",width = 1500,height = 1500,res=200)
plot(ApparrentCindex,lwd=4,cex.axis=2,legend.x=190,legend.y=1,legend.cex=1,col=c('red',"purple3",'blue',"green","purple",'orange','black',"yellow","pink","cyan","magenta","seagreen4","red4"))
dev.off()


###################################################################################

png("CRG/result/figs/03_roccompare_bsr_TCGA.tiff",width = 1500,height = 1500,res=200)
plot(troc.1$FP, troc.1$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
title(main='TCGA',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.0$FP, troc.0$TP, type="l",lwd=4,col="purple3", xlim=c(0,1), ylim=c(0,1))
lines(troc.01$FP, troc.01$TP, type="l",lwd=4,col="blue", xlim=c(0,1), ylim=c(0,1))
lines(troc.02$FP, troc.02$TP, type="l",lwd=4,col="green", xlim=c(0,1), ylim=c(0,1))
legend(0.6,0.2,c(paste('Risk=',round(troc.1$AUC,3)),
                  paste('JIN=',round(troc.0$AUC,3)),
                  paste('LI=',round(troc.01$AUC,3)),
                  paste('GENG=',round(troc.02$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red',"purple3",'blue',"green","purple",'orange','black',"yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()


png("CRG/result/figs/03_rocyear_risk_TCGA.tiff",width = 1500,height = 1500,res=200)
plot(troc.1$FP, troc.1$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
title(main='TCGA',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.13$FP, troc.13$TP, type="l",lwd=4,col="orange", xlim=c(0,1), ylim=c(0,1))
lines(troc.15$FP, troc.15$TP, type="l",lwd=4,col="blue", xlim=c(0,1), ylim=c(0,1))
legend(0.55,0.2,c(paste('AUC of 1-year =',round(troc.1$AUC,3)),
                  paste('AUC of 3-year =',round(troc.13$AUC,3)),
                  paste('AUC of 5-year =',round(troc.15$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()

################################################
lxdata1$Group <- ifelse(lxdata1$Risk>median(lxdata1$Risk),'High risk','Low risk')
png("CRG/result/figs/03_OSKM_risk_TCGA.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(time, status) ~ Group, data=lxdata1),
                surv.median.line = "hv",size=3,
                title = "TCGA",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("CRG/result/figs/03_RFSKM_risk_TCGA.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(RFS, RFS.E) ~ Group, data=lxdata1),size=3,
                surv.median.line = "hv",
                title = "TCGA",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("red", "blue",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of RFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()



png("CRG/result/figs/03_OSKM_risk_TCGAlater.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(time, status) ~ Group, data=lxdata1[lxdata1$Stage%in%c('III','IV'),]),size=3,
                surv.median.line = "hv",
                title = "TCGA Stage III~IV",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("pink", "purple",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("CRG/result/figs/03_RFSKM_risk_TCGAlater.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(RFS, RFS.E) ~ Group, data=lxdata1[lxdata1$Stage%in%c('III','IV'),]),size=3,
                surv.median.line = "hv",
                title = "TCGA Stage III~IV",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("lightblue", "yellowgreen",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of RFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

###############################################
library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
#风险评分box图
cb <- ggboxplot(lxdata1,x = 'Status',y = 'Risk',fill='Status')+theme_bw()+ylab('Riskscore')+xlab('')+ggtitle('TCGA')+ylim(-2,3)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_signif(size=1,textsize=12,comparisons = list(c("Dead", "Alive")),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
    axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "none",
    panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('CRG/result/figs/03_riskbox.tiff',cb,width = 5,height = 5,dpi=200)

#点图
pd <- lxdata1[order(lxdata1$Risk,decreasing = F),]
pd$Sample <- 1:nrow(pd)
sp <- ggplot(pd,aes(Sample,time,color=Status))+geom_point()+geom_vline(aes(xintercept = 251),lty=2,lwd=1)+
  theme_bw()+scale_colour_manual(values = c('blue','red'))+ylab('Survival month')+xlab('')+
  theme(axis.text.x = element_text(size=16,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
rp <- ggplot(pd,aes(Sample,Risk,color=Group))+geom_point()+geom_vline(aes(xintercept = 251),lty=2,lwd=1)+
  theme_bw()+scale_colour_manual(values = c('red','blue'))+ylab('Risk score')+xlab('')+ggtitle('TCGA')+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.justification=c(0,1),legend.position=c(0.01,0.99),
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('CRG/result/figs/03_survpoint.tiff',sp,width = 5,height =3,dpi=1000)
ggsave('CRG/result/figs/03_riskpoint.tiff',rp,width = 5,height = 3,dpi=1000)

library(pheatmap)
bsrgenes[10] <- gsub('_','-',bsrgenes[10])
gene <- bsrgenes
dat <- pd[,colnames(pd)%in%gene]
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

save_pheatmap_pdf(p, "CRG/result/figs/03_pheat.pdf")
save(pd,gene,file = 'CRG/extdata/tcga_plot.Rdata')
load('CRG/extdata/tcga_plot.Rdata')


