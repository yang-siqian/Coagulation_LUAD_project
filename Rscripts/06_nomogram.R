library(dplyr)
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

bsrgenes <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7_PLA2G4B','CR2')
JIN_genes <- c('SERPINA1','CFHR3','PPP1CB','P2RX1','PLCB3','PLCB4','PIK3R6','GP6','PIK3R1','GP1BA','PLA2G4F')
LI_genes <- c('CASP6','NLRP7','NOD1','NLRP1','NLRP2')
GENG_genes <- c('NLRP7','CASP1','PLCG1','NLRP1')
exp_tumor <- exp[,colnames(exp)%in%clin_rna$sampleID[clin_rna$group=='Tumor']]
exp_BJ_tumor <- exp_tumor[rownames(exp_tumor)%in%c(bsrgenes[-10],gsub('_','-',bsrgenes[10]),JIN_genes,LI_genes,GENG_genes),]%>%t%>%as.data.frame()
exp_BJ_tumor$sampleID <- rownames(exp_BJ_tumor)
exp_BJ_tumor <- do.call(data.frame,exp_BJ_tumor)
colnames(exp_BJ_tumor)[23] <- 'JMJD7_PLA2G4B'
lxdata <- merge(exp_BJ_tumor,clin_tumor,by.x = 'sampleID')
lxdata1 <- merge(exp_BJ_tumor,clin_tumor,by.x = 'sampleID')
lxdata1$time <- round(lxdata1$time,3)
lxdata1$Age <- ifelse(lxdata$Age=='young',1,2 )
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$Stage <- ifelse(lxdata$Stage=='I',1,ifelse(lxdata$Stage=='II',2 ,ifelse(lxdata$Stage=='III',3,ifelse(lxdata$Stage=='IV',4,0))))
lxdata1$Stage <- factor(lxdata1$Stage,levels=c(0,1,2,3,4),labels=c('UNKOWN','I','II','III','IV'))
lxdata1$Risk=0.09*lxdata1$F2-0.299*lxdata1$P2RX1-0.106*lxdata1$GUCY1A1+0.107*lxdata1$PLAUR-0.124*lxdata1$PRKCZ+0.093*lxdata1$F12+0.056*lxdata1$SERPINE1+
  0.137*lxdata1$COL1A2-0.022*lxdata1$C4BPA-0.031*lxdata1$JMJD7_PLA2G4B-0.017*lxdata1$CR2
lxdata1$Group <- ifelse(lxdata1$Risk>median(lxdata1$Risk),'High risk','Low risk')
######################riskscore unicox
library(survival)
library(survminer)

covariates <- c('Gender','Age','Stage','Smoking',"Tumor" ,  "Lymph_Node", "Metastasis", 'Risk')
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data =  lxdata1)})
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
###########riskscore  multicox
mulcox<- coxph(Surv(time, status) ~ Age+Tumor+Lymph_Node+Metastasis+Risk, data = lxdata1 );mulcox

forest <- ggforest(mulcox,  #coxph得到的Cox回归结果
                   data = lxdata1,  #数据集
                   main = 'Hazard ratio of all data',  #标题
                   cpositions = c(0.05, 0.15, 0.35),  #前三列数值的距离
                   fontsize = 2, #字体大小
                   refLabel = 'Reference', #相对变量的标签
                   noDigits = 3 #HR、95%CI的小数位
)
ggsave(filename="CRG/result/figs/08_muti_forest.png",forest, width=20, height=14)


library(pec)
library(survival)
dat <- lxdata1[,c('time', 'status', 'Age','Stage','Risk')]
dat$time <- dat$time *30
library(rms)
dd <- datadist(dat)
options(datadist='dd')

m2<- psm(Surv(time, status) ~ Age+Stage+Risk, data = dat,dist='lognormal' )
med <- Quantile(m2)
surv <- Survival(m2)
nomo <- nomogram(m2,fun=list(function(x)surv(365,x),function(x)surv(1095,x),function(x)surv(1825,x)), funlabel=c("1-year survival probability","3-year survival probability","5-year survival probability"))
png("CRG/result/figs/08_nomo.tiff",width = 2500,height = 2000,res=200)
plot(nomo,xfrac=0.6, lwd=4,cex.axis=2,cex.lab=2,font.lab=4,cex.var=2,)
dev.off()

cindex <- rcorrcens(Surv(time, status) ~ predict(m2), data = dat)

#0.722

f2 <- psm(Surv(time, status) ~ Age+Stage+Risk, data = dat,x=T,y=T,dist='lognormal' )
call <- calibrate(f2,cmethod = 'KM',method='boot',u=365,m=164,B=1000)
png("CRG/result/figs/08_cal_1y.tiff",width = 1500,height = 1500,res=200)
plot(call,lty=1,conf.int=T,errbar.col='blue',col='red',xlim=c(0,1),ylim=c(0,1), lwd=4,cex.axis=1.6,ann = F,subtitles=F)
title(xlab='Nomogram-predicted OS of 1-year',ylab = 'Actual 1-year OS',cex.lab=1.6,font.lab=2)
dev.off()

call2 <- calibrate(f2,cmethod = 'KM',method='boot',u=1095,m=164,B=1000)
png("CRG/result/figs/08_cal_3y.tiff",width = 1500,height = 1500,res=200)
plot(call2,lty=1,conf.int=T,errbar.col='blue',col='red',xlim=c(0,1),ylim=c(0,1), lwd=4,cex.axis=1.6,ann = F,subtitles=F)
title(xlab='Nomogram-predicted OS of 3-year',ylab = 'Actual 3-year OS',cex.lab=1.6,font.lab=2)
dev.off()

call3 <- calibrate(f2,cmethod = 'KM',method='boot',u=1825,m=164,B=1000)
png("CRG/result/figs/08_cal_5y.tiff",width = 1500,height = 1500,res=200)
plot(call3,lty=1,conf.int=T,errbar.col='blue',col='red',xlim=c(0,1),ylim=c(0,1), lwd=4,cex.axis=1.6,ann = F,subtitles=F)
title(xlab='Nomogram-predicted OS of 5-year',ylab = 'Actual 5-year OS',cex.lab=1.6,font.lab=2)
dev.off()

library(pec)
library(survival)
dat1 <- na.omit(dat)
dat1$time <- dat1$time/30
cox1 <- coxph(Surv(time, status) ~ Age+Stage+Risk, data = dat1,x=T,y=T)
cox2 <- coxph(Surv(time, status) ~ Age, data = dat1,x=T,y=T)
cox3 <- coxph(Surv(time, status) ~ Stage, data = dat1,x=T,y=T)
cox4 <- coxph(Surv(time, status) ~ Risk, data = dat1,x=T,y=T)
ApparrentCindex  <- pec::cindex(list("Age"=cox2,
                                     "Stage"=cox3,"Risk"=cox4,"Nomogram"=cox1),
                                formula=Surv(time,status)~Age+Stage+Risk,
                                data=dat1,
                                eval.times=seq(1,242,20))

png("CRG/result/figs/08_cindex.tiff",width = 1500,height = 1500,res=200)
plot(ApparrentCindex,lwd=4,cex.axis=2,legend.x=190,legend.y=1,legend.cex=1,col=c('red',"purple3",'blue',"green","purple",'orange','black',"yellow","pink","cyan","magenta","seagreen4","red4"))
dev.off()


##########AUC combinied riskscore age stage
####1-year
coxph(Surv(time, status) ~ age, data = lxdata1)
coxph(Surv(time, status) ~ stage, data = lxdata1)
coxph(Surv(time, status) ~ age+stage+Risk, data = lxdata1 )
lxdata1$age <- ifelse(lxdata$Age=='young',1,2)
lxdata1$stage<- ifelse(lxdata$Stage=='I',1, ifelse(lxdata$Stage=='II',2, ifelse(lxdata$Stage=='III',3, ifelse(lxdata$Stage=='IV',4,0))))
lxdata1$riskscore1 <- -0.3797 *as.numeric(lxdata1$age)
lxdata1$riskscore2 <- 0.51381*as.numeric(lxdata1$stage)
lxdata1$riskscore3 <- 0.42552  *lxdata1$age+0.36421 *lxdata1$stage+0.91900 *lxdata1$Risk

library(survivalROC)
#12,36,60
time=60
troc= survivalROC(Stime=lxdata1$time,  
                  status=lxdata1$status,   
                  marker = lxdata1$riskscore3,     
                  predict.time = time, method="KM")
troc.1= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$Risk,     
                    predict.time = time, method="KM")
troc.2= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$riskscore1,     
                    predict.time = time, method="KM")
troc.3= survivalROC(Stime=lxdata1$time,  
                    status=lxdata1$status,   
                    marker = lxdata1$riskscore2,     
                    predict.time = time, method="KM")

png("CRG/result/figs/08_1-yearroc_combined.tiff",width = 1500,height = 1500,res=200)
png("CRG/result/figs/08_3-yearroc_combined.tiff",width = 1500,height = 1500,res=200)
png("CRG/result/figs/08_5-yearroc_combined.tiff",width = 1500,height = 1500,res=200)

plot(troc$FP, troc$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)

title(main='1-Year',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
title(main='3-Year',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
title(main='5-Year',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)

abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.1$FP, troc.1$TP, type="l",lwd=4,col="blue", xlim=c(0,1), ylim=c(0,1))
lines(troc.2$FP, troc.2$TP, type="l",lwd=4,col="green", xlim=c(0,1), ylim=c(0,1))
lines(troc.3$FP, troc.3$TP, type="l",lwd=4,col="purple", xlim=c(0,1), ylim=c(0,1))
legend(0.58,0.3,c(paste('Nomogram=',round(troc$AUC,3)),
                  paste('Risk=',round(troc.1$AUC,3)),
                  paste('Age=',round(troc.2$AUC,3)),
                  paste('Stage=',round(troc.3$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','blue',"green","purple"),
       bty='n',seg.len = 1.5,cex = 1.2)
dev.off()

