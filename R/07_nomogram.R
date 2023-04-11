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
exp_tumor <- exp[,colnames(exp)%in%clin_tumor$sampleID]
exp_BJ_tumor <- exp_tumor[rownames(exp_tumor)%in%c(bsrgenes[-10],gsub('_','-',bsrgenes[10])),]%>%t%>%as.data.frame()
exp_BJ_tumor$sampleID <- rownames(exp_BJ_tumor)
exp_BJ_tumor <- do.call(data.frame,exp_BJ_tumor)
colnames(exp_BJ_tumor)[11] <- 'JMJD7_PLA2G4B'
lxdata <- merge(exp_BJ_tumor,clin_tumor,by.x = 'sampleID')
lxdata1 <- merge(exp_BJ_tumor,clin_tumor,by.x = 'sampleID')
lxdata1$OS <- round(lxdata1$OS,3)
lxdata1$Age <- ifelse(lxdata$Age=='young',1,2 )
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$Stage <- ifelse(lxdata$Stage=='I',1,ifelse(lxdata$Stage=='II',2 ,ifelse(lxdata$Stage=='III',3,ifelse(lxdata$Stage=='IV',4,0))))
lxdata1$Stage <- factor(lxdata1$Stage,levels=c(0,1,2,3,4),labels=c('UNKOWN','I','II','III','IV'))
mulcox<- coxph(Surv(OS, OS.E) ~ CR2+F2+F12+GUCY1A1+ITGB1+P2RX1+SERPINE1+PLAUR+PRKCZ+C4BPA+JMJD7_PLA2G4B, data = lxdata1 );mulcox
lxdata1$Risk=predict(mulcox,lxdata1)
lxdata1$Group <- ifelse(lxdata1$Risk>median(lxdata1$Risk),'High risk','Low risk')
######################riskscore unicox
library(survival)
library(survminer)

covariates <- c('Gender','Age','Stage','Smoking',"Tumor" ,  "Lymph_Node", "Metastasis", 'Risk')
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS, OS.E)~', x)))
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
# Fig6B
mulcox1<- coxph(Surv(OS, OS.E) ~ Age+Tumor+Lymph_Node+Risk, data = lxdata1 );mulcox1

forest <- ggforest(mulcox1,  #coxph得到的Cox回归结果
                   data = lxdata1,  #数据集
                   main = 'Hazard ratio of all data',  #标题
                   cpositions = c(0.05, 0.15, 0.35),  #前三列数值的距离
                   fontsize = 2, #字体大小
                   refLabel = 'Reference', #相对变量的标签
                   noDigits = 3 #HR、95%CI的小数位
)
ggsave(filename="figures/06B_muti_forest.png",forest, width=20, height=14)

# Fig6C
library(pec)
library(survival)
dat <- lxdata1[,c('OS', 'OS.E', 'Age','Stage','Risk')]
dat <- dat[dat$Stage!='UNKOWN',]
dat$Stage <- as.character(dat$Stage)
dat$time <- dat$OS *30
library(rms)
dd <- datadist(dat)
options(datadist='dd')

m2<- psm(Surv(OS, OS.E) ~ Age+Stage+Risk, data = dat,dist='lognormal' )
med <- Quantile(m2)
surv <- Survival(m2)
nomo <- nomogram(m2,fun=list(function(x)surv(365,x),function(x)surv(1095,x),function(x)surv(1825,x)), funlabel=c("1-year survival probability","3-year survival probability","5-year survival probability"))
png("figures/06C_nomo.png",width =3600,height = 2000,res=200)
plot(nomo,xfrac=0.2, lwd=4,cex.axis=0.9,cex.lab=2,font.lab=4,cex.var=2)
dev.off()

cindex <- rcorrcens(Surv(OS, OS.E) ~ predict(m2), data = dat)

#0.722

# Fig6D
f2 <- psm(Surv(time, OS.E) ~ Age+Stage+Risk, data = dat,x=T,y=T,dist='lognormal' )
call <- calibrate(f2,cmethod = 'KM',method='boot',u=365,m=164,B=1000)
png("figures/06E_cal_1y.png",width = 1500,height = 1500,res=200)
plot(call,lty=1,conf.int=T,errbar.col='blue',col='red',xlim=c(0,1),ylim=c(0,1), lwd=4,cex.axis=1.6,ann = F,subtitles=F)
title(xlab='Nomogram-predicted OS of 1-year',ylab = 'Actual 1-year OS',cex.lab=1.6,font.lab=2)
dev.off()

call2 <- calibrate(f2,cmethod = 'KM',method='boot',u=1095,m=164,B=1000)
png("figures/06E_cal_3y.png",width = 1500,height = 1500,res=200)
plot(call2,lty=1,conf.int=T,errbar.col='blue',col='red',xlim=c(0,1),ylim=c(0,1), lwd=4,cex.axis=1.6,ann = F,subtitles=F)
title(xlab='Nomogram-predicted OS of 3-year',ylab = 'Actual 3-year OS',cex.lab=1.6,font.lab=2)
dev.off()

call3 <- calibrate(f2,cmethod = 'KM',method='boot',u=1825,m=164,B=1000)
png("figures/06E_cal_5y.png",width = 1500,height = 1500,res=200)
plot(call3,lty=1,conf.int=T,errbar.col='blue',col='red',xlim=c(0,1),ylim=c(0,1), lwd=4,cex.axis=1.6,ann = F,subtitles=F)
title(xlab='Nomogram-predicted OS of 5-year',ylab = 'Actual 5-year OS',cex.lab=1.6,font.lab=2)
dev.off()

# Fig6D
library(pec)
library(survival)
dat1 <- na.omit(dat)
cox1 <- coxph(Surv(OS, OS.E) ~ Age+Stage+Risk, data = dat1,x=T,y=T)
cox2 <- coxph(Surv(OS, OS.E) ~ Age, data = dat1,x=T,y=T)
cox3 <- coxph(Surv(OS, OS.E) ~ Stage, data = dat1,x=T,y=T)
cox4 <- coxph(Surv(OS, OS.E) ~ Risk, data = dat1,x=T,y=T)
ApparrentCindex  <- pec::cindex(list("Age"=cox2,
                                     "Stage"=cox3,"Risk"=cox4,"Nomogram"=cox1),
                                formula=Surv(OS, OS.E)~Age+Stage+Risk,
                                data=dat1,
                                eval.times=seq(1,242,20))

png("figures/06D_cindex.tiff",width = 1500,height = 1500,res=200)
plot(ApparrentCindex,lwd=4,cex.axis=2,legend.x=190,legend.y=1,legend.cex=1,col=c('red',"purple3",'blue',"green","purple",'orange','black',"yellow","pink","cyan","magenta","seagreen4","red4"))
dev.off()

# Fig6F
##########AUC combinied riskscore age stage
####1-year
lxdata1$age <- ifelse(lxdata$Age=='young',1,2)
lxdata1$stage<- ifelse(lxdata$Stage=='I',1, ifelse(lxdata$Stage=='II',2, ifelse(lxdata$Stage=='III',3, ifelse(lxdata$Stage=='IV',4,0))))

model1=coxph(Surv(OS, OS.E) ~age, data = lxdata1)
model2=coxph(Surv(OS, OS.E) ~ stage, data = lxdata1)
model3=coxph(Surv(OS, OS.E) ~ age+stage+Risk, data = lxdata1 )
lxdata1$riskscore1 <- predict(model1,lxdata1)
lxdata1$riskscore2 <- predict(model2,lxdata1)
lxdata1$riskscore3 <- predict(model3,lxdata1)

library(survivalROC)
#12,36,60
time=60
troc= survivalROC(Stime=lxdata1$OS,  
                  status=lxdata1$OS.E,   
                  marker = lxdata1$riskscore3,     
                  predict.time = time, method="KM")
troc.1= survivalROC(Stime=lxdata1$OS,  
                    status=lxdata1$OS.E,  
                    marker = lxdata1$Risk,     
                    predict.time = time, method="KM")
troc.2= survivalROC(Stime=lxdata1$OS,  
                    status=lxdata1$OS.E,  
                    marker = lxdata1$riskscore1,     
                    predict.time = time, method="KM")
troc.3= survivalROC(Stime=lxdata1$OS,  
                    status=lxdata1$OS.E, 
                    marker = lxdata1$riskscore2,     
                    predict.time = time, method="KM")

png("figures/06F_1-yearroc_combined.png",width = 1500,height = 1500,res=200)
png("figures/06F_3-yearroc_combined.png",width = 1500,height = 1500,res=200)
png("figures/06F_5-yearroc_combined.png",width = 1500,height = 1500,res=200)

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
write.xlsx(lxdata1[,c('Gender',	'Age',	'Stage',
                      'Tumor',	'Lymph_Node',	'Metastasis',
                      'OS',	'OS.E',	'RFS',	'RFS.E',	'Risk',	'Group')],file = 'files/Fig6.xlsx',sheetName = 'Fig6_clin',row.names = F,append = T)
