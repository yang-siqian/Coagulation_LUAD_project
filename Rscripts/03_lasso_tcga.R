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
# gf=bitr(rownames(exp), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
gf <- read.csv('files/gf.csv',header=T)
exp <- merge(exp,gf, by.y='ENTREZID')
exp<-exp[!duplicated(exp$SYMBOL),]
rownames(exp)<-exp$SYMBOL[!duplicated(exp$SYMBOL)]
exp <- exp[,-which(colnames(exp)%in%c( "ENSEMBL", "ENTREZID","SYMBOL"))]
exp_scale <- scale(exp)%>%as.data.frame()
group <- ifelse(as.numeric(substr(colnames(exp),14,15))<10,'Tumor','Normal')
exp_CRGs <- exp_scale[rownames(exp_scale)%in%CRGs$gene,]

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
survdata<- clin_tumor[,c(1,20,21)]
CRGs_surv_dat <- merge(survdata,exp_CRGs_tumor,by.x='sampleID')
library(glmnet)

X.train <- as.matrix(CRGs_surv_dat[,4:ncol(CRGs_surv_dat)])
Y.train <- as.matrix(cbind(time=CRGs_surv_dat$OS,status=CRGs_surv_dat$OS.E))
lgs <- c()
for(i in 1:100){
  set.seed(i)
  fit.cv <- cv.glmnet(X.train,Y.train, family = "cox", nfolds = 10,alpha=1,keep=TRUE)
  idmin = match(fit.cv$lambda.min, fit.cv$lambda)
  lsocoef<-as.numeric(coef(fit.cv ,s="lambda.min"))
  lasso_genes <- colnames(X.train)[lsocoef!=0]
  lgs <- c(lgs,lasso_genes)
}
table(lgs)
lasso_genes=names(table(lgs))[table(lgs)==100]

lgs_freq <- as.data.frame(t(table(lgs)))
lgs_freq <- lgs_freq[,-1]
lgs_freq <- lgs_freq[order(lgs_freq$Freq,decreasing = T),]
lgs_freq$lgs <- factor(lgs_freq$lgs,levels=lgs_freq$lgs)

discrete_palettes <- rep('skyblue',28)
p<- ggplot(lgs_freq,aes(x=lgs,y=Freq,fill = lgs))+geom_bar(stat = "identity",color=c(rep('red',6),rep('skyblue',22)))+coord_flip()+theme_bw()+ylab('Frequency')+xlab('')+
  theme(legend.position = 'none',axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),text=element_text(size=18,face='bold'),panel.grid = element_blank())+scale_fill_manual(values = discrete_palettes)
  
ggsave('figures/02A_LASSO_frequency.png',p,width = 10,height = 12,dpi=200)
library(xlsx)
write.xlsx(lgs_freq,'figures/02_data.xlsx',row.names = F)


fit.cv$lambda.min
#0.04818326

# Fig2B
png("figures/02B_lassocoef.png",width = 1500,height = 1500,res=200)
plot(fit.cv, lwd=3,cex.axis=1.5,ann = F)
title(xlab="log(λ)", ylab="Partial likelihood Deviance",cex.lab=1.5,font.lab=2)
dev.off()



model <- glmnet(X.train,Y.train, family = "cox", alpha=1)

# Fig2A
png("figures/02A_coefficient.png",width = 1500,height = 1500,res=200)
plot(model, lwd=3,cex.axis=1.5,ann = F)
title(xlab="log(λ)", ylab="Coefficients",cex.lab=1.5,font.lab=2)
dev.off()


dmf <- read.xlsx(file = 'files/Fig2.xlsx',sheetName = 'Fig2D')
dmf <- dmf[1:11,1:6]
dmf <- dmf[order(dmf$coef,decreasing = T),]
dmf$gene <- factor(dmf$gene,levels=dmf$gene)
# Fig2C
library(ggalt)
dp <- ggplot(dmf,aes(y = gene, x=coef,xend=0))+
  geom_segment(aes(x=coef,xend=0,y=gene,yend=gene),color='grey',size=2)+
  geom_dumbbell(size_x=6, size_xend = 0,color="grey",colour_x = "black")+
  labs(x='Coefficient', y=NULL) +geom_vline(aes(xintercept = 0),lty=2,lwd=1)+
  theme_bw() + theme(axis.text.y = element_text(size=16,face='bold',color='black'),
                     axis.text.x = element_text(size=15,face='bold',color='black'),
                     text=element_text(size=16,face='bold'),axis.line = element_line(size=1),panel.grid = element_blank(),panel.border = element_rect(size=1.5))

ggsave('figures/02C_dot_coef_bsr.png',dp,width = 6,height = 6,dpi=1000)




####################################################
# Fig2D
lasso_genes <- c("C4BPA" , "COL1A2" ,"F2" ,"P2RX1" , "PLAUR" )
hubdat <-  exp_CRGs_tumor[,c('sampleID',lasso_genes)]
modeldat<- merge(hubdat,clin_tumor,by.x = 'sampleID')
lxdata <- modeldat
lxdata1 <- merge(exp_CRGs_tumor,clin_tumor,by.x = 'sampleID')
lxdata1$OS <- round(lxdata1$OS,3)
lxdata1$Age <- ifelse(lxdata$Age=='young',1,2 )
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$Stage <- ifelse(lxdata$Stage=='I',1,ifelse(lxdata$Stage=='II',2 ,ifelse(lxdata$Stage=='III',3,ifelse(lxdata$Stage=='IV',4,0))))
lxdata1$Stage <- factor(lxdata1$Stage,levels=c(0,1,2,3,4),labels=c('UNKOWN','I','II','III','IV'))

library(survival)
library(survminer)
library(rms)
library(survivalROC)
library(dplyr)
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
HR=gsub("-"," ",HR)
HR=as.data.frame(do.call(cbind,strsplit(HR," ")),stringsAsFactors=F)
names(HR)=rownames(res)
HR <- t(HR)
###############################################
library(My.stepwise)
unigenes <- rownames(res)[as.numeric(res$p.value)<0.05]
# unigenes <-c("C4BPA" , "COL1A2" ,"F2" ,"P2RX1" , "PLAUR" )
bsrdata <- CRGs_surv_dat[,c('OS','OS.E',unigenes)]
# colnames(bsrdata)[13] <- gsub('-','_',colnames(bsrdata)[13])
# unigenes[11] <- gsub('-','_',unigenes[11])
rownames(bsrdata) <- CRGs_surv_dat$sampleID

library(StepReg)
formula = Surv(OS, OS.E) ~ . - OS.E
bsr1 <- stepwiseCox(formula,
                    bsrdata,
                    include=NULL,
                    selection=c("score"),
                    select="AIC",
                    method=c("efron"),
                    sle=0.05,
                    sls=0.1,
                    weights=NULL,
                    best=1)

# library(ggpubr)
# coef <- round(as.numeric(bsr1$`Coefficients of the Selected Variables`$coef),3)
# logp <- round(-log10(as.numeric(bsr1$`Coefficients of the Selected Variables`$`Pr(>|z|)`)),1)
# gene <- bsr1$`Coefficients of the Selected Variables`$Variable
# gene[10] <- gsub('_','-',gene[10] )
# dmf <-cbind(gene,coef,logp)
# dmf <- as.data.frame(dmf)
# dmf <- dmf[order(dmf$logp,decreasing = T),]
# dmf$gene <- factor(dmf$gene,levels=dmf$gene) 
# pp <- ggbarplot(dmf,x='gene',y='logp',fill='grey')+coord_flip()+
#   labs(x='',y='-log10(p-value)', y=NULL) +geom_vline(aes(xintercept = 0),lty=2,lwd=1)+
#   theme_bw() + theme(axis.text.y = element_text(size=16,face='bold',color='black'),
#                      axis.text.x = element_text(size=15,face='bold',color='black'),
#                      text=element_text(size=16,face='bold'),axis.line = element_line(size=1),panel.grid = element_blank(),panel.border = element_rect(size=1.5))
# ggsave('figures/02C_bar_pvalue_bsr.tiff',pp,width = 6,height = 6,dpi=1000)
# 
# library(xlsx)
# write.xlsx(bsr1$`Coefficients of the Selected Variables`,file='files/Fig2.xlsx',sheetName = 'Fig2C',row.names = F)
###########riskscore  multicox
# Fig2D
# mulcox<- coxph(Surv(OS, OS.E) ~ CR2+F2+F12+GUCY1A1+ITGB1+P2RX1+SERPINE1+PLAUR+PRKCZ+C4BPA+JMJD7_PLA2G4B, data = lxdata1 );mulcox

mulcox<- coxph(Surv(OS, OS.E) ~ COL1A2+F2+P2RX1+PLAUR+C4BPA, data = lxdata1 );mulcox


forest <- ggforest(mulcox,  #coxph得到的Cox回归结果
                   data = lxdata1,  #数据集
                   main = 'Hazard ratio of all data',  #标题
                   cpositions = c(0.05, 0.15, 0.35),  #前三列数值的距离
                   fontsize = 2, #字体大小
                   refLabel = 'Reference', #相对变量的标签
                   noDigits = 3 #HR、95%CI的小数位
)
ggsave(filename="figures/02C_muti_forest.png",forest, width=28, height=14)

dmf <- read.xlsx(file = 'figures/02_data.xlsx',sheetName = 'Sheet2')
dmf <- dmf[order(dmf$coef,decreasing = T),]
dmf$gene <- factor(dmf$gene,levels=dmf$gene)
# Fig2C
library(ggalt)
dp <- ggplot(dmf,aes(y = gene, x=coef,xend=0))+
  geom_segment(aes(x=coef,xend=0,y=gene,yend=gene),color='grey',size=2)+
  geom_dumbbell(size_x=6, size_xend = 0,color="grey",colour_x = "black")+
  labs(x='Coefficient', y=NULL) +geom_vline(aes(xintercept = 0),lty=2,lwd=1)+
  theme_bw() + theme(axis.text.y = element_text(size=16,face='bold',color='black'),
                     axis.text.x = element_text(size=15,face='bold',color='black'),
                     text=element_text(size=16,face='bold'),axis.line = element_line(size=1),panel.grid = element_blank(),panel.border = element_rect(size=1.5))

ggsave('figures/02C_dot_coef_bsr.png',dp,width = 6,height = 6,dpi=1000)


# Fig2F
library(rms)
library(survivalROC)
library(dplyr)
bsrgenes <- unigenes
exp_tumor <- exp[,colnames(exp)%in%clin_tcga$sampleID[clin_tcga$group=='Tumor']]
exp_BJ_tumor <- exp_tumor[rownames(exp_tumor)%in%bsrgenes,]%>%t%>%as.data.frame()
exp_BJ_tumor$sampleID <- rownames(exp_BJ_tumor)
exp_BJ_tumor <- do.call(data.frame,exp_BJ_tumor)
lxdata1 <- merge(exp_BJ_tumor,clin_tumor,by.x = 'sampleID')
lxdata1$OS <- round(lxdata1$OS,3)
lxdata1$Age <- ifelse(lxdata$Age=='young',1,2 )
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$Stage <- ifelse(lxdata$Stage=='I',1,ifelse(lxdata$Stage=='II',2 ,ifelse(lxdata$Stage=='III',3,ifelse(lxdata$Stage=='IV',4,0))))
lxdata1$Stage <- factor(lxdata1$Stage,levels=c(0,1,2,3,4),labels=c('UNKOWN','I','II','III','IV'))
lxdata1$Risk=predict(mulcox,lxdata1)
lxdata1$CRRS=predict(mulcox,lxdata1)
mulcox1<- coxph(Surv(OS, OS.E) ~ CRRS+Stage, data = lxdata1 );mulcox1
lxdata1$CRRS_Stage=predict(mulcox1,lxdata1)


lxdata1$COL1A2_model=predict(univ_models$COL1A2,lxdata1)
lxdata1$F2_model=predict(univ_models$F2,lxdata1)
lxdata1$P2RX1_model=predict(univ_models$P2RX1,lxdata1)
lxdata1$PLAUR_model=predict(univ_models$PLAUR,lxdata1)
lxdata1$C4BPA_model=predict(univ_models$C4BPA,lxdata1)


troc.1= survivalROC(Stime=lxdata1$OS,  
                    status=lxdata1$OS.E,   
                    marker = lxdata1$Risk,     
                    predict.time = 12, method="KM")
troc.13= survivalROC(Stime=lxdata1$OS,  
                     status=lxdata1$OS.E,   
                     marker = lxdata1$Risk,     
                     predict.time = 36, method="KM")
troc.15= survivalROC(Stime=lxdata1$OS,  
                     status=lxdata1$OS.E,   
                     marker = lxdata1$Risk,     
                     predict.time = 60, method="KM")
troc.2= survivalROC(Stime=lxdata1$OS,  
                    status=lxdata1$OS.E,   
                    marker = lxdata1$F2_model,     
                    predict.time =12, method="KM")
troc.3= survivalROC(Stime=lxdata1$OS,  
                    status=lxdata1$OS.E,   
                    marker = lxdata1$P2RX1_model,     
                    predict.time = 12, method="KM")

troc.5= survivalROC(Stime=lxdata1$OS,  
                    status=lxdata1$OS.E,   
                    marker = lxdata1$PLAUR_model,     
                    predict.time = 12, method="KM")

troc.10= survivalROC(Stime=lxdata1$OS,
                     status=lxdata1$OS.E,
                     marker = lxdata1$C4BPA_model,
                     predict.time = 12, method="KM")


library(xlsx)
write.xlsx(lxdata1[,c(1,31,32,35:ncol(lxdata1))],file='files/Fig2.xlsx',sheetName = 'Fig2F',row.names = F,append = T)
png("figures/02F_rocmethods_bsr_TCGA.png",width = 1500,height = 1500,res=200)
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
                   paste('F12=',round(troc.7$AUC,3)),
                   paste('SERPINE1=',round(troc.8$AUC,3)),
                   paste('ITGB1=',round(troc.9$AUC,3)),
                   paste('C4BPA=',round(troc.10$AUC,3)),
                   paste('JMJD7-PLA2G4B=',round(troc.11$AUC,3)),
                   paste('CR2=',round(troc.12$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','blue',"green","purple",'orange','black',"yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)
dev.off()

# Fig2E
library(pec)
library(survival)
dat1 <- lxdata1[,c(1:12,31,32,35)]
write.xlsx(dat1,file = 'files/Fig2-1.xlsx',sheetName = 'Fig2E',row.names = F,append = T)
cox1 <- coxph(Surv(OS, OS.E) ~ Risk, data = dat1,x=T,y=T)
cox2 <- coxph(Surv(OS, OS.E) ~ F2, data = dat1,x=T,y=T)
cox3 <- coxph(Surv(OS, OS.E) ~ P2RX1, data = dat1,x=T,y=T)
cox4 <- coxph(Surv(OS, OS.E) ~ COL1A2, data = dat1,x=T,y=T)
cox5 <- coxph(Surv(OS, OS.E) ~ PLAUR, data = dat1,x=T,y=T)
cox10 <- coxph(Surv(OS, OS.E) ~ C4BPA, data = dat1,x=T,y=T)

ApparrentCindex  <- pec::cindex(list("Risk"=cox1,
                                     "F2"=cox2,
                                     "P2RX1"=cox3,
                                     "COL1A2"=cox4,
                                     "PLAUR"=cox5,
                                     "C4BPA"=cox10
                                     ),
                                formula=Surv(OS,OS.E)~Risk,
                                data=dat1,
                                eval.times=seq(1,242,20))

png("figures/02E_cindex.png",width = 1500,height = 1500,res=200)
plot(ApparrentCindex,lwd=4,cex.axis=2,legend.x=190,legend.y=1,legend.cex=1,col=c('red',"purple3",'blue',"green","purple",'orange','black',"yellow","pink","cyan","magenta","seagreen4","red4"))
dev.off()

# Fig3C
png("figures/03C_rocyear_risk_TCGA.png",width = 1500,height = 1500,res=200)
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
# Fig3D-E
dat2 <- lxdata1[,c(1,31:35)]
res.cut <- surv_cutpoint(dat2,
                         time='OS',
                         event = 'OS.E',
                         variables = 'Risk')
summary(res.cut)
dat2$group <- ifelse(dat2$Risk>0.192116,'High risk','Low risk')
write.xlsx(dat2[,c('sampleID','Risk','OS','OS.E','group')],file = 'figures/02_data.xlsx',sheetName = 'Sheet3',row.names = F,append = T)
png("figures/03D_OSKM_risk_TCGA.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(OS, OS.E) ~ group, data=dat2),
                surv.median.line = "hv",size=3,
                title = "TCGA-LUAD(n=502)",
                legend.labs = c('High-risk','Low-risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Overall survival',xlab='Time(months)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("figures/03E_RFSKM_risk_TCGA.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(RFS, RFS.E) ~ group, data=dat2),size=3,
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

###############################################
library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
# Fig3B
dat3 <- lxdata1[,c(1,21,35)]
dat3$group <- ifelse(dat3$Risk>median(dat3$Risk),'High risk','Low risk')
write.xlsx(dat3,file = 'files/Fig3.xlsx',sheetName = 'Fig3B',row.names = F,append = T)
#风险评分box图
cb <- ggboxplot(dat3,x = 'Status',y = 'Risk',fill='Status')+theme_bw()+ylab('Riskscore')+xlab('')+ggtitle('TCGA')+ylim(-2,3)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_signif(size=1,textsize=12,comparisons = list(c("Dead", "Alive")),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('figures/03B_riskbox.png',cb,width = 5,height = 5,dpi=200)

# Fig3A
#点图
dat4 <- lxdata1[order(lxdata1$Risk,decreasing = F),c(1:12,21,31,35)]
dat4$Sample <- 1:nrow(dat4)
dat4$group <- ifelse(dat4$Risk>median(dat4$Risk),'High risk','Low risk')
write.xlsx(dat4,file = 'files/Fig3.xlsx',sheetName = 'Fig3A',row.names = F,append = T)
sp <- ggplot(dat4,aes(Sample,OS,color=Status))+geom_point()+geom_vline(aes(xintercept = 251),lty=2,lwd=1)+
  theme_bw()+scale_colour_manual(values = c('blue','red'))+ylab('Survival month')+xlab('')+
  theme(axis.text.x = element_text(size=16,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
rp <- ggplot(dat4,aes(Sample,Risk,color=group))+geom_point()+geom_vline(aes(xintercept = 251),lty=2,lwd=1)+
  theme_bw()+scale_colour_manual(values = c('red','blue'))+ylab('Risk score')+xlab('')+ggtitle('TCGA')+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.justification=c(0,1),legend.position=c(0.01,0.99),
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('figures/03A_survpoint.png',sp,width = 5,height =3,dpi=1000)
ggsave('figures/03A_riskpoint.png',rp,width = 5,height = 3,dpi=1000)

library(pheatmap)
bsrgenes[10] <- gsub('_','-',bsrgenes[10])
gene <- bsrgenes
dat5 <- dat4[,colnames(dat4)%in%gene]
library(pheatmap)
group=dat4$group%>%as.data.frame()
colnames(group) <- 'group'
rownames(group) <- rownames(dat4)
group$group <- factor(group$group,levels=c('Low risk','High risk'))
data <- t(dat5)
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
save_pheatmap_pdf(p, "figures/03A_pheat.pdf")



