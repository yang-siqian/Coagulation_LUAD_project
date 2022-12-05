hubgenes[hubgenes%in%rownames(data)]
#c('F2','P2RX1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7-PLA2G4B','CR2')
#GUCY1A1 没有
########################################FUSCC
library(dplyr)
load('LUADdata/FUSCC/omics_clin.Rdata')
hubgenes <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7-PLA2G4B','CR2')
rna_data<-read.table('LUADdata/FUSCC/all.counts.id1.txt',row.names = 1,header = T)
rna_data<-rna_data[,c(6:409)]
names<-c("X1684LC","X1890N","X2310N","X3238LC", "X3244LC","X3246LC","X3255LC","X3284LC","X3338LC" 
         ,"X3439LC","X4873LC" ,"X4894LC", "X4912LC","X4933LC","X4941LC","X4942LC"  
         ,"X4948LC","X5069LC","X5081LC","X5091LC")
colnames(rna_data)<-c(sub("X","",colnames(rna_data)[colnames(rna_data)%in%names]),colnames(rna_data)[which(!colnames(rna_data)%in%names)])

load('LUADdata/FUSCC/clin_data.Rdata')

Barcode.T<-clin_data$RNA.LC
Barcode.N<-clin_data$RNA.N
Barcodes<-c(as.character(Barcode.T),as.character(Barcode.N))
Data.T<-rna_data[,which(colnames(rna_data)%in%Barcode.T)]
Data.N<-rna_data[,which(colnames(rna_data)%in%Barcode.N)]
rna_data1<-rna_data[,colnames(rna_data)%in%Barcodes]
data<-rna_data1[,colnames(rna_data1)[which(colnames(rna_data1)%in%Barcodes)]]
group_list<-ifelse(colnames(data)%in%as.character(Barcode.T),"Tumor","Normal")
modeldat <- as.data.frame(t(data[rownames(data)%in%hubgenes,]))
modeldat<- do.call(data.frame,modeldat )
modeldat<- log2(modeldat+1)
rownames(modeldat) <- colnames(data)
modeldat$RNA.LC<- rownames(modeldat)
modeldat <- merge(modeldat,clin_data,by.y='RNA.LC')
modeldat <-  modeldat[!is.na(modeldat$RFS.E),]
modeldat$Status <- ifelse(modeldat$OS.E =='0','Alive','Dead')
modeldat$Status  <-factor(modeldat$Status ,levels=c("Alive", "Dead"))
modeldat$Risk <- 0.09*modeldat$F2-0.299*modeldat$P2RX1+0.107*modeldat$PLAUR-0.124*modeldat$PRKCZ+0.093*modeldat$F12+0.056*modeldat$SERPINE1+
  0.137*modeldat$COL1A2-0.022*modeldat$C4BPA-0.031*modeldat$JMJD7.PLA2G4B-0.017*modeldat$CR2
lxdata <- modeldat[,c('age','Pathologic_stage','Risk','OS','OS.E')]
lxdata1 <- lxdata
lxdata1$Age <- ifelse(lxdata$age<70,1,2 )
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$Stage <- ifelse(lxdata$Pathologic_stage%in%c('IA','IB'),1,3)
lxdata1$Stage <- factor(lxdata1$Stage,levels=c(1,3),labels=c('I','III'))

lxdata1$riskscore1 <- -0.3797 *as.numeric(lxdata1$Age)
lxdata1$riskscore2 <- 0.51381*as.numeric(lxdata1$Stage)
lxdata1$riskscore3 <- 0.42552  *as.numeric(lxdata1$Age)+0.36421 *as.numeric(lxdata1$Stage)+0.91900 *lxdata1$Risk

library(survivalROC)
library(survival)
library(survminer)
#12,36,60
time=36
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

png("CRG/result/figs/09_1-yearroc_combined.tiff",width = 1500,height = 1500,res=200)
png("CRG/result/figs/09_3-yearroc_combined.tiff",width = 1500,height = 1500,res=200)

plot(troc$FP, troc$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)

title(main='1-Year',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
title(main='3-Year',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
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

library(pec)
library(survival)
dat <- lxdata1[,c('OS', 'OS.E', 'Age','Stage','Risk')]
dat$time <- dat$OS *30
library(rms)
dd <- datadist(dat)
options(datadist='dd')

f2 <- psm(Surv(time, OS.E) ~ Age+Stage+Risk, data = dat,x=T,y=T,dist='lognormal' )
call <- calibrate(f2,cmethod = 'KM',method='boot',u=365,m=30,B=1000)
png("CRG/result/figs/09_cal_1y.tiff",width = 1500,height = 1500,res=200)
plot(call,lty=1,conf.int=T,errbar.col='blue',col='red',xlim=c(0,1),ylim=c(0,1), lwd=4,cex.axis=1.6,ann = F,subtitles=F)
title(xlab='Nomogram-predicted OS of 1-year',ylab = 'Actual 1-year OS',cex.lab=1.6,font.lab=2)
dev.off()

call2 <- calibrate(f2,cmethod = 'KM',method='boot',u=1095,m=30,B=1000)
png("CRG/result/figs/09_cal_3y.tiff",width = 1500,height = 1500,res=200)
plot(call2,lty=1,conf.int=T,errbar.col='blue',col='red',xlim=c(0,1),ylim=c(0,1), lwd=4,cex.axis=1.6,ann = F,subtitles=F)
title(xlab='Nomogram-predicted OS of 3-year',ylab = 'Actual 3-year OS',cex.lab=1.6,font.lab=2)
dev.off()


