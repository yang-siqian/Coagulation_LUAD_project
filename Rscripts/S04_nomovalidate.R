library(dplyr)
gene <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','ITGB1','C4BPA','JMJD7-PLA2G4B','CR2')
load('Rdata/GSE72094.Rdata')
########################################GSE72094
gene[!gene%in%rownames(exp398)]
exp <- exp398 %>% t %>% as.data.frame()
exp <- exp[,colnames(exp)%in%gene]
exp$sample <- rownames(exp)
clin <- merge(clindata398,exp,by.x='sample')
clin$Status <- ifelse(clin$status=='Alive',0,1)
clin$Time <- round(as.numeric(clin$time)/30,2)
  
clin$Risk=0.083*clin$F2-0.256*clin$P2RX+0.127*clin$PLAUR-0.163*clin$PRKCZ+0.098*clin$F12+0.083*clin$SERPINE1+
  0.076*clin$ITGB1-0.029*clin$C4BPA-0.021*clin$CR2
lxdata <- clin[,c('sample','age','stage','time','Time','Status','Risk')]
lxdata1 <- lxdata
lxdata1$Age <- ifelse(lxdata$age<70,1,2)
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$stage<- ifelse(lxdata$stage%in%c('1','1A','1B'),1, ifelse(lxdata$stage%in%c('2A','2B'),2, ifelse(lxdata$stage%in%c('3A','3B'),3, ifelse(lxdata$stage=='4',4,0))))
lxdata1$stage <- factor(lxdata1$stage,levels=c(0,1,2,3,4),labels=c('UNKOWN','I','II','III','IV'))
lxdata1$riskscore1 <- 0.3797 *as.numeric(lxdata1$age)
lxdata1$riskscore2 <- 0.51381*as.numeric(lxdata1$stage)
lxdata1$riskscore3 <- 0.42282  *as.numeric(lxdata1$age)+0.37025 *as.numeric(lxdata1$stage)+0.90638 *lxdata1$Risk

# SupFig4B
library(survivalROC)
#12,36,60
for(i in c(12,36,60)){
  time=i
  troc= survivalROC(Stime=lxdata1$Time,  
                    status=lxdata1$Status,   
                    marker = lxdata1$riskscore3,     
                    predict.time = time, method="KM")
  troc.1= survivalROC(Stime=lxdata1$Time,  
                      status=lxdata1$Status,   
                      marker = lxdata1$Risk,     
                      predict.time = time, method="KM")
  troc.2= survivalROC(Stime=lxdata1$Time,  
                      status=lxdata1$Status,   
                      marker = lxdata1$riskscore1,     
                      predict.time = time, method="KM")
  troc.3= survivalROC(Stime=lxdata1$Time,  
                      status=lxdata1$Status,   
                      marker = lxdata1$riskscore2,     
                      predict.time = time, method="KM")
  
  png(paste0("figures/S04B_",i/12,"-yearroc_combined.png"),width = 1500,height = 1500,res=200)
  plot(troc$FP, troc$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
  
  title(main=paste0(i/12,'-Year'),xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
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
}





# SupFig4A
library(pec)
library(survival)
dat <- lxdata1[,c('time', 'Status', 'Age','stage','Risk')]
dat$time <- as.numeric(dat$time)
library(rms)
dd <- datadist(dat)
options(datadist='dd')
f2 <- psm(Surv(time, Status) ~ Age+stage+Risk, data = dat,x=T,y=T,dist='lognormal' )
call <- calibrate(f2,cmethod = 'KM',method='boot',u=365,m=100,B=1000)

png("figures/S04A_cal_1y.png",width = 1500,height = 1500,res=200)
plot(call,lty=1,conf.int=T,errbar.col='blue',col='red',xlim=c(0,1),ylim=c(0,1), lwd=4,cex.axis=1.6,ann = F,subtitles=F)
title(xlab='Nomogram-predicted OS of 1-year',ylab = 'Actual 1-year OS',cex.lab=1.6,font.lab=2)
dev.off()

call2 <- calibrate(f2,cmethod = 'KM',method='boot',u=1095,m=100,B=1000)
png("figures/S04A_cal_3y.png",width = 1500,height = 1500,res=200)
plot(call2,lty=1,conf.int=T,errbar.col='blue',col='red',xlim=c(0,1),ylim=c(0,1), lwd=4,cex.axis=1.6,ann = F,subtitles=F)
title(xlab='Nomogram-predicted OS of 3-year',ylab = 'Actual 3-year OS',cex.lab=1.6,font.lab=2)
dev.off()

call3 <- calibrate(f2,cmethod = 'KM',method='boot',u=1825,m=60,B=1000)
png("figures/S04A_cal_5y.png",width = 1500,height = 1500,res=200)
plot(call3,lty=1,conf.int=T,errbar.col='blue',col='red',xlim=c(0,1),ylim=c(0,1), lwd=4,cex.axis=1.6,ann = F,subtitles=F)
title(xlab='Nomogram-predicted OS of 5-year',ylab = 'Actual 5-year OS',cex.lab=1.6,font.lab=2)
dev.off()

write.xlsx(lxdata1,file='files/SupFig4.xlsx',sheetName = 'SupFig4_clin',row.names = F)
