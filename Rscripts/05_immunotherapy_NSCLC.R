#########GSE93157_lusc
library(survival)
library(survminer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(readxl)
CRRS <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','ITGB1','C4BPA','JMJD7-PLA2G4B','CR2')
clin1 <- read.table('rawdata/NSCLC/LUSC-GSE93157.clin.tsv',header=T,row.names = 1)
exp1 <- read.table('rawdata/NSCLC/LUSC-GSE93157.Response.tsv',header=T,row.names = 1)
clin1$OS <- clin1$overall_survival_.days./30
clin1$OS.E <- ifelse(clin1$vital_status=='Dead',1,0)
clin1$Response <- ifelse(clin1$response_NR%in%c('PD','SD'),'SD/PD','CR/PR')
exp1 <- exp1[exp1$GENE_SYMBOL%in%CRRS,]
exp1 <- na.omit(exp1)
rownames(exp1) <- exp1$GENE_SYMBOL
exp1 <- exp1[,-1] %>% t %>% as.data.frame()
exp1$patient_name <- rownames(exp1)
survdata <- merge(exp1,clin1,by.x='patient_name')
survdata$Risk <- 0.127*survdata$PLAUR+0.098*survdata$F12+0.076*survdata$ITGB1-0.029*survdata$C4BPA-0.021*survdata$CR2
survdata$group <- ifelse(survdata$Risk>median(survdata$Risk),'High risk','Low risk')

png("figures/05E_OSKM_risk_LUSC.png",width = 1500,height = 1500,res=200)
p1 <- ggsurvplot(survfit(Surv(OS, OS.E) ~ group, data=survdata),size=3,
                surv.median.line = "hv",
                title = "LUSC_PD1",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p1
dev.off()
png("figures/S07A_OSKM_risk_response_LUSC.png",width = 1500,height = 1500,res=200)
p1 <- ggsurvplot(survfit(Surv(OS, OS.E) ~ Response, data=survdata),size=3,
                surv.median.line = "hv",
                title = "LUSC_PD1",
                legend.labs = c('CR/PR','SD/PD'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p1
dev.off()

# Fig5C_LUSC
library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
library(reshape2)
a <-data.frame(group=c('High risk','Low risk'),CR_PR=c(t(table(survdata$Response[survdata$group=='High risk'])/length(survdata$Response[survdata$group=='High risk']))[1],
                                                       t(table(survdata$Response[survdata$group=='Low risk'])/length(survdata$Response[survdata$group=='Low risk']))[1]),
               SD_PD=c(t(table(survdata$Response[survdata$group=='High risk'])/length(survdata$Response[survdata$group=='High risk']))[2],
                       t(table(survdata$Response[survdata$group=='Low risk'])/length(survdata$Response[survdata$group=='Low risk']))[2]))
b=melt(a,id.vars='group')
b$variable <- gsub('_','/',b$variable)
# write.xlsx(b,file = 'files/Fig5.xlsx',sheetName = 'Fig5C_ACT',row.names = F,append = T)
P <- ggplot(data=b,aes(x=group,y=value,fill=variable)) + theme_bw()+ylab('Percent')+xlab('')+ggtitle('LUSC_PD1')+ylim(0,1)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_col(position = 'stack', width = 0.6)+
  geom_bar(position = "stack", stat = "identity", width = 0.6) +
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('figures/05E_barplot_LUSC.png',P,width = 5,height = 5,dpi=200)

# SupFig3B_LUSC
# write.xlsx(survdata1[,c(1,12:14)],file = 'files/SupFig3.xlsx',sheetName = 'SupFig3B_ACT',row.names = F,append = T)
library(survivalROC)
pd=survdata
troc.1= survivalROC(Stime=pd$OS,  
                    status=pd$OS.E,   
                    marker = pd$Risk,     
                    predict.time = 12, method="KM")
troc.13= survivalROC(Stime=pd$OS,  
                     status=pd$OS.E,   
                     marker = pd$Risk,       
                     predict.time = 36, method="KM")

png("figures/S07B_rocyear_risk_LUSC.png",width = 1500,height = 1500,res=200)
plot(troc.1$FP, troc.1$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
title(main='LUSC_PD1',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.13$FP, troc.13$TP, type="l",lwd=4,col="orange", xlim=c(0,1), ylim=c(0,1))
legend(0.55,0.2,c(paste('AUC of 1-year =',round(troc.1$AUC,3)),
                  paste('AUC of 3-year =',round(troc.13$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()

# SupFig3C_LUSC
# write.xlsx(survdata1[,c(1,12,17)],file = 'files/SupFig3.xlsx',sheetName = 'SupFig3C_ACT',row.names = F,append = T)
library(pROC)
roc1 <- roc(pd$Response,pd$Risk)
png("figures/S07C_roc_risk_LUSC.png",width = 1500,height = 1500,res=200)
plot(roc1,col="red",lwd=4,cex.axis=1.5,ann = F, xlim=c(1,0), ylim=c(0,1))
title(xlab="Specificity", ylab="Sensitivity",cex.lab=1.6,font.lab=2)
legend(0.55,0.2,c(paste('AUC =',round(roc1$auc,3)) ),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()





########GSE93157_LUAD
exp2 <- read.table('rawdata/NSCLC/nonsqNSCLC-GSE93157.Response.tsv',header=T,row.names = 1)
clin2 <- read.table('rawdata/NSCLC/nonsqNSCLC-GSE93157.clin.tsv',header=T,row.names = 1)
clin2$OS <- clin2$overall_survival_.days./30
clin2$OS.E <- ifelse(clin2$vital_status=='Dead',1,0)
clin2$Response <- ifelse(clin2$response_NR%in%c('PD','SD'),'SD/PD','CR/PR')
exp2 <- exp2[exp2$GENE_SYMBOL%in%CRRS,]
exp2 <- na.omit(exp2)
rownames(exp2) <- exp2$GENE_SYMBOL
exp2 <- exp2[,-1] %>% t %>% as.data.frame()
exp2$patient_name <- rownames(exp2)
survdata1 <- merge(exp2,clin2,by.x='patient_name')
survdata1$Risk <- 0.127*survdata1$PLAUR+0.098*survdata1$F12+0.076*survdata1$ITGB1-0.029*survdata1$C4BPA-0.021*survdata1$CR2
survdata1$group <- ifelse(survdata1$Risk>median(survdata1$Risk),'High risk','Low risk')
png("figures/05E_OSKM_risk_LUAD.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(OS, OS.E) ~ group, data=survdata1),size=3,
                surv.median.line = "hv",
                title = "LUAD_PD1",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("figures/S07A_OSKM_risk_response_LUAD.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(OS, OS.E) ~ Response, data=survdata1),size=3,
                surv.median.line = "hv",
                title = "LUAD_PD1",
                legend.labs = c('CR/PR','SD/PD'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

# Fig5C_LUAD
library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
library(reshape2)
a <-data.frame(group=c('High risk','Low risk'),CR_PR=c(t(table(survdata1$Response[survdata1$group=='High risk'])/length(survdata1$Response[survdata1$group=='High risk']))[1],
                                                       t(table(survdata1$Response[survdata1$group=='Low risk'])/length(survdata1$Response[survdata1$group=='Low risk']))[1]),
               SD_PD=c(t(table(survdata1$Response[survdata1$group=='High risk'])/length(survdata1$Response[survdata1$group=='High risk']))[2],
                       t(table(survdata1$Response[survdata1$group=='Low risk'])/length(survdata1$Response[survdata1$group=='Low risk']))[2]))
b=melt(a,id.vars='group')
b$variable <- gsub('_','/',b$variable)
# write.xlsx(b,file = 'files/Fig5.xlsx',sheetName = 'Fig5C_ACT',row.names = F,append = T)
P1 <- ggplot(data=b,aes(x=group,y=value,fill=variable)) + theme_bw()+ylab('Percent')+xlab('')+ggtitle('LUAD_PD1')+ylim(0,1)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_col(position = 'stack', width = 0.6)+
  geom_bar(position = "stack", stat = "identity", width = 0.6) +
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('figures/05E_barplot_LUAD.png',P1,width = 5,height = 5,dpi=200)

# SupFig3B_LUSC
# write.xlsx(survdata1[,c(1,12:14)],file = 'files/SupFig3.xlsx',sheetName = 'SupFig3B_ACT',row.names = F,append = T)
library(survivalROC)
pd=survdata1
troc.1= survivalROC(Stime=pd$OS,  
                    status=pd$OS.E,   
                    marker = pd$Risk,     
                    predict.time = 12, method="KM")
troc.13= survivalROC(Stime=pd$OS,  
                     status=pd$OS.E,   
                     marker = pd$Risk,       
                     predict.time = 36, method="KM")


png("figures/S07B_rocyear_risk_LUAD.png",width = 1500,height = 1500,res=200)
plot(troc.1$FP, troc.1$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
title(main='LUAD_PD1',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.13$FP, troc.13$TP, type="l",lwd=4,col="orange", xlim=c(0,1), ylim=c(0,1))
legend(0.55,0.2,c(paste('AUC of 1-year =',round(troc.1$AUC,3)),
                  paste('AUC of 3-year =',round(troc.13$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()

# SupFig3C_LUAD
# write.xlsx(survdata1[,c(1,12,17)],file = 'files/SupFig3.xlsx',sheetName = 'SupFig3C_ACT',row.names = F,append = T)
library(pROC)
roc1 <- roc(pd$Response,pd$Risk)
png("figures/S07C_roc_risk_LUAD.png",width = 1500,height = 1500,res=200)
plot(roc1,col="red",lwd=4,cex.axis=1.5,ann = F, xlim=c(1,0), ylim=c(0,1))
title(xlab="Specificity", ylab="Sensitivity",cex.lab=1.6,font.lab=2)
legend(0.55,0.2,c(paste('AUC =',round(roc1$auc,3)) ),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()

##########################################合并
survdata2 <- rbind(survdata,survdata1)
write.xlsx(survdata2[,c(1,22:26)],file = 'files/Fig5.xlsx',sheetName = 'Fig5C_NSCLC_surv',row.names = F,append = T)
png("figures/05E_OSKM_risk_NSCLC.png",width = 1500,height = 1500,res=200)
ggsurvplot(survfit(Surv(OS, OS.E) ~ group, data=survdata2),size=3,
           surv.median.line = "hv",
           title = "NSCLC_PD1",
           legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
           palette = c("#E69F00", "#56B4E9",'grey'),
           pval = T,pval.size=18,legend.title = "",
           ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
           font.main = c(30, "bold"),
           ylab='Survival probability of OS',xlab='Time(month)',
           font.x = c(30, "bold"),
           font.y = c(30, "bold"))
dev.off()

png("figures/S07A_OSKM_risk_response_NSCLC.png",width = 1500,height = 1500,res=200)
ggsurvplot(survfit(Surv(OS, OS.E) ~ Response, data=survdata2),size=3,
           surv.median.line = "hv",
           title = "NSCLC_PD1",
           legend.labs = c('CR/PR','SD/PD'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
           palette = c("#E69F00", "#56B4E9",'grey'),
           pval = T,pval.size=18,legend.title = "",
           ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
           font.main = c(30, "bold"),
           ylab='Survival probability of OS',xlab='Time(month)',
           font.x = c(30, "bold"),
           font.y = c(30, "bold"))
dev.off()

a <-data.frame(group=c('High risk','Low risk'),CR_PR=c(t(table(survdata2$Response[survdata2$group=='High risk'])/length(survdata2$Response[survdata2$group=='High risk']))[1],
                                                       t(table(survdata2$Response[survdata2$group=='Low risk'])/length(survdata2$Response[survdata2$group=='Low risk']))[1]),
               SD_PD=c(t(table(survdata2$Response[survdata2$group=='High risk'])/length(survdata2$Response[survdata2$group=='High risk']))[2],
                       t(table(survdata2$Response[survdata2$group=='Low risk'])/length(survdata2$Response[survdata2$group=='Low risk']))[2]))
b=melt(a,id.vars='group')
b$variable <- gsub('_','/',b$variable)
write.xlsx(b,file = 'files/Fig5.xlsx',sheetName = 'Fig5C_NSCLC',row.names = F,append = T)
P1 <- ggplot(data=b,aes(x=group,y=value,fill=variable)) + theme_bw()+ylab('Percent')+xlab('')+ggtitle('NSCLC_PD1')+ylim(0,1)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_col(position = 'stack', width = 0.6)+
  geom_bar(position = "stack", stat = "identity", width = 0.6) +
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('figures/05E_barplot_NSCLC.png',P1,width = 5,height = 5,dpi=200)

# SupFig3B_LUSC
write.xlsx(survdata1[,c(1,22,23,25)],file = 'files/SupFig3.xlsx',sheetName = 'SupFig3B_NSCLC',row.names = F,append = T)
library(survivalROC)
pd=survdata2
troc.1= survivalROC(Stime=pd$OS,  
                    status=pd$OS.E,   
                    marker = pd$Risk,     
                    predict.time = 12, method="KM")
troc.13= survivalROC(Stime=pd$OS,  
                     status=pd$OS.E,   
                     marker = pd$Risk,       
                     predict.time = 36, method="KM")


png("figures/S07B_rocyear_risk_NSCLC.png",width = 1500,height = 1500,res=200)
plot(troc.1$FP, troc.1$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
title(main='NSCLC_PD1',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.13$FP, troc.13$TP, type="l",lwd=4,col="orange", xlim=c(0,1), ylim=c(0,1))
legend(0.55,0.2,c(paste('AUC of 1-year =',round(troc.1$AUC,3)),
                  paste('AUC of 3-year =',round(troc.13$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()

# SupFig3C_NSCLC
library(xlsx)
write.xlsx(survdata2[,c(1,24,25)],file = 'files/SupFig3.xlsx',sheetName = 'SupFig3C_NSCLC',row.names = F,append = T)
library(pROC)
roc1 <- roc(pd$Response,pd$Risk)
png("figures/S07C_roc_risk_NSCLC.png",width = 1500,height = 1500,res=200)
plot(roc1,col="red",lwd=4,cex.axis=1.5,ann = F, xlim=c(1,0), ylim=c(0,1))
title(xlab="Specificity", ylab="Sensitivity",cex.lab=1.6,font.lab=2)
legend(0.55,0.2,c(paste('AUC =',round(roc1$auc,3)) ),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()


