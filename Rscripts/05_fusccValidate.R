
########################################FUSCC
library(dplyr)
load('rawdata/FUSCC/omics_clin.Rdata')
hubgenes <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7-PLA2G4B','CR2')
rna_data<-read.table('rawdata/FUSCC/all.counts.id1.txt',row.names = 1,header = T)
rna_data<-rna_data[,c(6:409)]
names<-c("X1684LC","X1890N","X2310N","X3238LC", "X3244LC","X3246LC","X3255LC","X3284LC","X3338LC" 
         ,"X3439LC","X4873LC" ,"X4894LC", "X4912LC","X4933LC","X4941LC","X4942LC"  
         ,"X4948LC","X5069LC","X5081LC","X5091LC")
colnames(rna_data)<-c(sub("X","",colnames(rna_data)[colnames(rna_data)%in%names]),colnames(rna_data)[which(!colnames(rna_data)%in%names)])

load('rawdata/FUSCC/clin_data.Rdata')

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
modeldat$CRRS <-predict(mulcox,modeldat)
modeldat$Age <- ifelse(modeldat$age>70,'old','young')
modeldat$Gender<- ifelse(modeldat$sex%in%c("Female  ","Female","Female   "),'F','M')
modeldat$Stage <- ifelse(modeldat$Pathologic_stage%in%c("IA","IB"),'I','III')
modeldat$T<- ifelse(modeldat$pT%in%c("T1a","T1b"),'T1','T2-4')
modeldat$N<- modeldat$pN
mulcox1<- coxph(Surv(OS, OS.E) ~ CRRS+Age+Gender+T+N+TP53+EGFR, data = modeldat)
forest1 <- ggforest(mulcox1,  
                    data = modeldat,  
                    main = 'FUSCC overall survival', 
                    cpositions = c(0.05, 0.15, 0.35), 
                    fontsize = 2, 
                    refLabel = 'Reference', 
                    noDigits = 3 
)
mulcox2<- coxph(Surv(RFS,RFS.E) ~ CRRS+Age+Gender+T+N+TP53+EGFR, data = modeldat)
forest2 <- ggforest(mulcox2,  
                    data = modeldat,  
                    main = 'FUSCC recurrence-free', 
                    cpositions = c(0.05, 0.15, 0.35), 
                    fontsize = 2, 
                    refLabel = 'Reference', 
                    noDigits = 3 
)
par(mfrow=c(2,1))
p=forest1+forest2
p
ggsave(filename="figures/03S_fuscc_muti_forest.png",p, width=30, height=14)



library(survival)
library(survminer)
res.cut <- surv_cutpoint(modeldat,
                         time='OS',
                         event = 'OS.E',
                         variables = 'CRRS')
modeldat$group<- ifelse(modeldat$CRRS>median(modeldat$CRRS),'High risk','Low risk')

res.cut <- surv_cutpoint(modeldat,
                         time='OS',
                         event = 'OS.E',
                         variables = 'P2RX1')
modeldat$group <- ifelse(modeldat$P2RX1>as.numeric(summary(res.cut)[1]),'High P2RX1','Low P2RX1')
png("figures/07G_OSKM_P2RX1_FUSCC.png",width = 1500,height = 1500,res=200)


PP <- ggboxplot(modeldat,x = 'group',y = "P2RX1",fill='group')+theme_classic()+ylab('P2RX1 Expression')+xlab('')+ggtitle('FUSCC')+ylim(min(modeldat$P2RX1)-0.2,max(modeldat$P2RX1)+2)+
  scale_fill_manual(values=c("#FFD700", "#BF3EFF"))+
  geom_signif(size=1,textsize=12,comparisons = list(c('Tumor','Normal')),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1))


ggsave('figures/07D_P2RX1_FUSCC.png',PP,width = 5,height =5,dpi=200)


p <- ggsurvplot(survfit(Surv(OS, OS.E) ~ group, data=modeldat),
                surv.median.line = "hv",size=3,
                title = "FUSCC",
                legend.labs = c('High P2RX1','Low P2RX1'),legend=c(0.8,0.87),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Overall survival',xlab='Time(months)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("CRG/result/figs/04_OSKM_risk_FUSCC.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(OS, OS.E) ~ group, data=modeldat),
                surv.median.line = "hv",size=3,
                title = "FUSCC",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.5),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("CRG/result/figs/04_RFSKM_risk_FUSCC.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(RFS, RFS.E) ~ group, data=modeldat),size=3,
                surv.median.line = "hv",
                title = "FUSCC",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("red", "blue",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of RFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()
############################
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
#风险评分box图
cb <- ggboxplot(modeldat,x = 'Status',y = 'Risk',fill='Status')+theme_bw()+ylab('Riskscore')+xlab('')+ggtitle('FUSCC')+ylim(-1.5,2)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_signif(size=1,textsize=12,comparisons = list(c("Dead", "Alive")),test = wilcox.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('CRG/result/figs/04_riskbox.tiff',cb,width = 5,height = 5,dpi=200)
#点图
pd <- modeldat[order(modeldat$Risk,decreasing = F),]
colnames(pd)[9] <- gsub('\\.','-',colnames(pd)[9])
pd$Sample <- 1:nrow(pd)
sp <- ggplot(pd,aes(Sample,OS,color=Status))+geom_point()+geom_vline(aes(xintercept = 50),lty=2,lwd=1)+
  theme_bw()+scale_colour_manual(values = c('blue','red'))+ylab('Survival month')+xlab('')+
  theme(axis.text.x = element_text(size=16,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
rp <- ggplot(pd,aes(Sample,Risk,color=group))+geom_point()+geom_vline(aes(xintercept = 50),lty=2,lwd=1)+
  theme_bw()+scale_colour_manual(values = c('red','blue'))+ylab('Riskscore')+xlab('')+ggtitle('FUSCC')+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.justification=c(0,1),legend.position=c(0.01,0.99),
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('CRG/result/figs/04_survpoint.tiff',sp,width = 5,height =3,dpi=1000)
ggsave('CRG/result/figs/04_riskpoint.tiff',rp,width = 5,height = 3,dpi=1000)


dat <- pd[,2:11]
rownames(dat) <- pd$RNA.LC
library(pheatmap)
group=pd$group%>%as.data.frame()
colnames(group) <- 'group'
rownames(group) <- pd$RNA.LC
group$group <- factor(group$group,levels=c('Low risk','High risk'))
data <- t(dat)
colnames(data) <- pd$RNA.LC
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

save_pheatmap_pdf(p, "CRG/result/figs/04_pheat_fuscc.pdf")

save(pd,file = 'CRG/extdata/fuscc_surv.Rdata')
load('CRG/extdata/fuscc_surv.Rdata')
############################################ROC
library(survivalROC)
troc.1= survivalROC(Stime=pd$OS,  
                  status=pd$OS.E,   
                  marker = pd$Risk,     
                  predict.time = 12, method="KM")
troc.12= survivalROC(Stime=pd$OS,  
                     status=pd$OS.E,   
                     marker = pd$Risk,     
                     predict.time = 24, method="KM")
troc.13= survivalROC(Stime=pd$OS,  
                     status=pd$OS.E,   
                     marker = pd$Risk,     
                    predict.time = 36, method="KM")

png("CRG/result/figs/04_rocyear_risk_fuscc.tiff",width = 1500,height = 1500,res=200)
plot(troc.1$FP, troc.1$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
title(main='FUSCC',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.12$FP, troc.12$TP, type="l",col="orange", xlim=c(0,1), ylim=c(0,1),lwd=4)
lines(troc.13$FP, troc.13$TP, type="l",col="blue", xlim=c(0,1), ylim=c(0,1),lwd=4)
legend(0.55,0.2,c(paste('AUC of 1-year =',round(troc.1$AUC,3)),paste('AUC of 2-year =',round(troc.12$AUC,3)),paste('AUC of 3-year =',round(troc.13$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd = 4,col = c('red','orange','blue'),bty='n',seg.len =1.5,cex = 1.2)
dev.off()

##################################################
pd$F2_model=0.108 *pd$F2
pd$P2RX1_model=-0.213 *pd$P2RX1
pd$PLAUR_model=0.145 *pd$PLAUR
pd$PRKCZ_model=-0.194 *pd$PRKCZ 
pd$F12_model=0.133 *pd$F12
pd$SERPINE1_model=0.105 *pd$SERPINE1
pd$COL1A2_model=0.083 *pd$COL1A2
pd$C4BPA_model=-0.084 *pd$C4BPA
pd$JMJD7_PLA2G4B_model=-0.144 *pd$`JMJD7-PLA2G4B`
pd$CR2_model=-0.109 *pd$CR2

troc.1= survivalROC(Stime=pd$OS,  
                    status=pd$OS.E,   
                    marker = pd$Risk,     
                    predict.time = 12, method="KM")
troc.2= survivalROC(Stime=pd$OS,  
                    status=pd$OS.E,   
                    marker = pd$F2_model,     
                    predict.time =12, method="KM")
troc.3= survivalROC(Stime=pd$OS,  
                    status=pd$OS.E,   
                    marker = pd$P2RX1_model,     
                    predict.time = 12, method="KM")
troc.5= survivalROC(Stime=pd$OS,  
                    status=pd$OS.E,   
                    marker = pd$PLAUR_model,     
                    predict.time = 12, method="KM")
troc.6= survivalROC(Stime=pd$OS,  
                    status=pd$OS.E,   
                    marker = pd$PRKCZ_model,     
                    predict.time = 12, method="KM")
troc.7= survivalROC(Stime=pd$OS,  
                    status=pd$OS.E,   
                    marker = pd$F12_model,     
                    predict.time = 12, method="KM")
troc.8= survivalROC(Stime=pd$OS,  
                    status=pd$OS.E,   
                    marker = pd$SERPINE1_model,     
                    predict.time = 12, method="KM")
troc.9= survivalROC(Stime=pd$OS,  
                    status=pd$OS.E,   
                    marker = pd$COL1A2_model,     
                    predict.time = 12, method="KM")
troc.10= survivalROC(Stime=pd$OS,  
                     status=pd$OS.E,   
                     marker = pd$C4BPA_model,     
                     predict.time = 12, method="KM")
troc.11= survivalROC(Stime=pd$OS,  
                     status=pd$OS.E,   
                     marker = pd$JMJD7_PLA2G4B_model,     
                     predict.time = 12, method="KM")
troc.12= survivalROC(Stime=pd$OS,  
                     status=pd$OS.E,   
                     marker = pd$CR2_model,     
                     predict.time = 12, method="KM")



png("CRG/result/figs/03_rocmethods_bsr_FUSCC.tiff",width = 1500,height = 1500,res=200)
plot(troc.1$FP, troc.1$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
title(main='FUSCC',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.2$FP, troc.2$TP, type="l",lwd=4,col="blue", xlim=c(0,1), ylim=c(0,1))
lines(troc.3$FP, troc.3$TP, type="l",lwd=4,col="green", xlim=c(0,1), ylim=c(0,1))
lines(troc.5$FP, troc.5$TP, type="l",lwd=4,col="orange", xlim=c(0,1), ylim=c(0,1))
lines(troc.6$FP, troc.6$TP, type="l",lwd=4,col="black", xlim=c(0,1), ylim=c(0,1))
lines(troc.7$FP, troc.7$TP, type="l",lwd=4,col="yellow", xlim=c(0,1), ylim=c(0,1))
lines(troc.8$FP, troc.8$TP, type="l",lwd=4,col="pink", xlim=c(0,1), ylim=c(0,1))
lines(troc.9$FP, troc.9$TP, type="l",lwd=4,col="cyan", xlim=c(0,1), ylim=c(0,1))
lines(troc.10$FP, troc.10$TP, type="l",lwd=4,col="magenta", xlim=c(0,1), ylim=c(0,1))
lines(troc.11$FP, troc.11$TP, type="l",lwd=4,col="seagreen4", xlim=c(0,1), ylim=c(0,1))
lines(troc.12$FP, troc.12$TP, type="l",lwd=4,col="red4", xlim=c(0,1), ylim=c(0,1))
legend(0.55,0.5,c(paste('Risk=',round(troc.1$AUC,3)),
                   paste('F2=',round(troc.2$AUC,3)),
                   paste('P2RX1=',round(troc.3$AUC,3)),
                   paste('PLAUR=',round(troc.5$AUC,3)),
                   paste('PRKCZ=',round(troc.6$AUC,3)),
                   paste('LF12=',round(troc.7$AUC,3)),
                   paste('SERPINE1=',round(troc.8$AUC,3)),
                   paste('COL1A2=',round(troc.9$AUC,3)),
                   paste('C4BPA=',round(troc.10$AUC,3)),
                   paste('JMJD7-PLA2G4B=',round(troc.11$AUC,3)),
                   paste('CR2=',round(troc.12$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','blue',"green",'orange','black',"yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()


