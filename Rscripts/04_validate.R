library(dplyr)
gene <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7-PLA2G4B','CR2')
          
load('LUADdata/GEO/GSE13213.Rdata')
load('LUADdata/GEO/GSE72094.Rdata')
load('LUADdata/GEO/GSE68465.Rdata')
load('LUADdata/GEO/GSE42127.Rdata')
load('LUADdata/GEO/GSE37745.Rdata')
load('LUADdata/GEO/GSE31210.Rdata')

table(gene%in%rownames(exp226))
# "F2"    "P2RX1" "PRKCZ" "C4BPA"
gene[gene%in%rownames(exp106)]
# F2"    "P2RX1" "PRKCZ" "C4BPA"
gene[gene%in%rownames(exp117)] #无
gene[gene%in%rownames(exp133)]
# "P2RX1"    "F12"      "SERPINE1" "C4BPA"  
gene[gene%in%rownames(exp442)]
# "F2"    "P2RX1" "PRKCZ" "F12"   "C4BPA" "CR2"  
gene[gene%in%rownames(exp443)]
# "F2"    "P2RX1" "PRKCZ" "C4BPA" "CR2"
gene[gene%in%rownames(omics[[1]])]
# "F2"       "P2RX1"    "PLAUR"    "PRKCZ"    "F12"      "SERPINE1" "COL1A2"   "C4BPA"    "CR2"  

########################################GSE72094
exp <- exp442 %>% t %>% as.data.frame()
exp <- exp[,colnames(exp)%in%gene]
exp$sample <- rownames(exp)
clin <- merge(clindata442,exp,by.x='sample')
clin$Status <- ifelse(clin$status=='NA',2,ifelse(clin$status=='Alive',0,1))
clin <- clin[clin$Status!=2,]
clin <- clin[clin$time!='NA',]
clin$Time <- round(as.numeric(clin$time)/30,2)
clin$riskscore <- 0.08691*clin$F2 -0.31064*clin$P2RX1 + -0.13589*clin$PRKCZ+ 0.09316*clin$F12 -0.02263*clin$C4BPA -0.01604*clin$CR2  
clin$Group <- ifelse(clin$riskscore>median(clin$riskscore),'High risk','Low risk')
library(survival)
library(survminer)
ggsurvplot(survfit(Surv(Time, Status) ~ Group, data=clin),surv.median.line = "hv",palette = c('red','blue','grey'),pval = T,pval.size=10)

library(ggstatsplot)
library(tidyverse)
library(ggpubr)
#风险评分box图
cb <- ggboxplot(clin,x = 'status',y = 'riskscore',fill='status')+theme_classic()+ylab('Riskscore')+xlab('')+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_signif(textsize=5,comparisons = list(c("Dead", "Alive")),test = wilcox.test,map_signif_level = F)+
  theme(axis.text.x = element_text(size=16,face='bold',color='black'),
        axis.text.y = element_text(size=16,face='bold',color='black'),
        text=element_text(size=16,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1))
ggsave('CRG/result/figs/10_riskbox442.tiff',cb,width = 5,height = 5,dpi=1000)

#点图
pd <- clin[order(clin$riskscore,decreasing = F),]
pd$Sample <- 1:nrow(pd)
sp <- ggplot(pd,aes(Sample,Time,color=status))+geom_point()+geom_vline(aes(xintercept = 199),lty=2,lwd=1)+
  theme_bw()+scale_colour_manual(values = c('blue','red'))+ylab('Survival month')+xlab('')+
  theme(axis.text.x = element_text(size=16,face='bold',color='black'),
        axis.text.y = element_text(size=16,face='bold',color='black'),
        text=element_text(size=16,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
rp <- ggplot(pd,aes(Sample,riskscore,color=Group))+geom_point()+geom_vline(aes(xintercept = 199),lty=2,lwd=1)+
  theme_bw()+scale_colour_manual(values = c('red','blue'))+ylab('Risk score')+xlab('')+
  theme(axis.text.x = element_text(size=16,face='bold',color='black'),
        axis.text.y = element_text(size=16,face='bold',color='black'),
        text=element_text(size=16,face='bold'),
        legend.title = element_blank(),legend.justification=c(0,1),legend.position=c(0.01,0.99),
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('CRG/result/figs/10_survpoint442.tiff',sp,width = 5,height =3,dpi=1000)
ggsave('CRG/result/figs/10_riskpoint442.tiff',rp,width = 5,height = 3,dpi=1000)

library(pheatmap)
dat <- pd[,8:13]
rownames(dat) <- pd$sample
library(pheatmap)
group=pd$Group%>%as.data.frame()
colnames(group) <- 'group'
rownames(group) <- pd$sample
group$group <- factor(group$group,levels=c('Low risk','High risk'))
data <- t(dat)
colnames(data) <- pd$sample
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

save_pheatmap_pdf(p, "CRG/result/figs/10_pheat442.pdf")

library(writexl)
writexl::write_xlsx(pd[,c('Time','Status','riskscore')],'CRG/result/files/01_survroc442.xlsx')
save(pd,file = 'CRG/extdata/GSE72094_surv442.Rdata')
load( 'CRG/extdata/GSE72094_surv442.Rdata')
library(survivalROC)
troc= survivalROC(Stime=pd$Time,  
                  status=pd$Status,   
                  marker = pd$riskscore,     
                  predict.time = 12, method="KM")
troc.3= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$riskscore,     
                    predict.time = 36, method="KM")
troc.5= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$riskscore,     
                    predict.time = 60, method="KM")
png("CRG/result/figs/10_roc442.tiff",width = 1500,height = 1500,res=200)
plot(troc$FP, troc$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=3,cex.axis=1.5,ann = F)
title(xlab="1-Specificity", ylab="Sensitivity",cex.lab=1.5,font.lab=4)
abline(0,1,col='gray',lty=2,lwd=2)
lines(troc.3$FP, troc.3$TP, type="l",col="blue", xlim=c(0,1), ylim=c(0,1),lwd=2)
lines(troc.5$FP, troc.5$TP, type="l",col="green", xlim=c(0,1), ylim=c(0,1),lwd=2)
legend(0.45,0.3,c(paste('AUC of 1-year =',round(troc$AUC,3)),paste('AUC of 3-year =',round(troc.3$AUC,3)),paste('AUC of 5-year =',round(troc.5$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd = 2,col = c('red','blue',"green"),bty='n',seg.len = 2,cex = 1.3)
dev.off()
# "F2"    "P2RX1" "PRKCZ" "F12"   "C4BPA" "CR2"  
pd$riskscore1 <- 0.10843*as.numeric(pd$F2)
pd$riskscore2 <- -0.21305*as.numeric(pd$P2RX1)
pd$riskscore5 <- -0.19426*as.numeric(pd$PRKCZ)
pd$riskscore9 <- 0.8969*as.numeric(pd$CR2)
pd$riskscore10 <- 0.13263 *as.numeric(pd$F12)
library(survivalROC)
troc= survivalROC(Stime=pd$Time,  
                  status=pd$Status,   
                  marker = pd$riskscore,     
                  predict.time = 12, method="KM")
troc.1= survivalROC(Stime=pd$Time,  
                    status=pd$Status,    
                    marker = pd$riskscore1,     
                    predict.time = 12, method="KM")
troc.2= survivalROC(Stime=pd$Time,  
                    status=pd$Status,    
                    marker = pd$riskscore2,     
                    predict.time = 12, method="KM")
troc.5= survivalROC(Stime=pd$Time,  
                    status=pd$Status,    
                    marker = pd$riskscore5,     
                    predict.time = 12, method="KM")
troc.9= survivalROC(Stime=pd$Time,  
                    status=pd$Status,  
                    marker = pd$riskscore9,     
                    predict.time = 12, method="KM")
troc.10= survivalROC(Stime=pd$Time,  
                     status=pd$Status,    
                     marker = pd$riskscore10,     
                     predict.time = 12, method="KM")
png("CRG/result/figs/10_rocall.tiff",width = 1500,height = 1500,res=200)
plot(troc$FP, troc$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=3,cex.axis=1.5,ann = F)
title(xlab="1-Specificity", ylab="Sensitivity",cex.lab=1.5,font.lab=4)
abline(0,1,col='gray',lty=2,lwd=2)
lines(troc.1$FP, troc.1$TP, type="l",lwd=2,col="blue", xlim=c(0,1), ylim=c(0,1))
lines(troc.2$FP, troc.2$TP, type="l",lwd=2,col="green", xlim=c(0,1), ylim=c(0,1))
lines(troc.5$FP, troc.5$TP, type="l",lwd=2,col="pink", xlim=c(0,1), ylim=c(0,1))
lines(troc.9$FP, troc.9$TP, type="l",lwd=2,col="orange", xlim=c(0,1), ylim=c(0,1))
lines(troc.10$FP, troc.10$TP,lwd=2, type="l",col="purple", xlim=c(0,1), ylim=c(0,1))
legend(0.45,0.3,c(paste('Risk=',round(troc$AUC,3)),
                  paste('F2=',round(troc.1$AUC,3)),
                  paste('P2RX1=',round(troc.2$AUC,3)),
                  paste('PRKCZ=',round(troc.5$AUC,3)),
                  paste('CR2=',round(troc.9$AUC,3)),
                  paste('F12=',round(troc.10$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty=1,lwd=2,col = c('red','blue',"green","pink","orange","purple"),
       bty='n',seg.len = 2,cex = 0.78)
dev.off()


#######################################GSE68465
exp <- exp443 %>% t %>% as.data.frame()
exp <- exp[,colnames(exp)%in%gene]
exp$sample <- rownames(exp)
clin <- merge(clindata443,exp,by.x='sample')
clin$Status <- ifelse(clin$status=='NA',2,ifelse(clin$status=='Alive',0,1))
clin <- clin[clin$Status!=2,]
clin <- clin[clin$time!='NA',]
clin$Time <- as.numeric(clin$time)
clin <- na.omit(clin)
clin$riskscore <- 0.08691*clin$F2 -0.31064*clin$P2RX1 + -0.13589*clin$PRKCZ -0.02263*clin$C4BPA -0.01604*clin$CR2  
clin$Group <- ifelse(clin$riskscore>median(clin$riskscore),'High risk','Low risk')
library(survival)
library(survminer)
ggsurvplot(survfit(Surv(Time, Status) ~ Group, data=clin),surv.median.line = "hv",palette = c('red','blue','grey'),pval = T,pval.size=10)

library(ggstatsplot)
library(tidyverse)
#风险评分box图
# ggbetweenstats(data = clin,x = Status,y = riskscore,plot.type = "box",type = "np",p.adjust.method = "fdr")
cb <- ggboxplot(clin,x = 'status',y = 'riskscore',fill='status')+theme_classic()+ylab('Riskscore')+xlab('')+geom_signif(comparisons = list(c("Dead", "Alive")),test = wilcox.test,map_signif_level = F)
ggsave('CRG/result/figs/11_box443.tiff',cb,width = 5,height = 5,dpi=1000)
#点图
pd <- clin[order(clin$riskscore,decreasing = F),]
pd$Sample <- 1:nrow(pd)
sp <- ggplot(pd,aes(Sample,Time,color=status))+geom_point()+geom_vline(aes(xintercept = 221),lty=2)+theme_bw()+scale_colour_manual(values = c('blue','red'))
rp <- ggplot(pd,aes(Sample,riskscore,color=Group))+geom_point()+geom_vline(aes(xintercept = 221),lty=2)+theme_bw()+scale_colour_manual(values = c('red','blue'))
ggsave('CRG/result/figs/11_survpoint443.tiff',sp,width = 5,height =2,dpi=1000)
ggsave('CRG/result/figs/11_riskpoint443.tiff',rp,width = 5,height = 3,dpi=1000)

dat <- pd[,8:12]
rownames(dat) <- pd$sample
library(pheatmap)
group=pd$Group%>%as.data.frame()
colnames(group) <- 'group'
pheatmap(t(dat) ,show_colnames = F,cluster_cols = F)
library(writexl)
writexl::write_xlsx(pd[,c('Time','Status','riskscore')],'CRG/result/files/02_survroc443.xlsx')
save(pd,file = 'CRG/extdata/GSE68465_surv443.Rdata')
library(survivalROC)
troc= survivalROC(Stime=pd$Time,  
                  status=pd$Status,   
                    marker = pd$riskscore,     
                  predict.time = 12, method="KM")
troc.3= survivalROC(Stime=pd$Time,  
                  status=pd$Status,   
                  marker = pd$riskscore,     
                  predict.time = 36, method="KM")
troc.5= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$riskscore,     
                    predict.time = 60, method="KM")
plot(troc$FP, troc$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1),   
     xlab="1-Specificity", 
     ylab="Sensitivity",lwd=2)
abline(0,1,col='gray',lty=2)
lines(troc.3$FP, troc.3$TP, type="l",col="blue", xlim=c(0,1), ylim=c(0,1),lwd=2)
lines(troc.5$FP, troc.5$TP, type="l",col="green", xlim=c(0,1), ylim=c(0,1),lwd=2)
legend(0.45,0.3,c(paste('AUC of 1-year =',round(troc$AUC,3)),paste('AUC of 3-year =',round(troc.3$AUC,3)),paste('AUC of 5-year =',round(troc.5$AUC,3))),
       x.intersp =2,y.intersp = 2,lty = 1,lwd = 2,col = c('red','blue','green'),bty='n',seg.len = 2,cex = 1)


#GSE31210
exp <- exp226 %>% t %>% as.data.frame()
exp <- exp[,colnames(exp)%in%gene]
exp$sample <- rownames(exp)
clin <- merge(clindata226,exp,by.x='sample')
clin$Status <- ifelse(clin$status=='NA',2,ifelse(clin$status=='alive',0,1))
clin <- clin[clin$Status!=2,]
clin <- clin[clin$time!='NA',]
clin$Time <- round(as.numeric(gsub(';.*','',clin$time)),2)
clin <- na.omit(clin)
clin$riskscore <- 0.08691*clin$F2 -0.31064*clin$P2RX1 -0.13589*clin$PRKCZ -0.02263*clin$C4BPA  
clin$Group <- ifelse(clin$riskscore>median(clin$riskscore),'High risk','Low risk')
library(survival)
library(survminer)
ggsurvplot(survfit(Surv(Time, Status) ~ Group, data=clin),surv.median.line = "hv",palette = c('red','blue','grey'),pval = T,pval.size=10)

library(ggstatsplot)
library(tidyverse)
#风险评分box图
# ggbetweenstats(data = clin,x = Status,y = riskscore,plot.type = "box",type = "np",p.adjust.method = "fdr")
cb <- ggboxplot(clin,x = 'status',y = 'riskscore',fill='status')+theme_classic()+ylab('Riskscore')+xlab('')+geom_signif(comparisons = list(c("Dead", "Alive")),test = wilcox.test,map_signif_level = F)
ggsave('CRG/result/figs/12_box226.tiff',cb,width = 5,height = 5,dpi=1000)
#点图
pd <- clin[order(clin$riskscore,decreasing = F),]
pd$Sample <- 1:nrow(pd)
sp <- ggplot(pd,aes(Sample,Time,color=status))+geom_point()+geom_vline(aes(xintercept = 113),lty=2)+theme_bw()+scale_colour_manual(values = c('blue','red'))
rp <- ggplot(pd,aes(Sample,riskscore,color=Group))+geom_point()+geom_vline(aes(xintercept = 113),lty=2)+theme_bw()+scale_colour_manual(values = c('red','blue'))
ggsave('CRG/result/figs/12_survpoint226.tiff',sp,width = 5,height =2,dpi=1000)
ggsave('CRG/result/figs/12_riskpoint226.tiff',rp,width = 5,height = 3,dpi=1000)

dat <- pd[,8:11]
rownames(dat) <- pd$sample
library(pheatmap)
group=pd$Group%>%as.data.frame()
colnames(group) <- 'group'
pheatmap(t(dat) ,show_colnames = F,cluster_cols = F)
library(writexl)
writexl::write_xlsx(pd[,c('Time','Status','riskscore')],'CRG/result/files/03_survroc226.xlsx')
save(pd,file = 'CRG/extdata/GSE31210_surv226.Rdata')
library(survivalROC)
troc= survivalROC(Stime=pd$Time,  
                  status=pd$Status,   
                  marker = pd$riskscore,     
                  predict.time = 12, method="KM")
troc.3= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$riskscore,     
                    predict.time = 36, method="KM")
troc.5= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$riskscore,     
                    predict.time = 60, method="KM")
plot(troc$FP, troc$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1),   
     xlab="1-Specificity", 
     ylab="Sensitivity",lwd=2)
abline(0,1,col='gray',lty=2)
lines(troc.3$FP, troc.3$TP, type="l",col="blue", xlim=c(0,1), ylim=c(0,1),lwd=2)
lines(troc.5$FP, troc.5$TP, type="l",col="green", xlim=c(0,1), ylim=c(0,1),lwd=2)
legend(0.45,0.3,c(paste('AUC of 1-year =',round(troc$AUC,3)),paste('AUC of 3-year =',round(troc.3$AUC,3)),paste('AUC of 5-year =',round(troc.5$AUC,3))),
       x.intersp =2,y.intersp = 2,lty = 1,lwd = 2,col = c('red','blue','green'),bty='n',seg.len = 2,cex = 1)

#GSE42127
exp <- exp133 %>% t %>% as.data.frame()
exp <- exp[,colnames(exp)%in%gene]
exp$sample <- rownames(exp)
clin <- merge(clindata133,exp,by.x='sample')
clin$Status <- ifelse(clin$status=='NA',2,ifelse(clin$status=='A',0,1))
clin <- clin[clin$Status!=2,]
clin <- clin[clin$time!='NA',]
clin$Time <- as.numeric(clin$time)
clin <- na.omit(clin)
clin <- as.data.frame(clin)
clin$riskscore <- 0.08691*as.numeric(clin$F12)-0.31064*as.numeric(clin$P2RX1) +0.05993*as.numeric(clin$SERPINE1) -0.02263*as.numeric(clin$C4BPA ) 
clin$Group <- ifelse(clin$riskscore>median(clin$riskscore),'High risk','Low risk')
library(survival)
library(survminer)
ggsurvplot(survfit(Surv(Time, Status) ~ Group, data=clin),surv.median.line = "hv",palette = c('red','blue','grey'),pval = T,pval.size=10)

library(ggstatsplot)
library(tidyverse)
#风险评分box图
# ggbetweenstats(data = clin,x = Status,y = riskscore,plot.type = "box",type = "np",p.adjust.method = "fdr")
cb <- ggboxplot(clin,x = 'status',y = 'riskscore',fill='status')+theme_classic()+ylab('Riskscore')+xlab('')+geom_signif(comparisons = list(c("Dead", "Alive")),test = wilcox.test,map_signif_level = F)
ggsave('CRG/result/figs/13_box133.tiff',cb,width = 5,height = 5,dpi=1000)
#点图
pd <- clin[order(clin$riskscore,decreasing = F),]
pd$Sample <- 1:nrow(pd)
sp <- ggplot(pd,aes(Sample,Time,color=status))+geom_point()+geom_vline(aes(xintercept = 66),lty=2)+theme_bw()+scale_colour_manual(values = c('blue','red'))
rp <- ggplot(pd,aes(Sample,riskscore,color=Group))+geom_point()+geom_vline(aes(xintercept = 66),lty=2)+theme_bw()+scale_colour_manual(values = c('red','blue'))
ggsave('CRG/result/figs/13_survpoint133.tiff',sp,width = 5,height =2,dpi=1000)
ggsave('CRG/result/figs/13_riskpoint133.tiff',rp,width = 5,height = 3,dpi=1000)

dat <- pd[,8:11]
rownames(dat) <- pd$sample
library(pheatmap)
group=pd$Group%>%as.data.frame()
colnames(group) <- 'group'
pheatmap(t(dat) ,show_colnames = F,cluster_cols = F)
library(writexl)
writexl::write_xlsx(pd[,c('Time','Status','riskscore')],'CRG/result/files/04_survroc133.xlsx')
save(pd,file = 'CRG/extdata/GSE42127_surv133.Rdata')
library(survivalROC)
troc= survivalROC(Stime=pd$Time,  
                  status=pd$Status,   
                  marker = pd$riskscore,     
                  predict.time = 12, method="KM")
troc.3= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$riskscore,     
                    predict.time = 36, method="KM")
troc.5= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$riskscore,     
                    predict.time = 60, method="KM")
plot(troc$FP, troc$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1),   
     xlab="1-Specificity", 
     ylab="Sensitivity",lwd=2)
abline(0,1,col='gray',lty=2)
lines(troc.3$FP, troc.3$TP, type="l",col="blue", xlim=c(0,1), ylim=c(0,1),lwd=2)
lines(troc.5$FP, troc.5$TP, type="l",col="green", xlim=c(0,1), ylim=c(0,1),lwd=2)
legend(0.45,0.3,c(paste('AUC of 1-year =',round(troc$AUC,3)),paste('AUC of 3-year =',round(troc.3$AUC,3)),paste('AUC of 5-year =',round(troc.5$AUC,3))),
       x.intersp =2,y.intersp = 2,lty = 1,lwd = 2,col = c('red','blue','green'),bty='n',seg.len = 2,cex = 1)

# GSE37745
exp <- exp106 %>% t %>% as.data.frame()
exp <- exp[,colnames(exp)%in%gene]
exp$sample <- rownames(exp)
clin <- merge(clindata106,exp,by.x='sample')
clin$Status <- ifelse(clin$status=='not known',2,ifelse(clin$status=='Alive',0,1))
clin <- clin[clin$Status!=2,]
clin <- clin[clin$time!='NA',]
clin$Time <- round(as.numeric(clin$time)/30,2)
clin <- na.omit(clin)
clin <- as.data.frame(clin)
clin$riskscore <- 0.08691*as.numeric(clin$F2)-0.31064*as.numeric(clin$P2RX1) -0.13589*as.numeric(clin$PRKCZ) -0.02263*as.numeric(clin$C4BPA ) 
clin$Group <- ifelse(clin$riskscore>median(clin$riskscore),'High risk','Low risk')
library(survival)
library(survminer)
ggsurvplot(survfit(Surv(Time, Status) ~ Group, data=clin),surv.median.line = "hv",palette = c('red','blue','grey'),pval = T,pval.size=10)

library(ggstatsplot)
library(tidyverse)
#风险评分box图
# ggbetweenstats(data = clin,x = Status,y = riskscore,plot.type = "box",type = "np",p.adjust.method = "fdr")
cb <- ggboxplot(clin,x = 'status',y = 'riskscore',fill='status')+theme_classic()+ylab('Riskscore')+xlab('')+geom_signif(comparisons = list(c("Dead", "Alive")),test = wilcox.test,map_signif_level = F)
ggsave('CRG/result/figs/14_box106.tiff',cb,width = 5,height = 5,dpi=1000)
#点图
pd <- clin[order(clin$riskscore,decreasing = F),]
pd$Sample <- 1:nrow(pd)
sp <- ggplot(pd,aes(Sample,Time,color=status))+geom_point()+geom_vline(aes(xintercept = 26),lty=2)+theme_bw()+scale_colour_manual(values = c('blue','red'))
rp <- ggplot(pd,aes(Sample,riskscore,color=Group))+geom_point()+geom_vline(aes(xintercept = 26),lty=2)+theme_bw()+scale_colour_manual(values = c('red','blue'))
ggsave('CRG/result/figs/14_survpoint106.tiff',sp,width = 5,height =2,dpi=1000)
ggsave('CRG/result/figs/14_riskpoint106.tiff',rp,width = 5,height = 3,dpi=1000)

dat <- pd[,7:10]
rownames(dat) <- pd$sample
library(pheatmap)
group=pd$Group%>%as.data.frame()
colnames(group) <- 'group'
pheatmap(t(dat) ,show_colnames = F,cluster_cols = F)
library(writexl)
writexl::write_xlsx(pd[,c('Time','Status','riskscore')],'CRG/result/files/05_survroc106.xlsx')
save(pd,file = 'CRG/extdata/GSE37745_surv106.Rdata')
library(survivalROC)
troc= survivalROC(Stime=pd$Time,  
                  status=pd$Status,   
                  marker = pd$riskscore,     
                  predict.time = 12, method="KM")
troc.3= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$riskscore,     
                    predict.time = 36, method="KM")
troc.5= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$riskscore,     
                    predict.time = 60, method="KM")
plot(troc$FP, troc$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1),   
     xlab="1-Specificity", 
     ylab="Sensitivity",lwd=2)
abline(0,1,col='gray',lty=2)
lines(troc.3$FP, troc.3$TP, type="l",col="blue", xlim=c(0,1), ylim=c(0,1),lwd=2)
lines(troc.5$FP, troc.5$TP, type="l",col="green", xlim=c(0,1), ylim=c(0,1),lwd=2)
legend(0.45,0.3,c(paste('AUC of 1-year =',round(troc$AUC,3)),paste('AUC of 3-year =',round(troc.3$AUC,3)),paste('AUC of 5-year =',round(troc.5$AUC,3))),
       x.intersp =2,y.intersp = 2,lty = 1,lwd = 2,col = c('red','blue','green'),bty='n',seg.len = 2,cex = 1)

