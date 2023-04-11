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
clin$Group <- ifelse(clin$Risk>median(clin$Risk),'High risk','Low risk')
# SupFig2C_GSE72094
library(xlsx)
write.xlsx(clin,file = 'files/SupFig2.xlsx',sheetName = 'SupFig2_GSE72094',row.names = F,append = T)
library(survival)
library(survminer)
png("figures/S02C_OSKM_risk_GSE72094.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(Time, Status) ~ Group, data=clin),
                surv.median.line = "hv",size=3,
                title = "GSE72094",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()


library(ggstatsplot)
library(tidyverse)
library(ggpubr)
#风险评分box图
# cb <- ggboxplot(clin,x = 'status',y = 'Risk',fill='status')+theme_bw()+ylab('Riskscore')+xlab('')+ggtitle('GSE72094')+ylim(-1,3)+
#   scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
#   geom_signif(size=1,textsize=12,comparisons = list(c("Dead", "Alive")),test = wilcox.test,map_signif_level = T)+
#   theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
#         axis.text.x = element_text(size=18,face='bold',color='black'),
#         axis.text.y = element_text(size=18,face='bold',color='black'),
#         text=element_text(size=25,face='bold'),
#         legend.title = element_blank(),legend.position = "none",
#         panel.grid = element_blank(),panel.border = element_rect(size=1.5))
# ggsave('figures/S02_riskbox.png',cb,width = 5,height = 5,dpi=200)

# SupFig2A_GSE72094
#点图
pd <- clin[order(clin$Risk,decreasing = F),]
pd$Sample <- 1:nrow(pd)
sp <- ggplot(pd,aes(Sample,Time,color=status))+geom_point()+geom_vline(aes(xintercept = 199),lty=2,lwd=1)+
  theme_bw()+scale_colour_manual(values = c('blue','red'))+ylab('Survival month')+xlab('')+
  theme(axis.text.x = element_text(size=16,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
rp <- ggplot(pd,aes(Sample,Risk,color=Group))+geom_point()+geom_vline(aes(xintercept = 199),lty=2,lwd=1)+
  theme_bw()+scale_colour_manual(values = c('red','blue'))+ylab('Riskscore')+xlab('')+ggtitle('GSE72094')+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),
        legend.title = element_blank(),legend.justification=c(0,1),legend.position=c(0.01,0.99),
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('figures/S02A_survpoint_GSE72094.png',sp,width = 5,height =3,dpi=1000)
ggsave('figures/S02A_riskpoint_GSE72094.png',rp,width = 5,height = 3,dpi=1000)

library(pheatmap)
dat <- pd[,8:16]
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

save_pheatmap_pdf(p, "figures/S02A_pheat_GSE72094.pdf")

# SupFig2B_GSE72094
library(survivalROC)
troc.1= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$Risk,     
                    predict.time = 12, method="KM")
troc.13= survivalROC(Stime=pd$Time,  
                     status=pd$Status,   
                     marker = pd$Risk,     
                     predict.time = 36, method="KM")
troc.15= survivalROC(Stime=pd$Time,  
                     status=pd$Status,   
                     marker = pd$Risk,     
                     predict.time = 60, method="KM")
png("figures/S02B_rocyear_risk_GSE72094.png",width = 1500,height = 1500,res=200)
plot(troc.1$FP, troc.1$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
title(main='GSE72094',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.13$FP, troc.13$TP, type="l",lwd=4,col="orange", xlim=c(0,1), ylim=c(0,1))
lines(troc.15$FP, troc.15$TP, type="l",lwd=4,col="blue", xlim=c(0,1), ylim=c(0,1))
legend(0.55,0.2,c(paste('AUC of 1-year =',round(troc.1$AUC,3)),
                  paste('AUC of 3-year =',round(troc.13$AUC,3)),
                  paste('AUC of 5-year =',round(troc.15$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()

##################################################
# Fig2G_GSE72094
pd$F2_model=0.108 *pd$F2
pd$P2RX1_model=-0.213 *pd$P2RX1
pd$PLAUR_model=0.145 *pd$PLAUR
pd$PRKCZ_model=-0.194 *pd$PRKCZ 
pd$F12_model=0.133 *pd$F12
pd$SERPINE1_model=0.105 *pd$SERPINE1
pd$ITGB1_model=0.189 *pd$ITGB1
pd$C4BPA_model=-0.084 *pd$C4BPA
pd$CR2_model=-0.109 *pd$CR2
write.xlsx(pd[,c(1,17:19,22:30)],file = 'files/Fig2.xlsx',sheetName = 'Fig2G_GSE72094',row.names = F,append = T)

troc.1= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$Risk,     
                    predict.time = 12, method="KM")
troc.2= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$F2_model,     
                    predict.time =12, method="KM")
troc.3= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$P2RX1_model,     
                    predict.time = 12, method="KM")
troc.5= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$PLAUR_model,     
                    predict.time = 12, method="KM")
troc.6= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$PRKCZ_model,     
                    predict.time = 12, method="KM")
troc.7= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$F12_model,     
                    predict.time = 12, method="KM")
troc.8= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$SERPINE1_model,     
                    predict.time = 12, method="KM")
troc.9= survivalROC(Stime=pd$Time,  
                    status=pd$Status,   
                    marker = pd$ITGB1_model,     
                    predict.time = 12, method="KM")
troc.10= survivalROC(Stime=pd$Time,  
                     status=pd$Status,   
                     marker = pd$C4BPA_model,     
                     predict.time = 12, method="KM")

troc.12= survivalROC(Stime=pd$Time,  
                     status=pd$Status,   
                     marker = pd$CR2_model,     
                     predict.time = 12, method="KM")



png("figures/02G_rocmethods_bsr_GSE72094.png",width = 1500,height = 1500,res=200)
plot(troc.1$FP, troc.1$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
title(main='GSE72094',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.2$FP, troc.2$TP, type="l",lwd=4,col="blue", xlim=c(0,1), ylim=c(0,1))
lines(troc.3$FP, troc.3$TP, type="l",lwd=4,col="green", xlim=c(0,1), ylim=c(0,1))
lines(troc.5$FP, troc.5$TP, type="l",lwd=4,col="orange", xlim=c(0,1), ylim=c(0,1))
lines(troc.6$FP, troc.6$TP, type="l",lwd=4,col="black", xlim=c(0,1), ylim=c(0,1))
lines(troc.7$FP, troc.7$TP, type="l",lwd=4,col="yellow", xlim=c(0,1), ylim=c(0,1))
lines(troc.8$FP, troc.8$TP, type="l",lwd=4,col="pink", xlim=c(0,1), ylim=c(0,1))
lines(troc.9$FP, troc.9$TP, type="l",lwd=4,col="cyan", xlim=c(0,1), ylim=c(0,1))
lines(troc.10$FP, troc.10$TP, type="l",lwd=4,col="magenta", xlim=c(0,1), ylim=c(0,1))
lines(troc.12$FP, troc.12$TP, type="l",lwd=4,col="red4", xlim=c(0,1), ylim=c(0,1))
legend(0.55,0.5,c(paste('Risk=',round(troc.1$AUC,3)),
                  paste('F2=',round(troc.2$AUC,3)),
                  paste('P2RX1=',round(troc.3$AUC,3)),
                  paste('PLAUR=',round(troc.5$AUC,3)),
                  paste('PRKCZ=',round(troc.6$AUC,3)),
                  paste('LF12=',round(troc.7$AUC,3)),
                  paste('SERPINE1=',round(troc.8$AUC,3)),
                  paste('ITGB1=',round(troc.9$AUC,3)),
                  paste('C4BPA=',round(troc.10$AUC,3)),
                  paste('CR2=',round(troc.12$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','blue',"green",'orange','black',"yellow","pink","cyan","magenta","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()


