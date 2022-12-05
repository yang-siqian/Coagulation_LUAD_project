###################黑色素瘤PRJEB23709 PD1
response <- read_xlsx('LUADdata/TCGA/Gide2019_PD1_Melanoma_RNASeq/mmc2.xlsx')
response$Group <- ifelse(response$`Best RECIST response`%in%c('SD','PD'),'SD/PD',ifelse(response$`Best RECIST response`%in%c('CR','PR'),'CR/PR','NE'))
response <- response[response$Group!='NE',]
colnames(response)[1] <- 'sample'
response$PFS <- round(response$`Progression Free Survival (Days)`/30,2)
response$OS <- round(response$`Overall Survival (Days)`/30,2)
response$OS.E <- ifelse(response$`Last Followup Status`=='Alive',0,1)
response$PFS.E <- ifelse(response$`Last Followup Status`=='Alive',0,1)
response$PFS.E[response$PFS<response$OS] <- 1
response$sample <-paste0('X',response$sample) 
Gide <- read.table('LUADdata/TCGA/Gide2019_PD1_Melanoma_RNASeq/ICB.Gide2019_Pembrolizumab-Nivolumab_Melanoma.self_subtract',header = T,row.names = 1)
Gide$ENTREZID <- rownames(Gide)
df <- bitr( rownames(Gide), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
Gide<- merge(Gide ,df,by.x='ENTREZID')
rownames(Gide) <- Gide$SYMBOL
bsrgenes <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7-PLA2G4B','CR2')
bsrgenes%in%df$SYMBOL
exp <- Gide[Gide$SYMBOL%in%bsrgenes,-c(1,ncol(Gide))]
exp1 <- exp %>% t %>% as.data.frame()
colnames(exp1) <- gsub('-','_',colnames(exp1))
exp1$Risk=-0.299*exp1$P2RX1-0.106*exp1$GUCY1A1+0.107*exp1$PLAUR-0.124*exp1$PRKCZ+0.093*exp1$F12+0.056*exp1$SERPINE1+
  0.137*exp1$COL1A2-0.031*exp1$JMJD7_PLA2G4B-0.017*exp1$CR2
exp1$sample <- rownames(exp1)
survdata <- merge(exp1,response,by.x='sample')
survdata$group <- ifelse(survdata$Risk>median(survdata$Risk),'High risk','Low risk')

png("CRG/result/figs/09_OSKM_risk_PRJEB23709_PD1.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(OS,OS.E) ~ group, data=survdata),size=3,
                surv.median.line = "hv",
                title = "PRJEB23709_PD1",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of PFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("CRG/result/figs/09_PFSKM_risk_PRJEB23709_PD1.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(PFS, PFS.E) ~ group, data=survdata),size=3,
                surv.median.line = "hv",
                title = "PRJEB23709_PD1",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of PFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("CRG/result/figs/09_OSKM_risk_response_PD1.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(OS,OS.E) ~ Group, data=survdata),size=3,
                surv.median.line = "hv",
                title = "PRJEB23709_PD1",
                legend.labs = c('CR/PR','SD/PD'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of PFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

png("CRG/result/figs/09_PFSKM_risk_response_PD1.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(PFS, PFS.E) ~ Group, data=survdata),size=3,
                surv.median.line = "hv",
                title = "PRJEB23709_PD1",
                legend.labs = c('CR/PR','SD/PD'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of PFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
#风险评分box图
survdata$Group <- factor(survdata$Group,levels = c("CR/PR", "SD/PD"))
cb <- ggboxplot(survdata,x = 'Group',y = 'Risk',fill='Group')+theme_bw()+ylab('Riskscore')+xlab('')+ggtitle('PRJEB23709_PD1')+ylim(-0.65,1)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_signif(size=1,textsize=12,comparisons = list(c("CR/PR", "SD/PD")),test = t.test,map_signif_level = T)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('CRG/result/figs/09_riskbox_PRJEB23709_PD1.tiff',cb,width = 5,height = 5,dpi=200)
library(reshape2)
table(survdata$Group[survdata$group=='High risk'])
table(survdata$Group[survdata$group=='Low risk'])
a <-data.frame(group=c('High risk','Low risk'),CR_PR=c(30,61.9),SD_PD=c(70,38.1))
b=melt(a,id.vars='group')
b$variable <- gsub('_','/',b$variable)
P <- ggplot(data=b,aes(x=group,y=value,fill=variable)) + theme_bw()+ylab('Percent')+xlab('')+ggtitle('PRJEB23709_PD1')+ylim(0,100)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_col(position = 'stack', width = 0.6)+
  geom_bar(position = "stack", stat = "identity", width = 0.6) +
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('CRG/result/figs/09_barplot_PRJEB23709_PD1.tiff',P,width = 6,height = 6,dpi=200)


library(survivalROC)
pd=survdata
troc.1= survivalROC(Stime=pd$PFS,  
                    status=pd$PFS.E,   
                    marker = pd$Risk,     
                    predict.time = 12, method="KM")
troc.12= survivalROC(Stime=pd$PFS,  
                     status=pd$PFS.E,   
                     marker = pd$Risk,       
                     predict.time = 24, method="KM")
troc.13= survivalROC(Stime=pd$PFS,  
                     status=pd$PFS.E,   
                     marker = pd$Risk,      
                     predict.time = 36, method="KM")
png("CRG/result/figs/09_rocyear_risk_PD1.tiff",width = 1500,height = 1500,res=200)
plot(troc.1$FP, troc.1$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
title(main='PRJEB23709_PD1',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.12$FP, troc.12$TP, type="l",lwd=4,col="orange", xlim=c(0,1), ylim=c(0,1))
lines(troc.13$FP, troc.13$TP, type="l",lwd=4,col="blue", xlim=c(0,1), ylim=c(0,1))
legend(0.55,0.2,c(paste('AUC of 1-year =',round(troc.1$AUC,3)),
                  paste('AUC of 2-year =',round(troc.12$AUC,3)),
                  paste('AUC of 3-year =',round(troc.13$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()


library(pROC)
roc1 <- roc(pd$Group,pd$Risk)
png("CRG/result/figs/09_roc_risk_PD1.tiff",width = 1500,height = 1500,res=200)
plot(roc1,col="red",lwd=4,cex.axis=1.5,ann = F, xlim=c(1,0), ylim=c(0,1))
title(xlab="Specificity", ylab="Sensitivity",cex.lab=1.6,font.lab=2)
legend(0.55,0.2,c(paste('AUC =',round(roc1$auc,3)) ),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)
dev.off()
