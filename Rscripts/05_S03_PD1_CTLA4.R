###################黑色素瘤PRJEB23709 PD1+CTLA4
response <- read_xlsx('rawdata/Melanoma/PRJEB23709/Gide2019_PD1+CTLA4_Melanoma_RNASeq/mmc3.xlsx')
response$Group <- ifelse(response$`Best RECIST response`%in%c('SD','PD'),'SD/PD',ifelse(response$`Best RECIST response`%in%c('CR','PR'),'CR/PR','NE'))
response <- response[response$Group!='NE',]
colnames(response)[1] <- 'sample'
response$PFS <- round(response$`Progression Free Survival (Days)`/30,2)
response$OS <- round(response$`Overall Survival (Days)`/30,2)
response$OS.E <- ifelse(response$`Last Followup Status`=='Alive',0,1)
response$PFS.E <- ifelse(response$`Last Followup Status`=='Alive',0,1)
response$PFS.E[response$PFS<response$OS] <- 1
response$sample <-paste0('X',response$sample) 
Gide <- read.table('rawdata/Melanoma/PRJEB23709/Gide2019_PD1+CTLA4_Melanoma_RNASeq/ICB.Gide2019_Pembrolizumab-Nivolumab+Ipilimumab_Melanoma.self_subtract',header = T,row.names = 1)
Gide$ENTREZID <- rownames(Gide)
df <- bitr( rownames(Gide), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
Gide<- merge(Gide ,df,by.x='ENTREZID')
rownames(Gide) <- Gide$SYMBOL
bsrgenes <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','ITGB1','C4BPA','JMJD7-PLA2G4B','CR2')
bsrgenes%in%df$SYMBOL
exp <- Gide[Gide$SYMBOL%in%bsrgenes,-c(1,ncol(Gide))]
exp1 <- exp %>% t %>% as.data.frame()
colnames(exp1) <- gsub('-','_',colnames(exp1))
exp1$Risk=0.083*exp1$F2-0.256*exp1$P2RX+0.127*exp1$PLAUR-0.163*exp1$PRKCZ+0.098*exp1$F12+0.083*exp1$SERPINE1+
  0.076*exp1$ITGB1-0.021*exp1$CR2-0.055*exp1$GUCY1A1-0.035*exp1$JMJD7_PLA2G4B
exp1$sample <- rownames(exp1)
survdata <- merge(exp1,response,by.x='sample')
survdata$group <- ifelse(survdata$Risk>median(survdata$Risk),'High risk','Low risk')

# Fig5D_PD1_CTLA4
write.xlsx(survdata[,c(1,12,28,31,32)],file = 'files/Fig5.xlsx',sheetName = 'Fig5D_PD1_CTLA4',row.names = F,append = T)
png("figures/05D_PFSKM_risk_PRJEB23709_PD1_CTLA4.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(PFS, PFS.E) ~ group, data=survdata),size=3,
                surv.median.line = "hv",
                title = "PRJEB23709_PD1_CTLA4",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of PFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()



# SupFig3A_PD1_CTLA4
write.xlsx(survdata[,c(1,28,31,27)],file = 'files/SupFig3.xlsx',sheetName = 'SupFig3A_PD1_CTLA4',row.names = F,append = T)
png("figures/S03A_PFSKM_risk_response_PD1_CTLA4.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(PFS, PFS.E) ~ Group, data=survdata),size=3,
                surv.median.line = "hv",
                title = "PRJEB23709_PD1_CTLA4",
                legend.labs = c('CR/PR','SD/PD'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of PFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()


# Fig5C_PD1_CTLA4
library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
library(reshape2)

a <-data.frame(group=c('High risk','Low risk'),CR_PR=c(t(table(survdata$Group[survdata$group=='High risk'])/length(survdata$Group[survdata$group=='High risk']))[1],
                                                       t(table(survdata$Group[survdata$group=='Low risk'])/length(survdata$Group[survdata$group=='Low risk']))[1]),
               SD_PD=c(t(table(survdata$Group[survdata$group=='High risk'])/length(survdata$Group[survdata$group=='High risk']))[2],
                       t(table(survdata$Group[survdata$group=='Low risk'])/length(survdata$Group[survdata$group=='Low risk']))[2]))
b=melt(a,id.vars='group')
b$variable <- gsub('_','/',b$variable)
write.xlsx(b,file = 'files/Fig5.xlsx',sheetName = 'Fig5C_PD1_CTLA4',row.names = F,append = T)
P <- ggplot(data=b,aes(x=group,y=value,fill=variable)) + theme_bw()+ylab('Percent')+xlab('')+ggtitle('PRJEB23709_PD1_CTLA4')+ylim(0,1)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_col(position = 'stack', width = 0.6)+
  geom_bar(position = "stack", stat = "identity", width = 0.6) +
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('figures/05C_barplot_PRJEB23709_PD1_CTLA4.png',P,width = 5.2,height = 5.2,dpi=200)

# SupFig3B_PD1_CTLA4
write.xlsx(survdata[,c(1,12,28,31)],file = 'files/SupFig3.xlsx',sheetName = 'SupFig3B_PD1_CTLA4',row.names = F,append = T)
library(survivalROC)
pd=survdata
troc.1= survivalROC(Stime=pd$PFS,  
                    status=pd$PFS.E,   
                    marker = pd$Risk,     
                    predict.time = 6, method="KM")
troc.12= survivalROC(Stime=pd$PFS,  
                     status=pd$PFS.E,   
                     marker = pd$Risk,       
                     predict.time = 12, method="KM")
troc.13= survivalROC(Stime=pd$PFS,  
                     status=pd$PFS.E,   
                     marker = pd$Risk,      
                     predict.time = 24, method="KM")
png("figures/S03B_rocyear_risk_PD1_CTLA4.png",width = 1500,height = 1500,res=200)
plot(troc.1$FP, troc.1$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
title(main='PRJEB23709_PD1_CTLA4',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.12$FP, troc.12$TP, type="l",lwd=4,col="orange", xlim=c(0,1), ylim=c(0,1))
lines(troc.13$FP, troc.13$TP, type="l",lwd=4,col="blue", xlim=c(0,1), ylim=c(0,1))
legend(0.55,0.2,c(paste('AUC of 6-month =',round(troc.1$AUC,3)),
                  paste('AUC of 1-year =',round(troc.12$AUC,3)),
                  paste('AUC of 2-year =',round(troc.13$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()

# SupFig3C_PD1_CTLA4
write.xlsx(survdata[,c(1,12,27)],file = 'files/SupFig3.xlsx',sheetName = 'SupFig3C_PD1_CTLA4',row.names = F,append = T)
library(pROC)
roc1 <- roc(pd$`Best RECIST response`,pd$Risk)
png("figures/S03C_roc_risk_PD1-CTLA4.png",width = 1500,height = 1500,res=200)
plot(roc1,col="red",lwd=4,cex.axis=1.5,ann = F, xlim=c(1,0), ylim=c(0,1))
title(xlab="Specificity", ylab="Sensitivity",cex.lab=1.6,font.lab=2)
legend(0.55,0.2,c(paste('AUC =',round(roc1$auc,3)) ),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)
dev.off()

