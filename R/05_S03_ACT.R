response <- read_xlsx('rawdata/Melanoma/GSE100797_ACT/luass.xlsx')
luass <- read.table('rawdata/Melanoma/GSE100797_ACT/ICB.Lauss2017_ACT_Melanoma.self_subtract',header = T,row.names = 1)
clin <- read.table('rawdata/Melanoma/GSE100797_ACT/ICB.Lauss2017_ACT_Melanoma.clinical',header = T,row.names = 1)
response$Group <- ifelse(response$RECIST%in%c('PD','SD'),'SD/PD',ifelse(response$RECIST%in%c('CR','PR'),'CR/PR','NE'))
response <- response[response$Group!='NE',]
colnames(response)[1] <- 'sample'
response$sample <- gsub('MM909_','p',response$sample)

luass$ENTREZID <- rownames(luass)
df <- bitr( rownames(luass), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
luass <- merge(luass ,df,by.x='ENTREZID')
rownames(luass) <- luass$SYMBOL
bsrgenes <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','ITGB1','C4BPA','JMJD7-PLA2G4B','CR2')
bsrgenes%in%df$SYMBOL
exp <- luass[luass$SYMBOL%in%bsrgenes,-c(1,27)]
exp1 <- exp %>% t %>% as.data.frame()
exp1$Risk=0.083*exp1$F2-0.256*exp1$P2RX+0.127*exp1$PLAUR-0.163*exp1$PRKCZ+0.098*exp1$F12+0.083*exp1$SERPINE1+
  0.076*exp1$ITGB1-0.029*exp1$C4BPA-0.021*exp1$CR2-0.055*exp1$GUCY1A1
exp1$sample <- rownames(exp1)
clin$sample <- rownames(clin)
survdata <- merge(exp1,clin,by.x='sample')
survdata$group <- ifelse(survdata$Risk>median(survdata$Risk),'High risk','Low risk')
response1 <- response[,c(1,15)]
survdata1 <- merge(survdata,response1,by.x='sample')
# Fig5D_ACT
write.xlsx(survdata1[,c(1,12:14,19)],file = 'files/Fig5.xlsx',sheetName = 'Fig5D_ACT',row.names = F,append = T)
png("figures/05D_PFSKM_risk_GSE100797_ACT.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(PFS, PFS.Event) ~ group, data=survdata),size=3,
                surv.median.line = "hv",
                title = "GSE100797_ACT",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of PFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

# SupFig3A_ACT
write.xlsx(survdata1[,c(1,13:14,20)],file = 'files/SupFig3.xlsx',sheetName = 'SupFig3A_ACT',row.names = F)
png("figures/S03A_PFSKM_risk_response_ACT.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(PFS, PFS.Event) ~ Group, data=survdata1),size=3,
                surv.median.line = "hv",
                title = "GSE100797_ACT",
                legend.labs = c('CR/PR','SD/PD'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of PFS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

# Fig5C_ACT
library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
library(reshape2)
a <-data.frame(group=c('High risk','Low risk'),CR_PR=c(t(table(survdata1$Group[survdata1$group=='High risk'])/length(survdata1$Group[survdata1$group=='High risk']))[1],
                                                       t(table(survdata1$Group[survdata1$group=='Low risk'])/length(survdata1$Group[survdata1$group=='Low risk']))[1]),
               SD_PD=c(t(table(survdata1$Group[survdata1$group=='High risk'])/length(survdata1$Group[survdata1$group=='High risk']))[2],
                       t(table(survdata1$Group[survdata1$group=='Low risk'])/length(survdata1$Group[survdata1$group=='Low risk']))[2]))
b=melt(a,id.vars='group')
b$variable <- gsub('_','/',b$variable)
write.xlsx(b,file = 'files/Fig5.xlsx',sheetName = 'Fig5C_ACT',row.names = F,append = T)
P <- ggplot(data=b,aes(x=group,y=value,fill=variable)) + theme_bw()+ylab('Percent')+xlab('')+ggtitle('GSE100797_ACT')+ylim(0,1)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_col(position = 'stack', width = 0.6)+
  geom_bar(position = "stack", stat = "identity", width = 0.6) +
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('figures/05C_barplot_GSE100797_ACT.png',P,width = 5,height = 5,dpi=200)

# SupFig3B_ACT
write.xlsx(survdata1[,c(1,12:14)],file = 'files/SupFig3.xlsx',sheetName = 'SupFig3B_ACT',row.names = F,append = T)
library(survivalROC)
pd=survdata1
troc.1= survivalROC(Stime=pd$PFS,  
                    status=pd$PFS.Event,   
                    marker = pd$Risk,     
                    predict.time = 12, method="KM")
troc.13= survivalROC(Stime=pd$PFS,  
                     status=pd$PFS.Event,   
                     marker = pd$Risk,       
                     predict.time = 36, method="KM")
troc.15= survivalROC(Stime=pd$PFS,  
                     status=pd$PFS.Event,   
                     marker = pd$Risk,      
                     predict.time = 60, method="KM")
png("figures/S03B_rocyear_risk_ACT.png",width = 1500,height = 1500,res=200)
plot(troc.1$FP, troc.1$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
title(main='GSE100797_ACT',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.13$FP, troc.13$TP, type="l",lwd=4,col="orange", xlim=c(0,1), ylim=c(0,1))
lines(troc.15$FP, troc.15$TP, type="l",lwd=4,col="blue", xlim=c(0,1), ylim=c(0,1))
legend(0.55,0.2,c(paste('AUC of 1-year =',round(troc.1$AUC,3)),
                  paste('AUC of 3-year =',round(troc.13$AUC,3)),
                  paste('AUC of 5-year =',round(troc.15$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()

# SupFig3C_ACT
write.xlsx(survdata1[,c(1,12,17)],file = 'files/SupFig3.xlsx',sheetName = 'SupFig3C_ACT',row.names = F,append = T)
library(pROC)
roc1 <- roc(pd$RECIST,pd$Risk)
png("figures/S03C_roc_risk_ACT.png",width = 1500,height = 1500,res=200)
plot(roc1,col="red",lwd=4,cex.axis=1.5,ann = F, xlim=c(1,0), ylim=c(0,1))
title(xlab="Specificity", ylab="Sensitivity",cex.lab=1.6,font.lab=2)
legend(0.55,0.2,c(paste('AUC =',round(roc1$auc,3)) ),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()


