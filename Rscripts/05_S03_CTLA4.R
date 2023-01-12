###############################################################黑色素瘤CTLA4
response <- read.table('rawdata/Melanoma/Nathanson_CTLA4/ICB.Nathanson2017_Ipilimumab_Melanoma_Post.clinical',header = T)
response$Group <- ifelse(response$Response=='0','SD/PD','CR/PR')
colnames(response)[1] <- 'sample'
Nathanson <- read.table('rawdata/Melanoma/Nathanson_CTLA4/ICB.Nathanson2017_Ipilimumab_Melanoma_Post.self_subtract',header = T,sep = '\t')
rownames(Nathanson) <- Nathanson$Entrez 
colnames(Nathanson)[1] <- "ENTREZID"
df <- bitr( rownames(Nathanson), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
Nathanson<- merge(Nathanson ,df,by.x='ENTREZID')
rownames(Nathanson) <- Nathanson$SYMBOL
bsrgenes <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7-PLA2G4B','CR2')
bsrgenes%in%df$SYMBOL
exp <- Nathanson[Nathanson$SYMBOL%in%bsrgenes,-c(1,ncol(Nathanson))]
exp1 <- exp %>% t %>% as.data.frame()
colnames(exp1) <- gsub('-','_',colnames(exp1))
exp1$Risk=0.09*exp1$F2-0.299*exp1$P2RX1-0.106*exp1$GUCY1A1+0.107*exp1$PLAUR-0.124*exp1$PRKCZ+0.093*exp1$F12+0.056*exp1$SERPINE1+
  0.137*exp1$COL1A2-0.022*exp1$C4BPA-0.031*exp1$JMJD7_PLA2G4B-0.017*exp1$CR2
exp1$sample <- rownames(exp1)
survdata <- merge(exp1,response,by.x='sample')
survdata$group <- ifelse(survdata$Risk>median(survdata$Risk),'High risk','Low risk')
survdata$OS <- survdata$OS*12

# Fig5D_CTLA4
write.xlsx(survdata[,c(1,13:15,19)],file = 'files/Fig5.xlsx',sheetName = 'Fig5D_CTLA4',row.names = F,append = T)
png("figures/05D_OSKM_risk_Nathanson_CTLA4.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(OS,OS.Event) ~ group, data=survdata),size=3,
                surv.median.line = "hv",
                title = "Nathanson_CTLA4",
                legend.labs = c('High risk','Low risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

# SupFig3A_CTLA4
write.xlsx(survdata[,c(1,14:15,18)],file = 'files/SupFig3.xlsx',sheetName = 'SupFig3A_CTLA4',row.names = F,append = T)
png("figures/S03A_OSKM_risk_response_CTLA4.tiff",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(OS,OS.Event) ~ Group, data=survdata),size=3,
                surv.median.line = "hv",
                title = "Nathanson_CTLA4",
                legend.labs = c('CR/PR','SD/PD'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E69F00", "#56B4E9",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_bw()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Survival probability of OS',xlab='Time(month)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()


# Fig5C_CTLA4
library(dplyr)
library(ggstatsplot)
library(tidyverse)
library(ggpubr)
library(reshape2)
table(survdata$Group[survdata$group=='High risk'])
table(survdata$Group[survdata$group=='Low risk'])
a <-data.frame(group=c('High risk','Low risk'),CR_PR=c(14.3,37.5),SD_PD=c(85.7,62.5))
b=melt(a,id.vars='group')
b$variable <- gsub('_','/',b$variable)
write.xlsx(b,file = 'files/Fig5.xlsx',sheetName = 'Fig5C_CTLA4',row.names = F,append = T)
P <- ggplot(data=b,aes(x=group,y=value,fill=variable)) + theme_bw()+ylab('Percent')+xlab('')+ggtitle('Nathanson_CTLA4')+ylim(0,100)+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  geom_col(position = 'stack', width = 0.6)+
  geom_bar(position = "stack", stat = "identity", width = 0.6) +
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=1.5))
ggsave('figures/05C_barplot_Nathanson_CTLA4.tiff',P,width = 5,height = 5,dpi=200)

# SupFig3B_CTLA4
write.xlsx(survdata[,c(1,13:15)],file = 'files/SupFig3.xlsx',sheetName = 'SupFig3B_CTLA4',row.names = F,append = T)
library(survivalROC)
pd=survdata
troc.1= survivalROC(Stime=pd$Time,  
                    status=pd$OS.Event,   
                    marker = pd$Risk,     
                    predict.time = 12, method="KM")
troc.2= survivalROC(Stime=pd$Time,  
                    status=pd$OS.Event,   
                    marker = pd$Risk,     
                    predict.time = 24, method="KM")
png("figures/S03B_rocyear_risk_CTLA4.tiff",width = 1500,height = 1500,res=200)
plot(troc.1$FP, troc.1$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
title(main='Nathanson_CTLA4',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.2$FP, troc.2$TP, type="l",lwd=4,col="orange", xlim=c(0,1), ylim=c(0,1))
legend(0.55,0.2,c(paste('AUC of 1-year =',round(troc.1$AUC,3)),
                  paste('AUC of 2-year =',round(troc.2$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()

# SupFig3C_CTLA4
write.xlsx(survdata[,c(1,13,18)],file = 'files/SupFig3.xlsx',sheetName = 'SupFig3C_CTLA4',row.names = F,append = T)
library(pROC)
roc1 <- roc(pd$Group,pd$Risk)
png("figures/S03C_roc_risk_CTLA4.tiff",width = 1500,height = 1500,res=200)
plot(roc1,col="red",lwd=4,cex.axis=1.5,ann = F, xlim=c(1,0), ylim=c(0,1))
title(xlab="Specificity", ylab="Sensitivity",cex.lab=1.6,font.lab=2)
legend(0.55,0.2,c(paste('AUC =',round(roc1$auc,3)) ),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)
dev.off()
