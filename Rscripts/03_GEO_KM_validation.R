library(dplyr)
###########################################FUSCC
clin <- modeldat[,c('RNA.LC',unigenes,'OS','OS.E')]
clin$Risk <- predict(mulcox,clin)
colnames(clin) <- c('sample',unigenes,'Time','Status','Risk')
res.cut <- surv_cutpoint(clin,
                         time='Time',
                         event = 'Status',
                         variables = 'Risk')
summary(res.cut)
clin$Group <- ifelse(clin$Risk> 1.911378,'High risk','Low risk')

library(xlsx)
write.xlsx(clin,file = 'figures/02_data.xlsx',sheetName = 'FUSCC',row.names = F,append = T)
library(survival)
library(survminer)
png("figures/02C_OSKM_risk_FUSCC.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(Time, Status) ~ Group, data=clin),
                surv.median.line = "hv",size=3,
                title = "FUSCC(n=99)",
                legend.labs = c('High-risk','Low-risk'),legend=c(0.8,1),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Overall survival',xlab='Time(months)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

clin <- modeldat[,c('RNA.LC',unigenes,'RFS','RFS.E')]
clin$Risk <- predict(mulcox,clin)
colnames(clin) <- c('sample',unigenes,'Time','Status','Risk')
res.cut <- surv_cutpoint(clin,
                         time='Time',
                         event = 'Status',
                         variables = 'Risk')
summary(res.cut)
clin$Group <- ifelse(clin$Risk> 0.7075708,'High risk','Low risk')

library(xlsx)
write.xlsx(clin,file = 'figures/02_data.xlsx',sheetName = 'FUSCC_RFS',row.names = F,append = T)
library(survival)
library(survminer)
png("figures/02C_RFSKM_risk_FUSCC.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(Time, Status) ~ Group, data=clin),
                surv.median.line = "hv",size=3,
                title = "FUSCC(n=99)",
                legend.labs = c('High-risk','Low-risk'),legend=c(0.8,1),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Recurrence-free survival',xlab='Time(months)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

troc.1= survivalROC(Stime=clin$Time,  
                    status=clin$Status,   
                    marker = clin$Risk,     
                    predict.time = 12, method="KM")
troc.2= survivalROC(Stime=clin$Time,  
                     status=clin$Status,   
                     marker = clin$Risk,     
                     predict.time = 24, method="KM")
troc.3= survivalROC(Stime=clin$Time,  
                     status=clin$Status,   
                     marker = clin$Risk,     
                     predict.time = 36, method="KM")

png("figures/04C_rocyear_risk_FUSCC.png",width = 1500,height = 1500,res=200)
plot(troc.1$FP, troc.1$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=4,cex.axis=1.5,ann = F)
title(main='',xlab="1-Specificity", ylab="Sensitivity",cex.main = 2.5, font.main= 2,cex.lab=1.6,font.lab=2)
abline(0,1,col='gray',lty=2,lwd=3)
lines(troc.2$FP, troc.2$TP, type="l",lwd=4,col="orange", xlim=c(0,1), ylim=c(0,1))
lines(troc.3$FP, troc.3$TP, type="l",lwd=4,col="blue", xlim=c(0,1), ylim=c(0,1))
legend(0.55,0.2,c(paste('AUC of 1-year =',round(troc.1$AUC,3)),
                  paste('AUC of 2-year =',round(troc.2$AUC,3)),
                  paste('AUC of 3-year =',round(troc.3$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()

clin99 <- clin
clin <- modeldat[,c('RNA.LC','age','sex','pT','pN','Pathologic_stage','TP53','EGFR')]
colnames(clin) <- c('sample','Age','Gender','pT','pN','Stage','TP53','EGFR')
clin <- merge(clin99,clin)
clin <- clin[,-c(2:6,10)]
clin$Age <- ifelse(clin$Age>70,'old','young')

colnames(clin) <- c('sample','Time','Status','CRRS','Age','Gender','pT','pN','Stage','TP53','EGFR')
covariates <- colnames(clin)[4:ncol(clin)]
formulas <- sapply(covariates,function(x) as.formula(paste('Surv(Time, Status)~', x)))
models <- lapply( formulas, function(x){coxph(x, data = clin,x=T,y=T)})
results <- lapply(models,function(x){ 
  Cindex <- x$concordance[6]
  std <- x$concordance[7]
  res <- c(Cindex,std)
  names(res)<-c("Cindex",'std')
  return(res)
})
res <- t(as.data.frame(results, check.names = FALSE))
res <-as.data.frame(res,stringsAsFactors=F)
res$signatures <- rownames(res)
res$signatures <- factor(res$signatures,levels=res$signatures)
discrete_palettes <- RColorBrewer::brewer.pal(n=nrow(res), name='Spectral')
p=ggplot(res,aes(x=signatures,y=Cindex,fill = signatures))+geom_bar(stat = "identity")+
  theme_bw()+ylab('C-index(Compared with CRRS)')+xlab('')+ggtitle('')+
  geom_errorbar(aes(ymin=Cindex-std,ymax=Cindex+std),width=0.2,size=1)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),legend.position = 'none',axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=18,face='bold'),panel.grid = element_blank())+scale_fill_manual(values = discrete_palettes)
p

ggsave('figures/04_cindex_fuscc_compare_clinic.png',p,width = 10,height = 10,dpi=200)

res$cohort=rep('FUSCC',length(res$signatures))
resdf <- rbind(resdf,res)
library(xlsx)
write.xlsx(resdf,'figures/03_data.xlsx',row.names = F,sheet='cindex_compare',append = T)
########################################GSE72094
load('Rdata/GSE72094.Rdata')
unigenes <-c("C4BPA" , "COL1A2" ,"F2" ,"P2RX1" , "PLAUR" )
gene <- unigenes[unigenes%in%rownames(exp398)]


exp <- scale(exp398) %>% t %>% as.data.frame()
exp <- exp[,colnames(exp)%in%gene]
exp$sample <- rownames(exp)
clin <- merge(clindata398,exp,by.x='sample')
clin$Status <- ifelse(clin$status=='Alive',0,1)
clin$Time <- round(as.numeric(clin$time)/30,2)
clin$Risk=predict(mulcox,clin)
res.cut <- surv_cutpoint(clin,
                         time='Time',
                         event = 'Status',
                         variables = 'Risk')
summary(res.cut)
clin$Group <- ifelse(clin$Risk>0.1155871,'High risk','Low risk')
# SupFig2C_GSE72094
library(xlsx)
write.xlsx(clin,file = 'figures/02_data.xlsx',sheetName = 'GSE72094',row.names = F,append = T)
library(survival)
library(survminer)
png("figures/02C_OSKM_risk_GSE72094.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(Time, Status) ~ Group, data=clin),
                surv.median.line = "hv",size=3,
                title = "GSE72094(n=398)",
                legend.labs = c('High-risk','Low-risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Overall survival',xlab='Time(months)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()


########################################GSE13213
load("Rdata/GSE13213.Rdata")
unigenes <-c("C4BPA" , "COL1A2" ,"F2" ,"P2RX1" , "PLAUR" )
gene <- unigenes[unigenes%in%rownames(exp117)]
boxplot( scale(exp117)[1:100,1:20])
exp <- scale(exp117) %>% t %>% as.data.frame()
exp <- exp[,colnames(exp)%in%gene]
exp$sample <- rownames(exp)
clin <- merge(clindata117,exp,by.x='sample')
clin$Status <- ifelse(clin$status=='Alive',0,1)
clin$Time <- round(as.numeric(clin$time)/30,2)
clin$P2RX1[is.na(clin$P2RX1)] <- 0

clin$Risk=predict(mulcox,clin)
res.cut <- surv_cutpoint(clin,
                         time='Time',
                         event = 'Status',
                         variables = 'Risk')
summary(res.cut)

clin$Group <- ifelse(as.numeric(clin$Risk)>-2.271746,'High risk','Low risk')
# SupFig2C_GSE13213
library(xlsx)
write.xlsx(clin,file = 'figures/02_data.xlsx',sheetName = 'GSE13213',row.names = F,append = T)
library(survival)
library(survminer)
png("figures/02C_OSKM_risk_GSE13213.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(Time, Status) ~ Group, data=clin),
                surv.median.line = "hv",size=3,
                title = "GSE13213(n=117)",
                legend.labs = c('High-risk','Low-risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Overall survival',xlab='Time(months)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()
########################################GSE31210
load("Rdata/GSE31210.Rdata")

gene <- unigenes[unigenes%in%rownames(exp226)]

exp <- scale(exp226) %>% t %>% as.data.frame()
exp <- exp[,colnames(exp)%in%gene]
exp$sample <- rownames(exp)
clin <- merge(clindata226,exp,by.x='sample')
clin$Status <- ifelse(clin$status=='Alive',0,1)
clin$time <- as.numeric(clin$time)
clin$Risk= predict(mulcox,clin)
clin$Risk=predict(mulcox,clin)
res.cut <- surv_cutpoint(clin,
                         time='time',
                         event = 'Status',
                         variables = 'Risk')
summary(res.cut)
clin$Group <- ifelse(as.numeric(clin$Risk)>-0.1973602,'High risk','Low risk')
# SupFig2C_GSE31210
library(xlsx)
write.xlsx(clin,file = 'figures/02_data.xlsx',sheetName = 'GSE31210',row.names = F,append = T)
library(survival)
library(survminer)
png("figures/02C_OSKM_risk_GSE31210.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(time, Status) ~ Group, data=clin),
                surv.median.line = "hv",size=3,
                title = "GSE31210(n=226)",
                legend.labs = c('High-risk','Low-risk'),legend=c(0.8,1),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Overall survival',xlab='Time(months)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

########################################GSE30219
load('../Coagulation_project/Rdata/GSE30219.Rdata')

unigenes <-c("C4BPA" , "COL1A2" ,"F2" ,"P2RX1" , "PLAUR" )
gene <- unigenes[unigenes%in%rownames(exp293)]


exp <- scale(exp293) %>% t %>% as.data.frame()
exp <- exp[,colnames(exp)%in%gene]
exp$sample <- rownames(exp)
clin <- merge(clindata293,exp,by.x='sample')
clin$Time <- as.numeric(clin$OS)
clin$Risk=predict(mulcox,clin)
res.cut <- surv_cutpoint(clin,
                         time='Time',
                         event = 'Status',
                         variables = 'Risk')
summary(res.cut)
clin$Group <- ifelse(clin$Risk>0.3178462,'High risk','Low risk')
# GSE30219
library(xlsx)
write.xlsx(clin,file = 'figures/02_data.xlsx',sheetName = 'GSE30219',row.names = F,append = T)
library(survival)
library(survminer)
png("figures/02C_OSKM_risk_GSE30219.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(Time, Status) ~ Group, data=clin),
                surv.median.line = "hv",size=3,
                title = "GSE30219(n=293)",
                legend.labs = c('High-risk','Low-risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Overall survival',xlab='Time(months)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

########################################GSE3141
load('../Coagulation_project/Rdata/GSE3141.Rdata')

unigenes <-c("C4BPA" , "COL1A2" ,"F2" ,"P2RX1" , "PLAUR" )
gene <- unigenes[unigenes%in%rownames(exp58)]


exp <- scale(exp58) %>% t %>% as.data.frame()
exp <- exp[,colnames(exp)%in%gene]
exp$sample <- rownames(exp)

clin <- merge(clindata58,exp,by.x='sample')
clin$Time <- as.numeric(clin$OS)
clin$Status <- as.numeric(clin$OS.E)
clin$Risk=predict(mulcox,clin)
res.cut <- surv_cutpoint(clin,
                         time='Time',
                         event = 'Status',
                         variables = 'Risk')
summary(res.cut)
clin$Group <- ifelse(clin$Risk>3.028203 ,'High risk','Low risk')
# GSE3141
library(xlsx)
write.xlsx(clin,file = 'figures/02_data.xlsx',sheetName = 'GSE3141',row.names = F,append = T)
library(survival)
library(survminer)
png("figures/02C_OSKM_risk_GSE3141.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(Time, Status) ~ Group, data=clin),
                surv.median.line = "hv",size=3,
                title = "GSE3141(n=58)",
                legend.labs = c('High-risk','Low-risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Overall survival',xlab='Time(months)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()
########################################GSE68465
load('../Coagulation_project/Rdata/GSE68465.Rdata')

unigenes <-c("C4BPA" , "COL1A2" ,"F2" ,"P2RX1" , "PLAUR" )
gene <- unigenes[unigenes%in%rownames(exp433)]


exp <- scale(exp433) %>% t %>% as.data.frame()
exp <- exp[,colnames(exp)%in%gene]
exp$sample <- rownames(exp)
clin <- merge(clindata433,exp,by.x='sample')
clin$Time <- as.numeric(clin$OS)

clin$Risk=predict(mulcox,clin)
res.cut <- surv_cutpoint(clin,
                         time='Time',
                         event = 'Status',
                         variables = 'Risk')
summary(res.cut)
clin$Group <- ifelse(clin$Risk>1.205299 ,'High risk','Low risk')
# GSE68465
library(xlsx)
write.xlsx(clin,file = 'figures/02_data.xlsx',sheetName = 'GSE68465',row.names = F,append = T)
library(survival)
library(survminer)
png("figures/02C_OSKM_risk_GSE68465.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(Time, Status) ~ Group, data=clin),
                surv.median.line = "hv",size=3,
                title = "GSE68465(n=433)",
                legend.labs = c('High-risk','Low-risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Overall survival',xlab='Time(months)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()

########################################meta-cohort
gene <-c("C4BPA" , "COL1A2" ,"F2" ,"P2RX1" , "PLAUR" )
TCGA <- lxdata1[,gene]%>%t%>%as.data.frame()
colnames(TCGA) <- lxdata1$sampleID
GSE13213 <- exp117[gene,]
GSE31210 <- exp226[gene,]
GSE70294 <- exp398[gene,]
GSE30219 <- exp293[gene,]
GSE68465 <- exp433[gene,]
GSE3141 <- exp58[gene,]
batch <- paste0('batch',rep(c(1,2,3,4,5,6,7),c(502,117,226,398,293,433,58)))
expr <- cbind(TCGA,GSE13213,GSE31210,GSE70294,GSE30219,GSE68465,GSE3141 )
library(sva)
expr_batch <- ComBat(dat=expr,batch = batch)

boxplot(expr_batch)
boxplot(expr)


exp <- expr_batch %>% t %>% as.data.frame()
exp$sample <- rownames(exp)
clin_TCGA <- lxdata1[,c('sampleID','OS','OS.E')]
colnames(clin_TCGA) <- c('sample','Time','Status')

clindata117$Time <- as.numeric(clindata117$time)/30
clindata117$Status<- ifelse(clindata117$status=='Dead',1,0)
clin117 <- clindata117[,c('sample','Time','Status')]

clindata226$Time <- as.numeric(clindata226$time)
clindata226$Status<- ifelse(clindata226$status=='Dead',1,0)
clin226<- clindata226[,c('sample','Time','Status')]

clindata398$Time <- as.numeric(clindata398$time)/30
clindata398$Status<- ifelse(clindata398$status=='Dead',1,0)
clin398 <- clindata398[,c('sample','Time','Status')]

clindata293$Time <- as.numeric(clindata293$OS)
clin293 <- clindata293[,c('sample','Time','Status')]

clindata433$Time <- as.numeric(clindata433$OS)
clin433 <- clindata433[,c('sample','Time','Status')]


clindata58$Time <- as.numeric(clindata58$OS)
clindata58$Status <- as.numeric(clindata58$OS.E)
clin58 <- clindata58[,c('sample','Time','Status')]


clin2027 <- rbind(clin_TCGA,clin117,clin226,clin398,clin293,clin433,clin58)
clin <- merge(clin2027,exp,by.x='sample')

save(expr_batch,clin2027,file = "Rdata/Meta.Rdata")
clin$Risk=predict(mulcox,clin)
res.cut <- surv_cutpoint(clin,
                         time='Time',
                         event = 'Status',
                         variables = 'Risk')
summary(res.cut)
clin$Group <- ifelse(clin$Risk>2.379753  ,'High risk','Low risk')
# GSE68465
library(xlsx)
write.xlsx(clin,file = 'figures/02_data.xlsx',sheetName = 'Meta',row.names = F,append = T)
library(survival)
library(survminer)
png("figures/02C_OSKM_risk_Meta.png",width = 1500,height = 1500,res=200)
p <- ggsurvplot(survfit(Surv(Time, Status) ~ Group, data=clin),
                surv.median.line = "hv",size=3,
                title = "Meta-Cohort(n=2027)",
                legend.labs = c('High-risk','Low-risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                palette = c("#E7B800", "#2E9FDF",'grey'),
                pval = T,pval.size=18,legend.title = "",
                ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0,face = 'bold')), 
                font.main = c(30, "bold"),
                ylab='Overall survival',xlab='Time(months)',
                font.x = c(30, "bold"),
                font.y = c(30, "bold"));p
dev.off()
###########################################TCGA

res.cut <- surv_cutpoint(clin,
                         time='Time',
                         event = 'Status',
                         variables = 'CRRS')
summary(res.cut)
clin$Group <- ifelse(clin$CRRS>0.192116  ,'High risk','Low risk')
library(xlsx)
write.xlsx(clin,file = 'figures/02_data.xlsx',sheetName = 'TCGA',row.names = F,append = T)


