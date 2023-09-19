immunotherapy <- function(i){
  
  if(i==1){
    response <- read_xlsx('rawdata/Melanoma/GSE100797_ACT/luass.xlsx')
    luass <- read.table('rawdata/Melanoma/GSE100797_ACT/ICB.Lauss2017_ACT_Melanoma.self_subtract',header = T,row.names = 1)
    clin <- read.table('rawdata/Melanoma/GSE100797_ACT/ICB.Lauss2017_ACT_Melanoma.clinical',header = T,row.names = 1)
    response$group <- ifelse(response$RECIST%in%c('PD','SD'),'SD/PD',ifelse(response$RECIST%in%c('CR','PR'),'CR/PR','NE'))
    response <- response[response$group!='NE',]
    colnames(response)[1] <- 'sample'
    response$sample <- gsub('MM909_','p',response$sample)
    
    luass$ENTREZID <- rownames(luass)
    df <- bitr( rownames(luass), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
    luass <- merge(luass ,df,by.x='ENTREZID')
    rownames(luass) <- luass$SYMBOL
    exp <- luass[luass$SYMBOL%in%unigenes,-c(1,27)]
    exp <- scale(exp) %>% t %>% as.data.frame()
    exp$Risk=predict(mulcox,exp)
    exp$sample <- rownames(exp)
    clin$sample <- rownames(clin)
    survdata <- merge(exp,clin,by.x='sample')
    survdata$group <- ifelse(survdata$Response==1,'CR/PR','SD/PD')
    res.cut <- surv_cutpoint(survdata,
                             time='PFS',
                             event = 'PFS.Event',
                             variables = 'Risk')
    
    survdata$Group <- ifelse(survdata$Risk>as.numeric(summary(res.cut)[1]),'High risk','Low risk')
    
  }else if(i==2){
    response <- read.table('rawdata/Melanoma/Nathanson_CTLA4/ICB.Nathanson2017_Ipilimumab_Melanoma_Post.clinical',header = T)
    response$group <- ifelse(response$Response=='0','SD/PD','CR/PR')
    colnames(response)[1] <- 'sample'
    Nathanson <- read.table('rawdata/Melanoma/Nathanson_CTLA4/ICB.Nathanson2017_Ipilimumab_Melanoma_Post.self_subtract',header = T,sep = '\t')
    rownames(Nathanson) <- Nathanson$Entrez 
    colnames(Nathanson)[1] <- "ENTREZID"
    df <- bitr( rownames(Nathanson), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
    Nathanson<- merge(Nathanson ,df,by.x='ENTREZID')
    rownames(Nathanson) <- Nathanson$SYMBOL
    exp <- Nathanson[Nathanson$SYMBOL%in%unigenes,-c(1,ncol(Nathanson))]
    exp <- scale(exp) %>% t %>% as.data.frame()
    exp$Risk=predict(mulcox,exp)
    exp$sample <- rownames(exp)
    survdata <- merge(exp,response,by.x='sample')
    res.cut <- surv_cutpoint(survdata,
                             time='OS',
                             event = 'OS.Event',
                             variables = 'Risk')
    
    survdata$Group <- ifelse(survdata$Risk>as.numeric(summary(res.cut)[1]),'High risk','Low risk')
    survdata$OS <- survdata$OS*12
  }else if(i==3){
    response <- read_xlsx('rawdata/Melanoma/PRJEB23709/Gide2019_PD1+CTLA4_Melanoma_RNASeq/mmc3.xlsx')
    response$group <- ifelse(response$`Best RECIST response`%in%c('PD','SD'),'SD/PD',ifelse(response$`Best RECIST response`%in%c('CR','PR'),'CR/PR','NE'))
    response <- response[response$group!='NE',]
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
    Gide<- merge(Gide,df,by.x='ENTREZID')
    rownames(Gide) <- Gide$SYMBOL
    exp <- Gide[Gide$SYMBOL%in%unigenes,-c(1,ncol(Gide))]
    exp <- scale(exp) %>% t %>% as.data.frame()
    exp$Risk=0.7403*exp$COL1A2+0.4458*exp$F2+0.7294*exp$PLAUR-1.3588*exp$P2RX1
    exp$sample <- rownames(exp)
    survdata <- merge(exp,response,by.x='sample')
    res.cut <- surv_cutpoint(survdata,
                             time='PFS',
                             event = 'PFS.E',
                             variables = 'Risk')
    
    survdata$Group <- ifelse(survdata$Risk>as.numeric(summary(res.cut)[1]),'High risk','Low risk')
  }else if(i==4){
    response <- read_xlsx('rawdata/Melanoma/PRJEB23709/Gide2019_PD1_Melanoma_RNASeq/mmc2.xlsx')
    response$group <- ifelse(response$`Best RECIST response`%in%c('SD','PD'),'SD/PD',ifelse(response$`Best RECIST response`%in%c('CR','PR'),'CR/PR','NE'))
    response <- response[response$group!='NE',]
    colnames(response)[1] <- 'sample'
    response$PFS <- round(response$`Progression Free Survival (Days)`/30,2)
    response$OS <- round(response$`Overall Survival (Days)`/30,2)
    response$OS.E <- ifelse(response$`Last Followup Status`=='Alive',0,1)
    response$PFS.E <- ifelse(response$`Last Followup Status`=='Alive',0,1)
    response$PFS.E[response$PFS<response$OS] <- 1
    response$sample <-paste0('X',response$sample) 
    Gide <- read.table('rawdata/Melanoma/PRJEB23709/Gide2019_PD1_Melanoma_RNASeq/ICB.Gide2019_Pembrolizumab-Nivolumab_Melanoma.self_subtract',header = T,row.names = 1)
    Gide$ENTREZID <- rownames(Gide)
    df <- bitr( rownames(Gide), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
    Gide<- merge(Gide ,df,by.x='ENTREZID')
    rownames(Gide) <- Gide$SYMBOL
    exp <- Gide[Gide$SYMBOL%in%unigenes,-c(1,ncol(Gide))]
    exp <- scale(exp) %>% t %>% as.data.frame()
    exp$Risk=0.7403*exp$COL1A2+0.7294*exp$PLAUR-1.3588*exp$P2RX1
    exp$sample <- rownames(exp)
    survdata <- merge(exp,response,by.x='sample')
    res.cut <- surv_cutpoint(survdata,
                             time='PFS',
                             event = 'PFS.E',
                             variables = 'Risk')
    
    survdata$Group <- ifelse(survdata$Risk>as.numeric(summary(res.cut)[1]),'High risk','Low risk')
  }else if(i==5){
    gset <- getGEO('GSE93157', GSEMatrix =TRUE, AnnotGPL=FALSE)
    pd <- pData(gset[[1]])
    exp <- exprs(gset[[1]])%>%as.data.frame()
    pd <- pd[pd$source_name_ch1%in%c('LUNG NON-SQUAMOUS CANCER'),c(2,48,55,56)]
    colnames(pd) <- c('sample','response','PFS','PFS.E')
    exp <- exp%>%t%>%as.data.frame()
    exp <- exp[,colnames(exp)%in%unigenes]
    exp$sample <- rownames(exp)
    clin <- merge(exp,pd,by.x='sample')
    clin$group <- ifelse(clin$response%in%c('PD','SD'),'SD/PD','CR/PR')
    clin$Risk <- -0.1753*clin$C4BPA+0.7294*clin$PLAUR
    clin$PFS <- as.numeric(clin$PFS)
    clin$PFS.E <- as.numeric(clin$PFS.E)
    
    res.cut <- surv_cutpoint(clin,
                             time='PFS',
                             event = 'PFS.E',
                             variables = 'Risk')
    
    clin$Group <- ifelse(clin$Risk>as.numeric(summary(res.cut)[1]),'High risk','Low risk')
    
  }else if(i==6){
    urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
    path <- paste(urld, "acc=GSE161537", "file=GSE161537_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
    tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
    dat <- log10(tbl + 1)%>%as.data.frame()
    df <-bitr(rownames(dat), fromType = "ENTREZID", toType = c( "SYMBOL" ), OrgDb = org.Hs.eg.db )
    dat$ENTREZID <- rownames(dat)
    exp <- merge(dat,df,by.x='ENTREZID')
    exp<-exp[!duplicated(exp$SYMBOL),]
    rownames(exp)<-exp$SYMBOL[!duplicated(exp$SYMBOL)]
    exp <- exp[,!colnames(exp)%in%c('ENTREZID','SYMBOL')]
    exp<- exp[rownames(exp)%in%unigenes,]
    exp <- scale(exp)%>%t%>%as.data.frame()
    exp$sample <- rownames(exp)
    gset <- getGEO('GSE161537', GSEMatrix =TRUE, AnnotGPL=FALSE)
    pd <- pData(gset[[1]])
    pd <- pd[,c(2,11,16,17,50)]
    colnames(pd) <- c('sample','type','OS','OS.E','response')
    pd$OS <- substr(pd$OS,13,30)
    pd$OS.E <- substr(pd$OS.E,13,14)
    pd$type <- substr(pd$type,12,14)
    exp$sample <- rownames(exp)
    clin <- merge(exp,pd,by.x='sample')
    clin$Risk=predict(mulcox,clin)
    clin$OS <- as.numeric(clin$OS)
    clin$OS.E <- as.numeric(clin$OS.E)
    res.cut <- surv_cutpoint(clin,
                             time='OS',
                             event = 'OS.E',
                             variables = 'Risk')
    
    clin$Group <- ifelse(clin$Risk>as.numeric(summary(res.cut)[1]),'High risk','Low risk')
    clin$group <- ifelse(clin$response%in%c('PD','SD'),'SD/PD',ifelse(clin$response%in%c('CR','PR'),'CR/PR','NE'))
    clin <- clin[clin$group!='NE',]
    survdata <- clin
  }else if(i==7){
    gset <- getGEO('GSE135222', GSEMatrix =TRUE, AnnotGPL=FALSE)
    pd <- pData(gset[[1]])
    urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
    path <- paste(urld, "acc=GSE135222", "file=GSE135222_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
    tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
    dat <- log10(tbl + 1)
    dat <- as.data.frame(dat)
    
    dat[1:4,1:4]
    gf <- read.csv('files/gf.csv',header=T)
    dat <- dat[gf$ENTREZID[gf$SYMBOL%in%unigenes][-4],]
    rownames(dat) <- gf$SYMBOL[gf$SYMBOL%in%unigenes][-4]
    exp <- dat%>%t%>%as.data.frame()
    exp$Risk <- predict(mulcox,exp)
    exp$sample <- rownames(exp)
    clin <- pd[,c(2,43,44)]
    colnames(clin) <- c('sample','PFS','PFS.E')
    clin <- merge(exp,clin,by.x='sample')
    clin$PFS <- as.numeric(clin$PFS)/30
    clin$PFS.E <- as.numeric(clin$PFS.E)
    res.cut <- surv_cutpoint(clin,
                             time='PFS',
                             event = 'PFS.E',
                             variables = 'Risk')
    
    clin$Group <- ifelse(clin$Risk>as.numeric(summary(res.cut)[1]),'High risk','Low risk')
    survdata <- clin
  }else if(i==8){
    gset <- getGEO('GSE14814', GSEMatrix =TRUE, AnnotGPL=FALSE)
    pd <- pData(gset[[1]])
    exp <- exprs(gset[[1]])%>%as.data.frame()
    gpl <- getGEO(gset[[1]]@annotation)
    gpl <- gpl@dataTable@table
    gpl <- gpl[,c(1,11)]
    exp$ID <- rownames(exp)
    exp <- merge(exp,gpl,by.x='ID')
    exp <- exp[exp$`Gene Symbol`!='',]
    exp <- exp[!duplicated(exp$`Gene Symbol`),]
    rownames(exp) <- exp$`Gene Symbol`
    exp <- exp[,-c(1,135)]
    exp <- scale(exp)%>%t%>%as.data.frame()
    exp <- exp[,colnames(exp)%in%unigenes]
    exp$sample <- rownames(exp)
    clin <- pd[pd$`Histology type:ch1`%in%c('ADC', 'SQCC' ),c(2,53:56)]
    clin$OS <- as.numeric(clin$`OS time:ch1`)*12
    clin$OS.E <- ifelse(clin$`OS status:ch1`=="Dead",1,0)
    clin <- clin[,c(1,5,6,7)]
    colnames(clin)[1:2] <- c('sample','type')
    clin <- merge(clin,exp,by.x='sample')
    clin$Risk <- predict(mulcox,clin)
    res.cut <- surv_cutpoint(clin,
                             time='OS',
                             event = 'OS.E',
                             variables = 'Risk')
    
    clin$Group <- ifelse(clin$Risk>as.numeric(summary(res.cut)[1]),'High risk','Low risk')
    survdata <- clin
  }
  return(survdata)
}


type=c('ACT','CTLA4','PD1_CTLA4','PD1','NSCLC_PD1')


clin <- immunotherapy(1)
png("figures/06E_OSKM_ACT.png",width = 1500,height = 1500,res=200)
p1 <- ggsurvplot(survfit(Surv(PFS, PFS.Event) ~ Group, data=clin),size=3,
                 surv.median.line = "hv",
                 title = "GSE100797_ACT",
                 legend.labs = c('High-risk','Low-risk'),
                 legend=c(0.8,0.96),font.legend=c(20,'bold'),
                 palette = c("#E69F00", "#56B4E9",'grey'),
                 pval = T,pval.size=18,legend.title = "",
                 ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                 font.main = c(30, "bold"),
                 ylab='Survival probability of PFS',xlab='Time(months)',
                 font.x = c(30, "bold"),
                 font.y = c(30, "bold"))
dev.off()
p5 <- ggboxplot(clin,x="group",y="Risk",fill = "group")+theme_classic()+ylab('CRRS')+xlab('')+ggtitle('GSE100797_ACT')+ylim(-6,5)+
  geom_signif(size=1,textsize=12,comparisons = list(c('CR/PR','SD/PD')),test = wilcox.test,map_signif_level = T)+
  theme( plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
         axis.text.x = element_text(size=18,face='bold',color='black'),
         axis.text.y = element_text(size=18,face='bold',color='black'),
         text=element_text(size=25,face='bold'),
         legend.title = element_blank(),legend.position = "none",
         axis.line = element_line(size=1))
ggsave('figures/06_ACT_boxplot.png',p5,width = 5,height = 5,dpi=200)


clin <- immunotherapy(2)
png("figures/06E_OSKM_CTLA4.png",width = 1500,height = 1500,res=200)
p2 <- ggsurvplot(survfit(Surv(OS,OS.Event) ~ Group, data=clin),size=3,
                 surv.median.line = "hv",
                 title = "Nathanson_CTLA4",
                 legend.labs =  c('High-risk','Low-risk'),
                 legend=c(0.88,1),font.legend=c(20,'bold'),
                 palette = c("#E69F00", "#56B4E9",'grey'),
                 pval = T,pval.size=18,legend.title = "",
                 ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                 font.main = c(30, "bold"),
                 ylab='Survival probability of OS',xlab='Time(months)',
                 font.x = c(30, "bold"),
                 font.y = c(30, "bold"))
dev.off()
p6 <- ggboxplot(clin,x="group",y="Risk",fill = "group")+theme_classic()+ylab('CRRS')+xlab('')+ggtitle('Nathanson_CTLA4')+ylim(-7,7)+
  geom_signif(size=1,textsize=12,comparisons = list(c('CR/PR','SD/PD')),test = wilcox.test,map_signif_level = T)+
  theme( plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
         axis.text.x = element_text(size=18,face='bold',color='black'),
         axis.text.y = element_text(size=18,face='bold',color='black'),
         text=element_text(size=25,face='bold'),
         legend.title = element_blank(),legend.position = "none",
         axis.line = element_line(size=1))
ggsave('figures/06_CTLA4_boxplot.png',p6,width = 5,height = 5,dpi=200)


clin <- immunotherapy(3)
png("figures/06E_OSKM_PD1_CTLA4.png",width = 1500,height = 1500,res=200)
p3 <- ggsurvplot(survfit(Surv(PFS, PFS.E) ~ Group, data=clin),size=3,
                 surv.median.line = "hv",
                 title = "PRJEB23709_PD1_CTLA4",
                 legend.labs =  c('High-risk','Low-risk'),
                 legend=c(0.8,0.96),font.legend=c(20,'bold'),
                 palette = c("#E69F00", "#56B4E9",'grey'),
                 pval = T,pval.size=18,legend.title = "",
                 ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                 font.main = c(30, "bold"),
                 ylab='Survival probability of PFS',xlab='Time(months)',
                 font.x = c(30, "bold"),
                 font.y = c(30, "bold"))
dev.off()
p7 <- ggboxplot(clin,x="group",y="Risk",fill = "group")+theme_classic()+ylab('CRRS')+xlab('')+ggtitle('PRJEB23709_PD1_CTLA4')+ylim(-6,5)+
  geom_signif(size=1,textsize=12,comparisons = list(c('CR/PR','SD/PD')),test = wilcox.test,map_signif_level = T)+
  theme( plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
         axis.text.x = element_text(size=18,face='bold',color='black'),
         axis.text.y = element_text(size=18,face='bold',color='black'),
         text=element_text(size=25,face='bold'),
         legend.title = element_blank(),legend.position = "none",
         axis.line = element_line(size=1))
ggsave('figures/06_PD1_CTLA4_boxplot.png',p7,width = 5,height = 5,dpi=200)


clin <- immunotherapy(4)
png("figures/06E_OSKM_PD1.png",width = 1500,height = 1500,res=200)
p4 <- ggsurvplot(survfit(Surv(PFS, PFS.E) ~ Group, data=clin),size=3,
                 surv.median.line = "hv",
                 title = "PRJEB23709_PD1",
                 legend.labs =  c('High-risk','Low-risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                 palette = c("#E69F00", "#56B4E9",'grey'),
                 pval = T,pval.size=18,legend.title = "",
                 ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                 font.main = c(30, "bold"),
                 ylab='Survival probability of PFS',xlab='Time(months)',
                 font.x = c(30, "bold"),
                 font.y = c(30, "bold"))
dev.off()
p8 <- ggboxplot(clin,x="group",y="Risk",fill = "group")+theme_classic()+ylab('CRRS')+xlab('')+ggtitle('PRJEB23709_PD1')+ylim(-6,5)+
  geom_signif(size=1,textsize=12,comparisons = list(c('CR/PR','SD/PD')),test = wilcox.test,map_signif_level = T)+
  theme( plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
         axis.text.x = element_text(size=18,face='bold',color='black'),
         axis.text.y = element_text(size=18,face='bold',color='black'),
         text=element_text(size=25,face='bold'),
         legend.title = element_blank(),legend.position = "none",
         axis.line = element_line(size=1))
ggsave('figures/06_PD1_boxplot.png',p8,width = 5,height = 5,dpi=200)


clin <- immunotherapy(5)
png("figures/06E_OSKM_NSCLC_PD1_GSE93157.png",width = 1500,height = 1500,res=200)
p9 <- ggsurvplot(survfit(Surv(PFS, PFS.E) ~ Group, data=clin),size=3,
                 surv.median.line = "hv",
                 title = "anti-PD-1_GSE93157",
                 legend.labs =  c('High-risk','Low-risk'),legend=c(0.8,0.96),font.legend=c(20,'bold'),
                 palette = c("#E69F00", "#56B4E9",'grey'),
                 pval = T,pval.size=18,legend.title = "",
                 ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                 font.main = c(30, "bold"),
                 ylab='Overall survival',xlab='Time(months)',
                 font.x = c(30, "bold"),
                 font.y = c(30, "bold"))
dev.off()
p10 <- ggboxplot(clin,x="group",y="Risk",fill = "group")+theme_classic()+ylab('CRRS')+xlab('')+ggtitle('anti-PD-1_GSE93157')+ylim(3,9)+
  geom_signif(size=1,textsize=12,comparisons = list(c('CR/PR','SD/PD')),test = wilcox.test,map_signif_level = T)+
  theme( plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
         axis.text.x = element_text(size=18,face='bold',color='black'),
         axis.text.y = element_text(size=18,face='bold',color='black'),
         text=element_text(size=25,face='bold'),
         legend.title = element_blank(),legend.position = "none",
         axis.line = element_line(size=1))

ggsave('figures/06_NSCLC_PD1_boxplot.png',p10,width = 5,height = 5,dpi=200)

library(pROC)
clin <- immunotherapy(5)
roc1 <- roc(clin$group,clin$Risk)
png("figures/06_roc_NSCLC.png",width = 1500,height = 1500,res=200)
plot(roc1,col="red",lwd=4,cex.axis=1.5,ann = F, xlim=c(1,0), ylim=c(0,1))
title(xlab="Specificity", ylab="Sensitivity",cex.lab=1.6,font.lab=2)
legend(0.55,0.2,c(paste('AUC =',round(roc1$auc,3)) ),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()

clin <- immunotherapy(6)
png("figures/06_PD-1_PD-L1_GSE161537.png",width = 1500,height = 1500,res=200)
p11 <- ggsurvplot(survfit(Surv(OS, OS.E) ~ Group, data=clin),size=3,
                 surv.median.line = "hv",
                 title = "anti-PD-1/PD-L1(GSE161537)",
                 legend.labs =  c('High-risk','Low-risk'),
                 legend=c(0.8,0.96),font.legend=c(20,'bold'),
                 palette = c("#E69F00", "#56B4E9",'grey'),
                 pval = T,pval.size=18,legend.title = "",
                 ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                 font.main = c(30, "bold"),
                 ylab='Overall survival',xlab='Time(months)',
                 font.x = c(30, "bold"),
                 font.y = c(30, "bold"));p11
dev.off()


p12 <- ggboxplot(clin,x="group",y="Risk",fill = "group")+theme_classic()+ylab('CRRS')+xlab('')+ggtitle('LUAD_PD1/PD-L1')+ylim(0,2)+
  geom_signif(size=1,textsize=12,comparisons = list(c('CR/PR','SD/PD')),test = wilcox.test,map_signif_level = T)+
  theme( plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
         axis.text.x = element_text(size=18,face='bold',color='black'),
         axis.text.y = element_text(size=18,face='bold',color='black'),
         text=element_text(size=25,face='bold'),
         legend.title = element_blank(),legend.position = "none",
         axis.line = element_line(size=1))

ggsave('figures/06_NSCLC_PD1_PD-L1_boxplot.png',p12,width = 5,height = 5,dpi=200)

library(pROC)
roc2 <- roc(clin$group,clin$Risk)
png("figures/06_roc_NSCLC_PD1_PD-L1.png",width = 1500,height = 1500,res=200)
plot(roc2,col="red",lwd=4,cex.axis=1.5,ann = F, xlim=c(1,0), ylim=c(0,1))
title(xlab="Specificity", ylab="Sensitivity",cex.lab=1.6,font.lab=2)
legend(0.55,0.2,c(paste('AUC =',round(roc2$auc,3)) ),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=4,col = c('red','orange','blue',"green","purple",'black',"purple3","yellow","pink","cyan","magenta","seagreen4","red4"),
       bty='n',seg.len = 1.5,cex = 1.2)

dev.off()


clin <- immunotherapy(7)
png("figures/06_PD-1_PD-L1_GSE135222.png",width = 1500,height = 1500,res=200)
p13 <- ggsurvplot(survfit(Surv(PFS, PFS.E) ~ Group, data=clin),size=3,
                  surv.median.line = "hv",
                  title = "anti-PD-1/PD-L1(GSE135222)",
                  legend.labs =  c('High-risk','Low-risk'),
                  legend=c(0.8,0.96),font.legend=c(20,'bold'),
                  palette = c("#E69F00", "#56B4E9",'grey'),
                  pval = T,pval.size=18,legend.title = "",
                  ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                  font.main = c(30, "bold"),
                  ylab='Progression-free survival',xlab='Time(months)',
                  font.x = c(30, "bold"),
                  font.y = c(30, "bold"));P13
dev.off()

clin <- immunotherapy(8)
png("figures/06_ACT_GSE14814.png",width = 1500,height = 1500,res=200)
p14 <- ggsurvplot(survfit(Surv(OS, OS.E) ~ Group, data=clin[clin$type=='ACT',]),size=3,
                  surv.median.line = "hv",
                  title = "ACT(GSE14814)",
                  #legend.labs =  c('High-risk','Low-risk'),
                  legend=c(0.8,0.96),font.legend=c(20,'bold'),
                  palette = c("#E69F00", "#56B4E9",'grey'),
                  pval = T,pval.size=18,legend.title = "",
                  ggtheme = theme_classic()+theme(plot.title=element_text(hjust = 0.5,face = 'bold')), 
                  font.main = c(30, "bold"),
                  ylab='Progression-free survival',xlab='Time(months)',
                  font.x = c(30, "bold"),
                  font.y = c(30, "bold"));P13
dev.off()


