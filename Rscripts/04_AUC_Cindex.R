##########################################LOAD DATA
clinf <- function(n){
  if(n==1){
    exp <- lxdata1[,c('sampleID',gene)]
    colnames(exp)[1] <- 'sample'
    clin <- merge(exp,clin_TCGA)
    clin$Risk <- predict(mulcox,clin)
  }else if(n==2){
    exp <- GSE13213%>%t%>%as.data.frame()
    exp$sample <- colnames(GSE13213)
    clin <- merge(exp,clin117)
    clin$Risk <- predict(mulcox,clin)
  }else if(n==3){
    exp <- GSE31210%>%t%>%as.data.frame()
    exp$sample <- colnames(GSE31210)
    clin <- merge(exp,clin226)
    clin$Risk <- predict(mulcox,clin)
  }else if(n==4){
    exp <- GSE70294%>%t%>%as.data.frame()
    exp$sample <- colnames(GSE70294)
    clin <- merge(exp,clin398)
    clin$Risk <- predict(mulcox,clin)
  }else if(n==5){
    exp <- GSE30219%>%t%>%as.data.frame()
    exp$sample <- colnames(GSE30219)
    clin <- merge(exp,clin293)
    clin$Risk <- predict(mulcox,clin)
  }else if(n==6){
    exp <- GSE68465%>%t%>%as.data.frame()
    exp$sample <- colnames(GSE68465)
    clin <- merge(exp,clin433)
    clin$Risk <- predict(mulcox,clin)
  }else if(n==7){
    exp <- GSE3141%>%t%>%as.data.frame()
    exp$sample <- colnames(GSE3141)
    clin <- merge(exp,clin58)
    clin$Risk <- predict(mulcox,clin)
  }else if(n==8){
    exp <- expr_batch%>%t%>%as.data.frame()
    exp$sample <- colnames(expr_batch)
    clin <- merge(exp,clin2027)
    clin$Risk <- predict(mulcox,clin)
  }
  return(clin)
}


####################################################cindex######################################
Cindex <- c()
std <- c()
for(i in 1:8){
  clin <- clinf(i)
  cox <- coxph(Surv(Time, Status) ~ Risk, data = clin,x=T,y=T)
  Cindex[i] <- cox$concordance[6]
  std[i] <- cox$concordance[7]
}


cohort <- c('TCGA-LUAD',
            'GSE13213',
            'GSE31210',
            'GSE70294',
            'GSE30219',
            'GSE68465',
            'GSE3141',
            'Meta-Cohort')

cdf <- data.frame(cohort=cohort,Cindex=Cindex,std=std)
discrete_palettes <- rep('skyblue',8)
p<- ggplot(cdf,aes(x=cohort,y=Cindex,fill = cohort))+geom_bar(stat = "identity")+coord_flip()+theme_bw()+ylab('C-index')+xlab('')+
  geom_errorbar(aes(ymin=Cindex-std,ymax=Cindex+std),width=0.2,size=1)+
  theme(legend.position = 'none',axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),text=element_text(size=18,face='bold'),panel.grid = element_blank())+scale_fill_manual(values = discrete_palettes)

ggsave('figures/03A_cindex.png',p,width = 10,height = 10,dpi=200)
library(xlsx)
write.xlsx(cdf,'figures/03_data.xlsx',row.names = F)
####################################################cindex-clinical######################################
clins <- function(n,c){
  if(n==1){
    TMB <- read.xlsx('files/02_tmb_msi.xlsx',sheetName = 'TMB')
    TMB <- TMB[,1:2]
    MSI <- read.xlsx('files/02_tmb_msi.xlsx',sheetName = 'MSI')
    MSI <- MSI[,1:2]
    exp <- lxdata1[,c('sampleID',gene,'Age','Gender','Tumor','Lymph_Node','Metastasis','Stage')]
    colnames(exp) <- c('sample',gene,'Age','Gender','pT','pN','pM','Stage')
    clin <- merge(exp,clin_TCGA)
    clin$CRRS <- predict(mulcox,clin)
    mulcox1<- coxph(Surv(Time, Status) ~ CRRS+Stage, data = clin)
    clin$CRRS_Stage <- predict(mulcox1,clin)
    clin <- merge(clin,TMB,by.x='sample')
    clin <- merge(clin,MSI,by.x='sample')
    if(c=='CRRS'){
      clin <- clin[,c('sample','Time','Status','CRRS','Age','Gender','pT','pN','pM','Stage','TMB','MSI')]
    }else if(c=='CRRS_Stage'){
      clin <- clin[,c('sample','Time','Status','CRRS_Stage','Age','Gender','pT','pN','pM','Stage','TMB','MSI')]
    }
    
  }else if(n==2){
    exp <- GSE13213%>%t%>%as.data.frame()
    exp$sample <- colnames(GSE13213)
    clin <- merge(exp,clin117)
    clin$CRRS  <- predict(mulcox,clin)
    clin <- clin[,c('sample','Time','Status','CRRS')]
    clin.tmp <- clindata117[,c('sample','age','gender','stage')]
    clin.tmp$age <- ifelse(clin.tmp$age>70 ,'old', 'young')
    clin <- merge(clin,clin.tmp)
    colnames(clin) <- c('sample','Time','Status','CRRS','Age','Gender','Stage')
    if(c=='CRRS_Stage'){
      mulcox1<- coxph(Surv(Time, Status) ~ CRRS+Stage, data = clin)
      clin$CRRS_Stage <- predict(mulcox1,clin)
      clin <- clin[,c('sample','Time','Status','CRRS_Stage','Age','Gender','Stage')]
    }
    
  }else if(n==3){
    exp <- GSE31210%>%t%>%as.data.frame()
    exp$sample <- colnames(GSE31210)
    clin <- merge(exp,clin226)
    clin$CRRS  <- predict(mulcox,clin)
    clin <- clin[,c('sample','Time','Status','CRRS')]
    clin.tmp <- clindata226[,c('sample','age','gender','stage')]
    clin.tmp$age <- ifelse(clin.tmp$age>70 ,'old', 'young')
    clin <- merge(clin,clin.tmp)
    colnames(clin) <- c('sample','Time','Status','CRRS','Age','Gender','Stage')
    if(c=='CRRS_Stage'){
      mulcox1<- coxph(Surv(Time, Status) ~ CRRS+Stage, data = clin);mulcox1
      clin$CRRS_Stage <- predict(mulcox1,clin)
      clin <- clin[,c('sample','Time','Status','CRRS_Stage','Age','Gender','Stage')]
    }
    
  }else if(n==4){
    exp <- GSE70294%>%t%>%as.data.frame()
    exp$sample <- colnames(GSE70294)
    clin <- merge(exp,clin398)
    clin$CRRS  <- predict(mulcox,clin)
    clin <- clin[,c('sample','Time','Status','CRRS')]
    clin.tmp <- clindata398[,c('sample','age','gender','stage')]
    clin.tmp$age <- ifelse(clin.tmp$age>70 ,'old', 'young')
    clin <- merge(clin,clin.tmp)
    colnames(clin) <- c('sample','Time','Status','CRRS','Age','Gender','Stage')
    if(c=='CRRS_Stage'){
      mulcox1<- coxph(Surv(Time, Status) ~ CRRS+Stage, data = clin);mulcox1
      clin$CRRS_Stage <- predict(mulcox1,clin)
      clin <- clin[,c('sample','Time','Status','CRRS_Stage','Age','Gender','Stage')]
    }
    
  }else if(n==5){
    exp <- GSE30219%>%t%>%as.data.frame()
    exp$sample <- colnames(GSE30219)
    clin <- merge(exp,clin293)
    clin$CRRS  <- predict(mulcox,clin)
    clin <- clin[,c('sample','Time','Status','CRRS')]
    clin.tmp <- clindata293[,c('sample','age','gender','pT','pN','pM')]
    clin.tmp$age <- ifelse(clin.tmp$age>70 ,'old', 'young')
    clin <- merge(clin,clin.tmp)
    colnames(clin) <- c('sample','Time','Status','CRRS','Age','Gender','pT','pN','pM')
    if(c=='CRRS_Stage'){
      mulcox1<- coxph(Surv(Time, Status) ~ CRRS+pT+pN+pM, data = clin);mulcox1
      clin$CRRS_Stage <- predict(mulcox1,clin)
      clin <- clin[,c('sample','Time','Status','CRRS_Stage','Age','Gender','pT','pN','pM')]
    }
    
  }else if(n==6){
    exp <- GSE68465%>%t%>%as.data.frame()
    exp$sample <- colnames(GSE68465)
    clin <- merge(exp,clin433)
    clin$CRRS  <- predict(mulcox,clin)
    clin <- clin[,c('sample','Time','Status','CRRS')]
    clin.tmp <- clindata433[,c('sample','age','gender','stage')]
    clin.tmp$age <- ifelse(clin.tmp$age>70 , 'old','young')
    clin <- merge(clin,clin.tmp)
    colnames(clin) <- c('sample','Time','Status','CRRS','Age','Gender','Stage')
    if(c=='CRRS_Stage'){
      mulcox1<- coxph(Surv(Time, Status) ~ CRRS+Stage, data = clin);mulcox1
      clin$CRRS_Stage <- predict(mulcox1,clin)
      clin <- clin[,c('sample','Time','Status','CRRS_Stage','Age','Gender','Stage')]
    }
    
  }else if(n==7){
    exp <- GSE3141%>%t%>%as.data.frame()
    exp$sample <- colnames(GSE3141)
    clin <- merge(exp,clin58)
    clin$CRRS  <- predict(mulcox,clin)
    clin <- clin[,c('sample','Time','Status','CRRS')]
    
  }else if(n==8){
    exp <- expr_batch%>%t%>%as.data.frame()
    exp$sample <- colnames(expr_batch)
    clin <- merge(exp,clin2027)
    clin$CRRS  <- predict(mulcox,clin)
    clin <- clin[,c('sample','Time','Status','CRRS')]
    
  }
  return(clin)
}
#clin <- clins(1,'CRRS')
resdf <- c()
par(mfrow=c(3,2))
# c='CRRS'
c='CRRS_Stage'
for(i in 1:6){
  clin <- clins(i,c)
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
  res$cohort=rep(cohort[i],length(res$signatures))
  discrete_palettes <- RColorBrewer::brewer.pal(n=nrow(res), name='Spectral')
  p=ggplot(res,aes(x=signatures,y=Cindex,fill = signatures))+geom_bar(stat = "identity")+
    theme_bw()+ylab('C-index(Compared with CRRS+Stage)')+xlab('')+ggtitle(cohort[i])+
    geom_errorbar(aes(ymin=Cindex-std,ymax=Cindex+std),width=0.2,size=1)+
    theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),legend.position = 'none',axis.text.x = element_text(size=18,face='bold',color='black'),
          axis.text.y = element_text(size=18,face='bold',color='black'),
          text=element_text(size=18,face='bold'),panel.grid = element_blank())+scale_fill_manual(values = discrete_palettes)
  resdf <- rbind(resdf,res)
  
  if(i==1){
    p1=p
    pl=p1
  }
  else if(i==2){
    p2=p
    pl <- p1+p2}
  else if(i==3){
    p3=p
    pl <- p1+p2+p3}
  else if(i==4){
    p4=p
    pl <- p1+p2+p3+p4}
  else if(i==5){
    p5=p
    pl <- p1+p2+p3+p4+p5}
  else if(i==6){
    p6=p
    pl <- p1+p2+p3+p4+p5+p6}

}
pl
ggsave('figures/03A_cindex_CRRS_Stage_compare.png',pl,width = 30,height = 15,dpi=200)

library(xlsx)
write.xlsx(resdf,'figures/03_data.xlsx',row.names = F,sheet='cindex_compare_CRRS+Stage',append = T)

####################################################AUC####################################################
clinc <- function(n){
  if(n==1){
    TMB <- read.xlsx('files/02_tmb_msi.xlsx',sheetName = 'TMB')
    TMB <- TMB[,1:2]
    MSI <- read.xlsx('files/02_tmb_msi.xlsx',sheetName = 'MSI')
    MSI <- MSI[,1:2]
    exp <- lxdata1[,c('sampleID',gene,'Age','Gender','Tumor','Lymph_Node','Metastasis','Stage','RFS','RFS.E')]
    colnames(exp) <- c('sample',gene,'Age','Gender','pT','pN','pM','Stage','RFS','RFS.E')
    clin <- merge(exp,clin_TCGA)
    clin$CRRS <- predict(mulcox,clin)
    clin <- merge(clin,TMB,by.x='sample')
    clin <- merge(clin,MSI,by.x='sample')
    clin <- clin[,c('sample','Time','Status','RFS','RFS.E','CRRS','Age','Gender','pT','pN','pM','Stage','TMB','MSI')]
    clin$T <- ifelse(clin$pT%in%c('T1','T2'),'T1-2','T3-4')
    clin$N <- ifelse(clin$pN%in%c('N0','N1'),'N0-1','N2-3')
    clin$M <- clin$pM
    clin$Stage <- ifelse(clin$Stage%in%c('I','II'),'I-II','III-IV')
    mulcox1<- coxph(Surv(Time, Status) ~ CRRS+Age+Gender+T+N+M+Stage+TMB+MSI, data = clin)
    forest1 <- ggforest(mulcox1,  
                        data = clin,  
                        main = 'TCGA overall survival', 
                        cpositions = c(0.05, 0.15, 0.35), 
                        fontsize = 2, 
                        refLabel = 'Reference', 
                        noDigits = 3 
    )
    mulcox2<- coxph(Surv(RFS,RFS.E) ~ CRRS++Age+Gender+pT+pN+pM+Stage+TMB+MSI, data = clin)
    forest2 <- ggforest(mulcox2,  
                        data = clin,  
                        main = 'TCGA recurrence-free survival', 
                        cpositions = c(0.05, 0.15, 0.35), 
                        fontsize = 2, 
                        refLabel = 'Reference', 
                        noDigits = 3 
    )
    par(mfrow=c(2,1))
    p=forest1+forest2
    p
    ggsave(filename="figures/03S_tcga_muti_forest.png",p, width=40, height=14)
    
  }else if(n==2){
    library(dplyr)
    load('rawdata/FUSCC/omics_clin.Rdata')
    hubgenes <- lasso_genes
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
                        main = 'FUSCC recurrence-free survival', 
                        cpositions = c(0.05, 0.15, 0.35), 
                        fontsize = 2, 
                        refLabel = 'Reference', 
                        noDigits = 3 
    )
    par(mfrow=c(2,1))
    p=forest1+forest2
    p
    ggsave(filename="figures/03S_fuscc_muti_forest.png",p, width=40, height=14)
    
    
  }
  
  

}


####################################################AUC####################################################
auc.1 <- c()
auc.2 <- c()
auc.3 <- c()
for(x in 1:8){
  clin <- clinf(x)
  for(i in c(12,36,60)){
    roc= survivalROC(Stime=clin$Time,  
                     status=clin$Status,   
                     marker = clin$Risk,     
                     predict.time = i, method="KM")
    
    if(i==12){ auc.1 <- c(auc.1,roc$AUC)}else if(i==36){auc.2 <- c(auc.2,roc$AUC)}else if(i==60){auc.3 <- c(auc.3,roc$AUC)}
    
  }
}


adf <- data.frame(cohort=cohort,Year1=auc.1,Year3=auc.2,Year5=auc.3)
adf$Year1 <- round(adf$Year1 ,3)
adf$Year3 <- round(adf$Year3 ,3)
adf$Year5 <- round(adf$Year5 ,3)

library(xlsx)
write.xlsx(adf,'figures/03_data.xlsx',row.names = F,sheetName = 'AUC',append = T)
discrete_palettes1 <- rep(RColorBrewer::brewer.pal(n=10, name='Spectral')[2],8)
discrete_palettes2 <- rep(RColorBrewer::brewer.pal(n=10, name='Spectral')[6],8)
discrete_palettes3 <- rep(RColorBrewer::brewer.pal(n=10, name='Spectral')[8],8)
par(mfrow=c(3,1))
g1<-ggplot(adf,aes(x=cohort,y=Year1,fill = cohort))+geom_bar(stat = "identity")+coord_flip()+theme_bw()+ylab('')+xlab('')+ggtitle('1-Year')+
  geom_text(aes(label = Year1), size = 5, hjust = 1, vjust = 0, position =     "stack")+ 
  theme(plot.title=element_text(hjust = 0.5,face = 'bold'), legend.position = 'none',axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),text=element_text(size=18,face='bold'),panel.grid = element_blank())+scale_fill_manual(values = discrete_palettes1)
g2<-ggplot(adf,aes(x=cohort,y=Year3,fill = cohort))+geom_bar(stat = "identity")+coord_flip()+theme_bw()+ylab('AUC')+xlab('')+ggtitle('3-Year')+
  geom_text(aes(label = Year3), size = 5, hjust = 1, vjust = 0, position =     "stack")+ 
  theme(plot.title=element_text(hjust = 0.5,face = 'bold'),legend.position = 'none',axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_blank(),text=element_text(size=18,face='bold'),panel.grid = element_blank())+scale_fill_manual(values = discrete_palettes2)
g3<-ggplot(adf,aes(x=cohort,y=Year5,fill = cohort))+geom_bar(stat = "identity")+coord_flip()+theme_bw()+ylab('')+xlab('')+ggtitle('5-Year')+
  geom_text(aes(label = Year5), size = 5, hjust = 1, vjust = 0, position =     "stack")+ 
  theme(plot.title=element_text(hjust = 0.5,face = 'bold'),legend.position = 'none',axis.text.x = element_text(size=18,face='bold',color='black'),
        axis.text.y = element_blank(),text=element_text(size=18,face='bold'),panel.grid = element_blank())+scale_fill_manual(values = discrete_palettes3)
g=g1+g2+g3


# define your own tags
ggsave('figures/03A_auc.png',g,width = 20,height = 10,dpi=200)




