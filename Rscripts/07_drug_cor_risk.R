library(readxl) 
library(impute)
library(limma)
rt <- read.table( "LUADdata/Drug/drug.txt",sep= "\t",header=T,check.names=F)
rt <- as.matrix(rt) 
rownames(rt) <- rt[,1] 
drug <- rt[-nrow(rt), 2:ncol(rt)] 

dimnames <- list(rownames(drug),colnames(drug)) 
data <- matrix(as.numeric(as.matrix(drug)),nrow=nrow(drug),dimnames=dimnames)
mat<- impute.knn(data) #通过impute.knn函数来评估并补齐药物数据NA
drug<- mat$data 
drug<- avereps(drug)


exp<- read.table( "LUADdata/Drug/geneExp.txt", sep= "\t", header=T, row.names = 1, check.names=F) 
# dim( exp) 
# exp[ 1: 4, 1: 4]
genelist <- c('F2','P2RX1','GUCY1A1','PLAUR','PRKCZ','F12','SERPINE1','COL1A2','C4BPA','JMJD7-PLA2G4B','CR2')
genelist <- intersect(genelist,row.names( exp))
exp<- exp[genelist,]
exp1 <- exp %>% t %>% as.data.frame()
exp1 <- do.call(data.frame,exp1)
rownames(exp1 ) <- colnames(exp)
colnames(exp1) <- gsub('\\.','_',colnames(exp1))
exp1$Risk <- 0.09*exp1$F2-0.299*exp1$P2RX1+0.107*exp1$PLAUR-0.124*exp1$PRKCZ+0.093*exp1$F12+0.056*exp1$SERPINE1+
  0.137*exp1$COL1A2-0.022*exp1$C4BPA-0.031*exp1$JMJD7_PLA2G4B-0.017*exp1$CR2
exp2 <- exp1$Risk%>%t%>%as.data.frame()
colnames(exp2) <- rownames(exp1)
rownames(exp2) <- 'Risk'
outTab <- data.frame()
for(m in 1:nrow(exp2)){
  x <- as.numeric(exp2[m,]) 
  for(i in 1:nrow(drug)){ 
    y <- as.numeric(drug[i,]) 
    corT <- cor.test(x,y,method= "pearson") 
    cor <- corT$estimate 
    pvalue <- corT$p.value 
    if( pvalue < 0.05){ 
      outVector <- cbind(rownames(exp2)[m],rownames(drug)[i],cor,pvalue) 
      outTab <- rbind(outTab,outVector)
    } 
    print(i)
  } 
}

outTab <- outTab[order(as.numeric( as.vector(outTab$pvalue))),] 
write.table(outTab, file= "CRG/result/files/RiskdrugCor.txt", sep= "\t", row.names=F, quote=F)



######################################################################################################
outTab1 <- outTab[outTab$pvalue<0.01,]
dd<- outTab1[, 2] 
dd_c1 <- drug[rownames(drug)%in%dd,] %>% t %>% as.data.frame()
dd_c1$sample <- rownames(dd_c1)
dd_c2 <- exp2 %>% t %>% as.data.frame()
dd_c2$sample <- rownames(dd_c2)
dd_c2$group <- ifelse(dd_c2$Risk>median(dd_c2$Risk),'High risk','Low risk')
dd_all <- merge(dd_c1,dd_c2,by.x='sample')

dd1 <- df_p_val1$Celltype[df_p_val1$p.signif!='ns']
dd_all1 <- dd_all[,c(1,which(colnames(dd_all)%in%dd1),ncol(dd_all)-1,ncol(dd_all))]
pd <-dd_all1 %>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(dd_all1)[-c(1,ncol(dd_all1)-1,ncol(dd_all1))]),names_to = 'Celltype',values_to = 'Composition')
library(tidyverse)
library(rstatix)
library(ggtext)
library(ggstatsplot)
library(tidyverse)
library(tidyr)
pd <-dd_all %>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(dd_all)[-c(1,ncol(dd_all)-1,ncol(dd_all))]),names_to = 'Celltype',values_to = 'Composition')
pd $Celltype <- factor(pd $Celltype)
df_p_val1 <- pd  %>% 
  group_by(Celltype) %>% 
  wilcox_test(formula = Composition ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='Celltype')
Gene='Riskscore'
PP<- ggboxplot(pd ,x = 'Celltype', y='Composition',fill='group')+theme_bw()+ggtitle(paste0( "Drug sensitivity of ", Gene))+xlab('')+ylab('IC50')+
  scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
  stat_pvalue_manual(size=14,remove.bracket=T,y.position = 6,vjust=1,df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(angle = 45,hjust = 0.8,vjust = 0.9,size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))

ggsave('CRG/result/figs/07_drug_sensitivity.tiff',PP,width = 16,height = 10,dpi=200)

######################################################################################################

library(ggplot2) 
library(ggpubr)

outTab2 <- outTab1[outTab1$V2%in%dd1,]
outTab2$pval <- ifelse(as.numeric(outTab2$pvalue)< 0.001,"p<0.001" ,paste0('p=',round(as.numeric(outTab2$pvalue),3)))
outTab2$cor <- round(as.numeric(outTab2$cor),3)
outTab2$label <- paste0('cor=',outTab2$cor,' ',outTab2$pval)
outTab2$x=rep(1,nrow(outTab2))
outTab2$y=rep(-2.5,nrow(outTab2))
colnames(outTab2)[2] <- 'Celltype'

cp <- ggplot( data= pd, aes(x = Risk, y = Composition))+ geom_point(size= 1)+ylab("IC50")+xlab( "Riskscore")+theme_bw() +
  stat_smooth(method= "lm",se=FALSE, formula=y~x)+facet_wrap(.~Celltype, ncol=6, labeller = label_wrap_gen(multi_line=TRUE))+
  theme( axis.text.x = element_text(size = 10,face = 'bold',colour = 'black'),
         axis.text.y = element_text(size = 10,face = 'bold',colour = 'black'),
         text=element_text(size=18,face='bold'),
         strip.text.x = element_text(size = 12,face = 'bold',colour = 'black'), 
         legend.position = "none") +geom_text(data=outTab2, mapping=aes(x=x,y=y,label=label),fontface='bold',colour='black',size=5)

ggsave('CRG/result/figs/07_drug_cor.tiff',cp,width = 16,height = 10,dpi=200)

######################################################################################################
chemotherapy <- c('Pemetrexed', 'Docetaxel','Paclitaxel', 'Gemcitabine', 'Cisplatin', 'Carboplatin')
dd<- chemotherapy
dd_c1 <- drug[rownames(drug)%in%dd,] %>% t %>% as.data.frame()
dd_c1$sample <- rownames(dd_c1)
dd_c2 <- exp2 %>% t %>% as.data.frame()
dd_c2$sample <- rownames(dd_c2)
dd_c2$group <- ifelse(dd_c2$Risk>median(dd_c2$Risk),'High risk','Low risk')
dd_all <- merge(dd_c1,dd_c2,by.x='sample')

dd1 <- df_p_val1$Celltype[df_p_val1$p.signif!='ns']
dd_all1 <- dd_all[,c(1,which(colnames(dd_all)%in%dd1),ncol(dd_all)-1,ncol(dd_all))]
pd <-dd_all1 %>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(dd_all1)[-c(1,ncol(dd_all1)-1,ncol(dd_all1))]),names_to = 'Celltype',values_to = 'Composition')
library(tidyverse)
library(rstatix)
library(ggtext)
library(ggstatsplot)
library(tidyverse)
library(tidyr)
pd <-dd_all %>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(dd_all)[-c(1,ncol(dd_all)-1,ncol(dd_all))]),names_to = 'Celltype',values_to = 'Composition')
pd $Celltype <- factor(pd $Celltype)
df_p_val1 <- pd  %>% 
  group_by(Celltype) %>% 
  wilcox_test(formula = Composition ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='Celltype')
Gene='Riskscore'
PP<- ggboxplot(pd ,x = 'Celltype', y='Composition',fill='group')+theme_bw()+ggtitle(paste0( "Chemotherapy drugs sensitivity of ", Gene))+xlab('')+ylab('IC50')+
  scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
  stat_pvalue_manual(size=14,remove.bracket=T,y.position = 6,vjust=1,df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
        axis.text.x = element_text(angle = 45,hjust = 0.8,vjust = 0.9,size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))

ggsave('CRG/result/figs/07_chemotherapydrug_sensitivity.tiff',PP,width = 16,height = 10,dpi=200)






