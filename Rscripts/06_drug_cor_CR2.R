library(readxl) 
rt1 <- read_xlsx("LUADdata/Drug/DTP_NCI60_ZSCORE.xlsx")
table(rt1$ `FDA status`)
# Clinical trial   FDA approved 
#  574            218 
rt1 <- rt1[rt1$`FDA status` %in% c( "FDA approved", "Clinical trial"),] 
rt1 <- rt1[,-c(2:5)] 
write.table(rt1, file = "LUADdata/Drug/drug.txt",sep = "\t",row.names = F,quote = F)

rt2<- read_xls(path = "LUADdata/Drug/RNA__RNA_seq_composite_expression.xls") 
rt2 <- rt2[,-c( 2:6)]
write.table(rt2, file = "LUADdata/Drug/geneExp.txt",sep = "\t",row.names = F,quote = F)
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
genelist <- 'CR2'
genelist <- intersect(genelist,row.names( exp))
exp<- exp[genelist,]
outTab <- data.frame()
for(m in 1:nrow(exp)){
  x <- as.numeric(exp[m,]) 
  for(i in 1:nrow(drug)){ 
    y <- as.numeric(drug[i,]) 
    corT <- cor.test(x,y,method= "pearson") 
    cor <- corT$estimate 
    pvalue <- corT$p.value 
    if( pvalue < 0.05){ 
      outVector <- cbind(rownames(exp)[m],rownames(drug)[i],cor,pvalue) 
      outTab <- rbind(outTab,outVector)
    } 
    print(i)
  } 
}

outTab <- outTab[order(as.numeric( as.vector(outTab$pvalue))),] 
write.table(outTab, file= "CRG/result/files/CR2drugCor.txt", sep= "\t", row.names=F, quote=F)




dd<- outTab[, 2] 
dd_c1 <- drug[rownames(drug)%in%dd,] %>% t %>% as.data.frame()
dd_c1$sample <- rownames(dd_c1)
dd_c2 <- exp %>% t %>% as.data.frame()
dd_c2$sample <- rownames(dd_c2)
dd_c2$group <- ifelse(dd_c2$CR2>median(dd_c2$CR2),'High','Low')
dd_all <- merge(dd_c1,dd_c2,by.x='sample')

dd1 <- df_p_val1$Celltype[df_p_val1$p.signif!='ns']
dd_all1 <- dd_all[,c(1,which(colnames(dd_all)%in%dd1),ncol(dd_all)-1,ncol(dd_all))]
pd <-dd_all1 %>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(dd_all1)[-c(1,ncol(dd_all1)-1,ncol(dd_all1))]),names_to = 'Celltype',values_to = 'Composition')

library(tidyr)
pd <-dd_all %>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(dd_all)[-c(1,ncol(dd_all)-1,ncol(dd_all))]),names_to = 'Celltype',values_to = 'Composition')
pd$Celltype <- factor(pd $Celltype)
df_p_val1 <- pd  %>% 
  group_by(Celltype) %>% 
  wilcox_test(formula = Composition ~ group) %>% 
  add_significance(p.col = 'p',cutpoints = c(0,0.001,0.01,0.05,1),symbols = c('***','**','*','ns')) %>% 
  add_xy_position(x='Celltype')
Gene='CR2'
PP<- ggboxplot(pd ,x = 'Celltype', y='Composition',fill='group')+theme_bw()+ggtitle(paste0( "Drug sensitivity of ", Gene))+xlab('')+ylab('IC50')+
  scale_fill_manual(values=c("#00AFBB", "#E7B800"))+
  stat_pvalue_manual(size=14,remove.bracket=T,y.position = 6,vjust=1,df_p_val1,label = '{p.signif}',tip.length = 0)+
  theme(plot.title = element_text(hjust = 0.5,vjust = 1.5,face = 'bold',size = 25),
    axis.text.x = element_text(angle = 45,hjust = 0.8,vjust = 0.9,size=18,face='bold',color='black'),
        axis.text.y = element_text(size=18,face='bold',color='black'),
        text=element_text(size=25,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))

ggsave('CRG/result/figs/02_drug_sensitivity.tiff',PP,width = 8,height = 8,dpi=200)


library(ggplot2) 
library(ggpubr)

outTab1 <- outTab[outTab$V2%in%dd1,]
outTab1$pval <- ifelse(as.numeric(outTab1$pvalue)< 0.001,"p<0.001" ,paste0('p=',round(as.numeric(outTab1$pvalue),3)))
outTab1$cor <- round(as.numeric(outTab1$cor),3)
outTab1$label <- paste0('cor=',outTab1$cor,'\n',outTab1$pval)
outTab1$x=rep(1,nrow(outTab1))
outTab1$y=rep(-2,nrow(outTab1))
colnames(outTab1)[2] <- 'Celltype'

cp <- ggplot( data= pd, aes(x = CR2, y = Composition))+ geom_point(size= 1)+ylab("IC50")+xlab( "Expression of CR2")+theme_bw() +
  stat_smooth(method= "lm",se=FALSE, formula=y~x)+facet_wrap(.~Celltype, ncol=4, labeller = label_wrap_gen(multi_line=TRUE))+
theme( axis.text.x = element_text(size = 10,face = 'bold',colour = 'black'),
       axis.text.y = element_text(size = 10,face = 'bold',colour = 'black'),
       text=element_text(size=18,face='bold'),
       strip.text.x = element_text(size = 12,face = 'bold',colour = 'black'), 
        legend.position = "none") +geom_text(data=outTab1, mapping=aes(x=x,y=y,label=label),fontface='bold',colour='black',size=5)

ggsave('CRG/result/figs/02_drug_cor.tiff',cp,width = 10,height = 6,dpi=200)





