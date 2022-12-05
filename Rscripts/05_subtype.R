load( "CRG/extdata/featureselect_data.Rdata")
data <- as.matrix(t(exp_CRGs_tumor[,-ncol(exp_CRGs_tumor)]))
result <- dist(t(data), method = "euclidean")
result_hc <- hclust(d = result, method = "ward.D2")
out.id=cutree(result_hc,k=2)
library(FactoMineR)
dat.pca <- PCA(exp_CRGs_tumor[,-ncol(exp_CRGs_tumor)], graph = FALSE)
plot(dat.pca,choix="ind")
group <- as.data.frame(out.id)
colnames(group) <- 'group'
group$Group <- ifelse(group$group=='1','cluster1','cluster2')
library(factoextra)
fviz_pca_ind(dat.pca,
             # 只显示点，不显示文字
             geom.ind = "point", 
             # 用不同颜色表示分组
             col.ind = group$Group, 
             palette = c("#00AFBB", "#E7B800"),
             # 是否圈起来
             addEllipses = TRUE,
             legend.title = "Groups"
)

library(writexl)
groupFile <- cbind(rownames(group),group$Group)
colnames(groupFile ) <- c('sampleID','group')
groupFile <- as.data.frame(groupFile )
write_xlsx(groupFile,path = 'CRG/extdata/GSEA_groupfile.xlsx')
library(readxl)
groupFile <- readxl::read_xlsx('CRG/extdata/GSEA_groupfile.xlsx')
colnames(groupFile ) <- c('sampleID','group')
##############KM
load('CRG/extdata/DEG_TN_TCGA.Rdata')
expFile <- exp[rownames(exp)%in%DEGgenes,colnames(exp)%in%groupFile$sample]
write.table(expFile,'CRG/extdata/GSEA_exp.txt',sep = '\t')
clin_cluster <- as.data.frame(clin_rna[clin_rna$sampleID%in%groupFile$sampleID,])
c <- c()
for(i in clin_cluster$sampleID){
  a <- groupFile$group[which(groupFile$sampleID==i)]
  c<- c(c,a)
}
clin_cluster$group <- c
clin_cluster$Time <- round(clin_cluster$time,2)
clin_cluster <- clin_cluster[clin_cluster$time!=0,]
library(survival)
library(survminer)
cs <- ggsurvplot(survfit(Surv(Time, status) ~ group, data=clin_cluster),surv.median.line = "hv",palette = c('red','blue','grey'),pval = T,pval.size=10)
s <- clin_cluster[,c(8,6,7)]
write_xlsx(s,'CRG/extdata/clin_cluster.xlsx')
##############TMB
Mut_tcga=read.table('LUADdata/TCGA/LUAD_mutation.txt.gz',header = T,fill = TRUE)
Mut_tcga$Hugo_Symbol <- Mut_tcga$gene
Mut_tcga$Chromosome <- Mut_tcga$chr
Mut_tcga$Start_Position <- Mut_tcga$start
Mut_tcga$End_Position <- Mut_tcga$end
Mut_tcga$Reference_Allele <- Mut_tcga$reference
Mut_tcga$Tumor_Seq_Allele2 <- Mut_tcga$alt
Mut_tcga$Variant_Classification<- Mut_tcga$effect
Mut_tcga$Variant_Type<- Mut_tcga$DNA_VAF
Mut_tcga$Tumor_Sample_Barcode<- Mut_tcga$sample
Mut_tcga <- Mut_tcga[Mut_tcga$gene%in%CRGs$gene,]
library(maftools)
laml = read.maf(maf = Mut_tcga)
tmb <- tmb(maf = laml,logScale=F)
tmb <- tmb[tmb$Tumor_Sample_Barcode%in%groupFile$sampleID,]
d <- c()
for(i in tmb$Tumor_Sample_Barcode){
  a <- groupFile$group[which(groupFile$sampleID==i)]
  d<- c(d,a)
}
tmb$group <- d
tmb$total_perMB <- as.numeric(tmb$total_perMB)
library(ggpubr)
library(rstatix)
stat.test <- t.test(tmb$total_perMB ~ tmb$group) 

tb <- ggplot(tmb,aes(group,total_perMB,fill=group))+geom_boxplot()+theme_classic()+ylab('TMB')+xlab('')+
  geom_signif(textsize=5,comparisons = list(c("cluster1", "cluster2")),test = wilcox.test,map_signif_level = F)+
  theme(axis.text.x = element_text(size=16,face='bold',color='black'),
        axis.text.y = element_text(size=16,face='bold',color='black'),
        text=element_text(size=16,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1))
ggsave('CRG/result/figs/16_TMB.tiff',tb,width = 5,height = 5,dpi=1000)

##############Estimate
library(estimate)
library(GSVA)
library(GSEABase)
library(limma)


gene_set<-read_xlsx('~/IMOA/LUAD/data/immueGenelist.xlsx')
gene_set <- gene_set[,c(1,2)]
colnames(gene_set) <- c('gene','type') 

# length(intersect(var_195_names,gene_set$gene))
exp_cluster1 <- exp[rownames(exp)%in%gene_set$gene,colnames(exp)%in%groupFile$sampleID[groupFile$group=='cluster1']]
exp_cluster2 <- exp[rownames(exp)%in%gene_set$gene,colnames(exp)%in%groupFile$sampleID[groupFile$group=='cluster2']]
immueExpr <- exp_cluster1
immueExpr <- exp_cluster2
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
ssGSEA_Score <- gsva(as.matrix(immueExpr),list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}

norm_ssGSEA_Score1 <- normalization(ssGSEA_Score)%>%t%>%as.data.frame()
norm_ssGSEA_Score1$group <- rep('cluster1',nrow(norm_ssGSEA_Score1))
norm_ssGSEA_Score2 <- normalization(ssGSEA_Score)%>%t%>%as.data.frame()
norm_ssGSEA_Score2$group <- rep('cluster2',nrow(norm_ssGSEA_Score2))
norm_ssGSEA_Score <- rbind(norm_ssGSEA_Score1,norm_ssGSEA_Score2)


library(tidyr)
fd <-norm_ssGSEA_Score %>% 
  rownames_to_column('sample')%>%
  pivot_longer(cols=c(colnames(norm_ssGSEA_Score)[-29]),names_to = 'Celltype',values_to = 'Composition')
fd
im <- ggboxplot(fd,x = 'Celltype', y='Composition',fill='group')+theme_classic()+theme(axis.text.x = element_text(angle = 90,hjust = 0.8,vjust = 0.9))+ylab('')+xlab('')
+geom_signif(comparisons = list(c("cluster1", "cluster2")),test = t.test,map_signif_level = T)
ggsave('CRG/result/figs/16_immue.tiff',im,width = 8,height = 6,dpi=1000)

p <- c()
for(i in 1:28){
  a <- wilcox.test(norm_ssGSEA_Score[,i]~norm_ssGSEA_Score$group)
  p <- c(p,a$p.value)
}
pval <- as.data.frame(cbind(colnames(norm_ssGSEA_Score)[1:28],p))
pval$p <- as.numeric(pval$p)
pval$group <- ifelse(pval$p>=0.05,"ns",ifelse(pval$p<0.05&&pval$p>=0.01,'*',ifelse(pval$p<0.01&&pval$p>=0.001,"**","***")))
library(writexl)
write_xlsx(pval,path = 'CRG/result/files/07_immue_pval.xlsx')

write.table(exp_cluster1 ,file="CRG/extdata/exp_cluster1.txt",sep = '\t',quote = F)
filterCommonGenes(input.f="CRG/extdata/exp_cluster1.txt",
                  output.f="CRG/extdata/exp_cluster1.gct",
                  id="GeneSymbol")

estimateScore(input.ds = "CRG/extdata/exp_cluster1.gct", 
              output.ds="CRG/extdata/estimateScore_cluster1.gct",
              platform="illumina")
scores1 <- read.table("CRG/extdata/estimateScore_cluster1.gct",skip = 2,header = T)
scores1 <- scores1 %>% t %>% as.data.frame()
colnames(scores1) <- scores1[1,]
scores1 <- scores1[-c(1,2),]
scores1$sample <- gsub('\\.','-',rownames(scores1))
scores1$group <- rep("cluster1",nrow(scores1))
write.table(exp_cluster2 ,file="CRG/extdata/exp_cluster2.txt",sep = '\t',quote = F)
filterCommonGenes(input.f="CRG/extdata/exp_cluster2.txt",
                  output.f="CRG/extdata/exp_cluster2.gct",
                  id="GeneSymbol")

estimateScore(input.ds = "CRG/extdata/exp_cluster2.gct", 
              output.ds="CRG/extdata/estimateScore_cluster2.gct",
              platform="illumina")
scores2 <- read.table("CRG/extdata/estimateScore_cluster2.gct",skip = 2,header = T)
scores2 <- scores2 %>% t %>% as.data.frame()
colnames(scores2) <- scores2[1,]
scores2 <- scores2[-c(1,2),]
scores2$sample <- gsub('\\.','-',rownames(scores2))
scores2$group <- rep("cluster2",nrow(scores2))

scores <- rbind(scores1,scores2)
scores$StromalScore <- as.numeric(scores$StromalScore)
scores$ImmuneScore <- as.numeric(scores$ImmuneScore)
scores$ESTIMATEScore <- as.numeric(scores$ESTIMATEScore)
ss <- ggboxplot(scores,x="group",y="StromalScore",fill = "group")+theme_classic()+ylab('StromalScore')+xlab('')+
  geom_signif(textsize=5,comparisons = list(c("cluster1", "cluster2")),test = t.test,map_signif_level = F,)+ 
  theme(axis.text.x = element_text(size=16,face='bold',color='black'),
        axis.text.y = element_text(size=16,face='bold',color='black'),
        text=element_text(size=16,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1))
is <- ggboxplot(scores,x='group',y='ImmuneScore',fill = "group")+theme_classic()+ylab('ImmuneScore')+xlab('')+
  geom_signif(textsize=5,comparisons = list(c("cluster1", "cluster2")),test = t.test,map_signif_level = F)+ 
  theme(axis.text.x = element_text(size=16,face='bold',color='black'),
        axis.text.y = element_text(size=16,face='bold',color='black'),
        text=element_text(size=16,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1))
es <- ggboxplot(scores,x='group',y='ESTIMATEScore',fill = "group")+theme_classic()+ylab('ESTIMATEScore')+xlab('')+
  geom_signif(textsize=5,comparisons = list(c("cluster1", "cluster2")),test = t.test,map_signif_level = F)+ 
  theme(axis.text.x = element_text(size=16,face='bold',color='black'),
        axis.text.y = element_text(size=16,face='bold',color='black'),
        text=element_text(size=16,face='bold'),
        legend.title = element_blank(),legend.position = "none",
        axis.line = element_line(size=1))
ggsave('CRG/result/figs/16_stromalscore.tiff',ss,width = 5,height = 5,dpi=1000)
ggsave('CRG/result/figs/16_immuescore.tiff',is,width = 5,height = 5,dpi=1000)
ggsave('CRG/result/figs/16_estimate.tiff',es,width = 5,height = 5,dpi=1000)
##############################MHC
HLA <- c("HLA-G","HLA-F","HLA-E","HLA-DRB5","HLA-DRB1",
         "HLA-DRA","HLA-DQB2","HLA-DQB1","HLA-DQA2","HLA-DQA1",
         "HLA-DPB1","HLA-DPA1","HLA-DOB","HLA-DOA","HLA-DMB","HLA-DMA","HLA-C","HLA-B","HLA-A")
HLA_c1 <- exp[rownames(exp)%in%HLA,colnames(exp)%in%groupFile$sampleID[groupFile$group=='cluster1']]
HLA_c2 <- exp[rownames(exp)%in%HLA,colnames(exp)%in%groupFile$sampleID[groupFile$group=='cluster2']]
HLA_c1 <- HLA_c1 %>% t %>% as.data.frame()
HLA_c1$sample <- rownames(HLA_c1)
HLA_c1$group <- rep("cluster1",nrow(HLA_c1))
HLA_c2 <- HLA_c2 %>% t %>% as.data.frame()
HLA_c2$sample <- rownames(HLA_c2)
HLA_c2$group <- rep("cluster2",nrow(HLA_c2))
HLA_all <- rbind(HLA_c1,HLA_c2)

library(tidyr)
pd <-HLA_all %>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(HLA_all)[-c(20,21)]),names_to = 'Celltype',values_to = 'Composition')

pp <- ggboxplot(pd,x = 'Celltype', y='Composition',fill='group')+theme_bw()+coord_flip()+ylab('')+xlab('')+
  theme(axis.text.x = element_text(size=10,face='bold',color='black'),
        axis.text.y = element_text(size=13,face='bold',color='black'),
        text=element_text(size=14,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/16_HLA.tiff',pp,width = 6,height = 8,dpi=1000)
p <- c()
for(i in 1:19){
  a <- wilcox.test(HLA_all[,i]~HLA_all$group)
  p <- c(p,a$p.value)
}
pval <- as.data.frame(cbind(colnames(HLA_all)[1:19],p))
pval$p <- as.numeric(pval$p)
pval$group <- ifelse(pval$p>=0.05,"ns",ifelse(pval$p<0.05&&pval$p>=0.01,'*',ifelse(pval$p<0.01&&pval$p>=0.001,"**","***")))
library(writexl)
write_xlsx(pval,path = 'CRG/result/files/08_HLA_pval.xlsx')
##################################TCAF
TCAF <- c('TNFRSF9','TNFRSF8','TNFRSF4','TNFRSF29','TNFRSF18','TNFRSF14',
          'ICOS','CD40LG','CD28','CD27','CD226','CD2')
TCAF_c1 <- exp[rownames(exp)%in%TCAF,colnames(exp)%in%groupFile$sampleID[groupFile$group=='cluster1']]
TCAF_c2 <- exp[rownames(exp)%in%TCAF,colnames(exp)%in%groupFile$sampleID[groupFile$group=='cluster2']]
TCAF_c1 <- TCAF_c1 %>% t %>% as.data.frame()
TCAF_c1$sample <- rownames(TCAF_c1)
TCAF_c1$group <- rep("cluster1",nrow(TCAF_c1))
TCAF_c2 <- TCAF_c2 %>% t %>% as.data.frame()
TCAF_c2$sample <- rownames(TCAF_c2)
TCAF_c2$group <- rep("cluster2",nrow(TCAF_c2))
TCAF_all <- rbind(TCAF_c1,TCAF_c2)
library(tidyr)
library(ggplot2)
td <-TCAF_all %>% 
  rownames_to_column('Sample')%>%
  pivot_longer(cols=c(colnames(TCAF_all)[-c(12,13)]),names_to = 'Celltype',values_to = 'Composition')

tp <- ggboxplot(td,x = 'Celltype', y='Composition',fill='group')+theme_bw()+coord_flip()+ylab('')+xlab('')+
  theme(axis.text.x = element_text(size=10,face='bold',color='black'),
        axis.text.y = element_text(size=13,face='bold',color='black'),
        text=element_text(size=14,face='bold'),
        legend.title = element_blank(),legend.position = "top",
        panel.grid = element_blank(),panel.border = element_rect(size=2))
ggsave('CRG/result/figs/16_TCAF.tiff',tp,width = 6,height = 8,dpi=1000)
p <- c()
for(i in 1:11){
  a <- wilcox.test(TCAF_all[,i]~TCAF_all$group)
  p <- c(p,a$p.value)
}
pval <- as.data.frame(cbind(colnames(TCAF_all)[1:11],p))
pval$p <- as.numeric(pval$p)
pval$group <- ifelse(pval$p>=0.05,"ns",ifelse(pval$p<0.05&&pval$p>=0.01,'*',ifelse(pval$p<0.01&&pval$p>=0.001,"**","***")))
library(writexl)
write_xlsx(pval,path = 'CRG/result/files/09_TCAF_pval.xlsx')

save(TCAF_all,HLA_all,scores,norm_ssGSEA_Score,groupFile,clin_cluster,result_hc,tmb,file='CRG/extdata/Subtype.Rdata')
load('CRG/extdata/Subtype.Rdata')

##############表达热图和临床特征
library(pheatmap)
load('LUADdata/TCGA/clinic.Rdata')
clin_tcga <- CLINIC[,c("sampleID","gender","age_at_initial_pathologic_diagnosis","pathologic_stage","pathologic_T","pathologic_N","pathologic_M","vital_status","tobacco_smoking_history_indicator")]
colnames(clin_tcga) <- c("sampleID",'Gender','Age','Stage',"Tumor","Lymph Node","Metastasis",'Status','Smoking')
clin_tcga <- as.data.frame(clin_tcga)
rownames(clin_tcga) <- clin_tcga$sampleID
clin_pheat<- clin_tcga[clin_tcga$sampleID%in%clin_cluster$sampleID,]
group <- merge(clin_cluster,clin_pheat,by.x='sampleID')
group <-group[,-c(2:6,8)]
group <- as.data.frame(group)
rownames(group) <- group$sampleID
colnames(group)[2] <- 'Cluster'
group <- group[order(group$Cluster),]
group$Gender <- ifelse(group$Gender=='MALE','Male','Female')
group$Stage <- ifelse(group$Stage%in%c('Stage I','Stage IA','Stage IB'),'I',
                      ifelse(group$Stage%in%c('Stage II','Stage IIA','Stage IIB'),"II",
                             ifelse(group$Stage%in%c('Stage IIIA','Stage IIIB'),"III",
                                    ifelse(group$Stage%in%c('Stage IV'),'IV','UNKOWN'))))
group$Tumor <- ifelse(group$Tumor%in%c('T1', 'T1a', 'T1b'),'T1',ifelse(group$Tumor%in%c('T2', 'T2a', 'T2b'),'T2',"T3-4"))
group$`Lymph Node` <- ifelse(group$`Lymph Node`%in%c('N0'),'N0',ifelse(group$`Lymph Node`%in%c('N1'),'N1',"N2-3"))
group$Metastasis <- ifelse(group$Metastasis%in%c('M0'),'M0','M1')
group$Smoking <- ifelse(group$Smoking%in%c(''),'UNKOWN',ifelse(group$Smoking%in%c('Lifelong Non-smoker'),'Never smoker','Current smoker'))
group$Status <- ifelse(group$Status%in%c('LIVING'),'Alive','Dead')
group$Age <- ifelse(group$Age>70,'old','young')
group <- na.omit(group )
Group <- group[,-1]
exp_pheat <- exp_CRGs_tumor[exp_CRGs_tumor$sampleID%in%group$sampleID,]
exp_pheat$sampleID <- rownames(exp_pheat)
data <- merge(group,exp_pheat,by.x = 'sampleID')
p <- c()
for(i in 11:ncol(data)){
  a <- wilcox.test(data[,i]~data$Cluster)
  p <- c(p,a$p.value)
}
pval <- as.data.frame(cbind(colnames(data)[11:ncol(data)],p))
pval$p <- as.numeric(pval$p)
pval$group <- ifelse(pval$p>=0.05,"ns",ifelse(pval$p<0.05&&pval$p>=0.01,'*',ifelse(pval$p<0.01&&pval$p>=0.001,"**","***")))
table(pval$group=="***") #178
top30_gene <- pval$V1[order(pval$p)][1:30]
library(writexl)
write_xlsx(pval,path = 'CRG/result/files/02_expCRGcluster_pval.xlsx')
rownames(data) <- data$sampleID
dat <- data[order(data$Cluster),colnames(data)%in%top30_gene]
dat1 <- as.data.frame(t(dat))

p <- pheatmap(dat1,show_colnames = F,cluster_cols =F,annotation_col = Group,cellheight = 10,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
save_pheatmap_pdf <- function(x, filename, width=7, height=12) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p, "CRG/result/figs/17_pheatcluster.pdf")

dat2 <- as.data.frame(t(data[,11:ncol(data)]))
colnames(dat2) <- data$sampleID
Group2 <- as.data.frame(Group[,1])
rownames(Group2) <- rownames(Group)
colnames(Group2) <- 'Group'
p2 <- pheatmap(dat2,show_colnames = F,show_rownames = F,cluster_cols =T,annotation_col = Group2,cellheight = 1.5,
              color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p2, "CRG/result/figs/17_pheatclusterall.pdf")

#sankey
sdf <- Group[,c(1,4,8)]
df1 <- as.data.frame(table(sdf[c(1,2)]))
colnames(df1)[1:2] <- c('Var1','Var2')
df2 <- as.data.frame(table(sdf[c(2,3)]))
colnames(df2)[1:2] <- c('Var1','Var2')
df1 <- rbind(df1,df2)

df1 <- df1[df1$Freq!="0",]
colnames(df1) <- c("source","target","value")
nodes <- data.frame(name=unique(c(df1$source,df1$target)))
df1$IDsource <- match(df1$source, nodes$name)-1
df1$IDtarget <- match(df1$target, nodes$name)-1
#用网站画比较好看
# library(networkD3)
# Sankey.p <- sankeyNetwork(Links = df1, Nodes = nodes,
#                           Source = "IDsource", Target = "IDtarget",
#                           Value = "value", NodeID = "name", 
#                           fontSize = 14, # 字体大小
#                           nodeWidth = 30, # node节点的宽度
#                           sinksRight=T)# 最后子节点文字标注方向

# htmlwidgets::saveWidget(Sankey.p, file='CRG/result/figs/17_sankey.html')

library(writexl)
write_xlsx(df1,'CRG/extdata/sankey.xlsx')
write_xlsx(sdf,'CRG/extdata/sankey_df.xlsx')

load('CRG/extdata/tcga_plot.Rdata')
sdf2 <- sdf 
g <- c()
for(i in rownames(sdf2)){
  a <- pd$Group[rownames(pd)==i]
  g <- c(g,a)
}
sdf2$Stage <- g
colnames(sdf2)[2] <- 'Risk'
df3<- as.data.frame(table(sdf2[c(1,2)]))
colnames(df3)[1:2] <- c('Var1','Var2')
df4 <- as.data.frame(table(sdf2[c(2,3)]))
colnames(df4)[1:2] <- c('Var1','Var2')
df3 <- rbind(df3,df4)

df3 <- df3[df3$Freq!="0",]
colnames(df3) <- c("source","target","value")
nodes <- data.frame(name=unique(c(df3$source,df3$target)))
df3$IDsource <- match(df3$source, nodes$name)-1
df3$IDtarget <- match(df3$target, nodes$name)-1

library(writexl)
write_xlsx(df3,'CRG/extdata/sankey2.xlsx')
write_xlsx(sdf2,'CRG/extdata/sankey_df2.xlsx')
########################F2 multicox
library(survival)
library(survminer)
pd$sampleID <- rownames(pd)
pd <- pd[pd$sampleID%in%data$sampleID,]
pd2 <- pd[,c(19,1,2)]
data2 <- data[,c(1,3,4,5,6,7,8,10,43)]
F2multidata <- merge(data2,pd2,by.x = 'sampleID')
colnames(F2multidata)[6] <- 'Lymph_Node'
mulcox<- coxph(Surv(time, status) ~ F2+Age+Stage+Tumor+Lymph_Node+Metastasis, data =F2multidata );mulcox
rcorrcens(Surv(time, status) ~ predict(mulcox), data = F2multidata)
# c-index 0.317

forest <- ggforest(mulcox,  #coxph得到的Cox回归结果
                   data = F2multidata,  #数据集
                   main = 'Hazard ratio of all data',  #标题
                   cpositions = c(0.05, 0.15, 0.35),  #前三列数值的距离
                   fontsize = 1.4, #字体大小
                   refLabel = 'Reference', #相对变量的标签
                   noDigits = 3 #HR、95%CI的小数位
)
ggsave(filename="CRG/result/figs/08_F2muti_forest.png",forest, width=17, height=10)
#####################################F2 KM TCGA
F2_surv_tcgadata <- F2multidata[,9:11]
F2_surv_tcgadata$Group <- ifelse(F2_surv_tcgadata$F2>median(F2_surv_tcgadata$F2),'High F2','Low F2')
library(writexl)
write_xlsx(F2_surv_tcgadata,path = 'CRG/extdata/surv_F2_TCGA.xlsx')
library(readxl)
F2_surv_tcgadata <- read_xlsx('CRG/extdata/surv_F2_TCGA.xlsx')
p <- ggsurvplot(survfit(Surv(time, status) ~ Group, data=F2_surv_tcgadata),surv.median.line = "hv",palette = c('red','blue','grey'),pval = T,pval.size=10)
######################riskscore unicox
library(survival)
library(survminer)
pd$sampleID <- rownames(pd)
pd <- pd[pd$sampleID%in%data$sampleID,]
pd1 <- pd[,c(19,1,2,15)]
data1 <- data[,c(1,3,4,5,6,7,8,10)]
unidata <- merge(data1,pd1,by.x = 'sampleID')
colnames(unidata)[6] <- 'Lymph_Node'
covariates <- c('Gender','Age','Stage','Smoking',"Tumor" ,  "Lymph_Node", "Metastasis", 'riskscore')
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data =  unidata)})
#提取HR，95%置信区间和p值
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         #获取HR
                         HR <-signif(x$coef[2], digits=2);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(x)
                       })
###########riskscore  multicox
mulcox<- coxph(Surv(time, status) ~ Age+Stage+Tumor+Lymph_Node+riskscore, data = unidata );mulcox
rcorrcens(Surv(time, status) ~ predict(mulcox), data = unidata)
# c-index 0.279

forest <- ggforest(mulcox,  #coxph得到的Cox回归结果
                   data = unidata,  #数据集
                   main = 'Hazard ratio of all data',  #标题
                   cpositions = c(0.05, 0.15, 0.35),  #前三列数值的距离
                   fontsize = 1.4, #字体大小
                   refLabel = 'Reference', #相对变量的标签
                   noDigits = 3 #HR、95%CI的小数位
)
ggsave(filename="CRG/result/figs/18_muti_forest.png",forest, width=17, height=10)

##列线数据
lxdata <- unidata
# lxdata$Age <- group$Age
# lxdata$Age <- as.numeric(lxdata$Age)
lxdata$riskscore <- as.numeric(lxdata$riskscore)
lxdata$Stage <- as.factor(lxdata$Stage)
lxdata$time <- as.numeric(lxdata$riskscore)*30
library(writexl)
write_xlsx(lxdata,'CRG/result/files/18_nomogramdata.xlsx')
library(readxl)
lxdata <- read_xlsx('CRG/result/files/18_nomogramdata.xlsx')
lxdata1 <- read_xlsx('CRG/result/files/18_nomogramdata.xlsx')
lxdata1$Age <- ifelse(lxdata$Age=='young',1,2 )
lxdata1$Age <- factor(lxdata1$Age,levels=c(1,2),labels=c('young','old'))
lxdata1$Stage <- ifelse(lxdata$Stage=='I',1,ifelse(lxdata$Stage=='II',2 ,ifelse(lxdata$Stage=='III',3,ifelse(lxdata$Stage=='IV',4,0))))
lxdata1$Stage <- factor(lxdata1$Stage,levels=c(0,1,2,3,4),labels=c('UNKOWN','I','II','III','IV'))
lxdata1$time <- lxdata1$time*30
dat <- lxdata1[,c(3,4,9,10,11)]
library(rms)
dd <- datadist(dat)
options(datadist='dd')

m2<- psm(Surv(time, status) ~ Age+Stage+riskscore, data = dat,dist='lognormal' )
med <- Quantile(m2)
surv <- Survival(m2)
nomo <- nomogram(m2,fun=list(function(x)surv(365,x),function(x)surv(1095,x),function(x)surv(1825,x)), funlabel=c("1-year survival probability","3-year survival probability","5-year survival probability"))
png("CRG/result/figs/18_nomo.tiff",width = 2000,height = 2000,res=200)
plot(nomo,xfrac=0.6, lwd=3,cex.axis=1.3,cex.lab=1.5,font.lab=4)

cindex <- rcorrcens(Surv(time, status) ~ predict(m2), data = dat)
plot(cindex )
#0.722
f2 <- psm(Surv(time, status) ~ Age+Stage+riskscore, data = dat,x=T,y=T,dist='lognormal' )
call <- calibrate(f2,cmethod = 'KM',method='boot',u=365,m=164,B=1000)
png("CRG/result/figs/18_cal_1y.tiff",width = 1500,height = 1500,res=200)
plot(call,lty=1,conf.int=T,errbar.col='blue',col='red',xlim=c(0,1),ylim=c(0,1), lwd=3,cex.axis=1.5,ann = F,subtitles=F)
title(xlab='Nomogram-predicted OS of 1-year',ylab = 'Actual 1-year OS',cex.lab=1.5,font.lab=4)

call2 <- calibrate(f2,cmethod = 'KM',method='boot',u=1095,m=164,B=1000)
png("CRG/result/figs/18_cal_3y.tiff",width = 1500,height = 1500,res=200)
plot(call2,lty=1,conf.int=T,errbar.col='blue',col='red',xlim=c(0,1),ylim=c(0,1), lwd=3,cex.axis=1.5,ann = F,subtitles=F)
title(xlab='Nomogram-predicted OS of 3-year',ylab = 'Actual 3-year OS',cex.lab=1.5,font.lab=4)

call3 <- calibrate(f2,cmethod = 'KM',method='boot',u=1825,m=164,B=1000)
png("CRG/result/figs/18_cal_5y.tiff",width = 1500,height = 1500,res=200)
plot(call3,lty=1,conf.int=T,errbar.col='blue',col='red',xlim=c(0,1),ylim=c(0,1), lwd=3,cex.axis=1.5,ann = F,subtitles=F)
title(xlab='Nomogram-predicted OS of 5-year',ylab = 'Actual 5-year OS',cex.lab=1.5,font.lab=4)

library(pec)
library(survival)
dat1 <- dat
dat1$time <- dat1$time/30
cox1 <- coxph(Surv(time, status) ~ Age+Stage+riskscore, data = dat1,x=T,y=T)
cox2 <- coxph(Surv(time, status) ~ Age, data = dat1,x=T,y=T)
cox3 <- coxph(Surv(time, status) ~ Stage, data = dat1,x=T,y=T)
cox4 <- coxph(Surv(time, status) ~ riskscore, data = dat1,x=T,y=T)
ApparrentCindex  <- pec::cindex(list("Age"=cox2,
                                     "Stage"=cox3,"risk score"=cox4,"Combined model"=cox1),
                                formula=Surv(time,status)~Age+Stage+riskscore,
                                data=dat1,
                                eval.times=seq(1,242,20))
print(ApparrentCindex)
png("CRG/result/figs/18_cindex.tiff",width = 1500,height = 1500,res=200)
plot(ApparrentCindex,ann = F,lwd=3,cex.axis=1.5,xlab='Time(month)',ylab = 'Concordance index')


##########AUC combinied riskscore age stage
####1-year
coxph(Surv(time, status) ~ age, data = lxdata )
coxph(Surv(time, status) ~ stage, data = lxdata )
coxph(Surv(time, status) ~ age+stage+riskscore, data = lxdata )
lxdata$age <- ifelse(lxdata$Age=='young',1,2)
lxdata$stage<- ifelse(lxdata$Stage=='I',1, ifelse(lxdata$Stage=='II',2, ifelse(lxdata$Stage=='III',3, ifelse(lxdata$Stage=='IV',4,0))))
lxdata$riskscore1 <- -0.3797 *as.numeric(lxdata$age)
lxdata$riskscore2 <- 0.48652*as.numeric(lxdata$stage)
lxdata$riskscore3 <-  0.43012 *lxdata$age+0.36518*lxdata$stage+0.92599*lxdata$riskscore

library(survivalROC)
troc= survivalROC(Stime=lxdata$time,  
                  status=lxdata$status,   
                  marker = lxdata$riskscore3,     
                  predict.time = 12, method="KM")
troc.1= survivalROC(Stime=lxdata$time,  
                    status=lxdata$status,   
                    marker = lxdata$riskscore,     
                    predict.time = 12, method="KM")
troc.2= survivalROC(Stime=lxdata$time,  
                    status=lxdata$status,   
                    marker = lxdata$riskscore1,     
                    predict.time = 12, method="KM")
troc.3= survivalROC(Stime=lxdata$time,  
                    status=lxdata$status,   
                    marker = lxdata$riskscore2,     
                    predict.time = 12, method="KM")

png("CRG/result/figs/18_1-yearroc_combined.tiff",width = 1500,height = 1500,res=200)
plot(troc$FP, troc$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=3,cex.axis=1.5,ann = F)
title(xlab="1-Specificity", ylab="Sensitivity",cex.lab=1.5,font.lab=4)
abline(0,1,col='gray',lty=2,lwd=2)
lines(troc.1$FP, troc.1$TP, type="l",lwd=2,col="blue", xlim=c(0,1), ylim=c(0,1))
lines(troc.2$FP, troc.2$TP, type="l",lwd=2,col="green", xlim=c(0,1), ylim=c(0,1))
lines(troc.3$FP, troc.3$TP, type="l",lwd=2,col="purple", xlim=c(0,1), ylim=c(0,1))
legend(0.45,0.3,c(paste('Combined model=',round(troc$AUC,3)),
                  paste('Riskscore=',round(troc.1$AUC,3)),
                  paste('Age=',round(troc.2$AUC,3)),
                  paste('Stage=',round(troc.3$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=2,col = c('red','blue',"green","purple"),
       bty='n',seg.len = 2,cex = 1.2)
dev.off()
####3-year
troc= survivalROC(Stime=lxdata$time,  
                  status=lxdata$status,   
                  marker = lxdata$riskscore3,     
                  predict.time = 36, method="KM")
troc.1= survivalROC(Stime=lxdata$time,  
                    status=lxdata$status,   
                    marker = lxdata$riskscore,     
                    predict.time = 36, method="KM")
troc.2= survivalROC(Stime=lxdata$time,  
                    status=lxdata$status,   
                    marker = lxdata$riskscore1,     
                    predict.time = 36, method="KM")
troc.3= survivalROC(Stime=lxdata$time,  
                    status=lxdata$status,   
                    marker = lxdata$riskscore2,     
                    predict.time = 36, method="KM")

png("CRG/result/figs/18_3-yearroc_combined.tiff",width = 1500,height = 1500,res=200)
plot(troc$FP, troc$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=3,cex.axis=1.5,ann = F)
title(xlab="1-Specificity", ylab="Sensitivity",cex.lab=1.5,font.lab=4)
abline(0,1,col='gray',lty=2,lwd=2)
lines(troc.1$FP, troc.1$TP, type="l",lwd=2,col="blue", xlim=c(0,1), ylim=c(0,1))
lines(troc.2$FP, troc.2$TP, type="l",lwd=2,col="green", xlim=c(0,1), ylim=c(0,1))
lines(troc.3$FP, troc.3$TP, type="l",lwd=2,col="purple", xlim=c(0,1), ylim=c(0,1))
legend(0.45,0.3,c(paste('Combined model=',round(troc$AUC,3)),
                  paste('Riskscore=',round(troc.1$AUC,3)),
                  paste('Age=',round(troc.2$AUC,3)),
                  paste('Stage=',round(troc.3$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=2,col = c('red','blue',"green","purple"),
       bty='n',seg.len = 2,cex = 1.2)
dev.off()
####5-year
troc= survivalROC(Stime=lxdata$time,  
                  status=lxdata$status,   
                  marker = lxdata$riskscore3,     
                  predict.time = 60, method="KM")
troc.1= survivalROC(Stime=lxdata$time,  
                    status=lxdata$status,   
                    marker = lxdata$riskscore,     
                    predict.time = 60, method="KM")
troc.2= survivalROC(Stime=lxdata$time,  
                    status=lxdata$status,   
                    marker = lxdata$riskscore1,     
                    predict.time = 60, method="KM")
troc.3= survivalROC(Stime=lxdata$time,  
                    status=lxdata$status,   
                    marker = lxdata$riskscore2,     
                    predict.time = 60, method="KM")

png("CRG/result/figs/18_5-yearroc_combined.tiff",width = 1500,height = 1500,res=200)
plot(troc$FP, troc$TP, type="l",col="red", xlim=c(0,1), ylim=c(0,1), lwd=3,cex.axis=1.5,ann = F)
title(xlab="1-Specificity", ylab="Sensitivity",cex.lab=1.5,font.lab=4)
abline(0,1,col='gray',lty=2,lwd=2)
lines(troc.1$FP, troc.1$TP, type="l",lwd=2,col="blue", xlim=c(0,1), ylim=c(0,1))
lines(troc.2$FP, troc.2$TP, type="l",lwd=2,col="green", xlim=c(0,1), ylim=c(0,1))
lines(troc.3$FP, troc.3$TP, type="l",lwd=2,col="purple", xlim=c(0,1), ylim=c(0,1))
legend(0.45,0.3,c(paste('Combined model=',round(troc$AUC,3)),
                  paste('Riskscore=',round(troc.1$AUC,3)),
                  paste('Age=',round(troc.2$AUC,3)),
                  paste('Stage=',round(troc.3$AUC,3))),
       x.intersp =1,y.intersp =1,text.font = 4,lty = 1,lwd=2,col = c('red','blue',"green","purple"),
       bty='n',seg.len = 2,cex = 1.2)
dev.off()

###########GSVA
library(tidyverse)
library(clusterProfiler)
library(msigdbr)  #install.packages("msigdbr")
library(GSVA) 
library(GSEABase)
library(pheatmap)
library(limma)
library(BiocParallel)
CRGs <- read_xlsx('CRG/result/files/CRGs_209.xlsx')
KEGG_df_all <-  msigdbr(species = "Homo sapiens",category = 'C2',subcategory = "CP:KEGG")
KEGG_df <- dplyr::select(KEGG_df_all,gs_name,gs_exact_source,gene_symbol)
kegg_list <- split(KEGG_df$gene_symbol, KEGG_df$gs_name) ##按照gs_name给gene_symbol分组

gsva_mat <- gsva(expr=as.matrix(exp_tumor), 
                 gset.idx.list=kegg_list , 
                 kcdf="Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                 verbose=T, 
                 parallel.sz = parallel::detectCores())#调用所有核
gsva_mat <- gsva_mat[,colnames(gsva_mat)%in%groupFile$sampleID]
write_xlsx(as.data.frame(gsva_mat),"CRG/extdata/gsva_go_matrix.xlsx")



exp="cluster1"
ctr="cluster2"

design <- model.matrix(~0+factor(gl$group))
colnames(design) <- levels(factor(gl$group))
rownames(design) <- colnames(gsva_mat)
contrast.matrix <- makeContrasts(contrasts=paste0(exp,'-',ctr),  #"exp/ctrl"
                                 levels = design)

fit1 <- lmFit(gsva_mat,design)                 #拟合模型
fit2 <- contrasts.fit(fit1, contrast.matrix) #统计检验
efit <- eBayes(fit2) 
summary(decideTests(efit,lfc=1, p.value=0.05)) #统计查看差异结果
tempOutput <- topTable(efit, coef=paste0(exp,'-',ctr), n=Inf)
degs <- na.omit(tempOutput) 
write.csv(degs,"CRG/extdata/gsva_degs.results.csv")
padj_cutoff=0.05
log2FC_cutoff=log2(2)

keep <- rownames(degs[degs$logFC>0.03, ])
length(keep)
str_to_title(gsub('_',' ',keep))
dat <- gsva_mat[keep,] 
dat <- as.data.frame(dat)

gl=groupFile
gl <- gl[,-1]
gl <- as.data.frame(gl)
rownames(gl) <- groupFile$sampleID
gl <- gl[order(gl$group),]
gl<- as.data.frame(gl)
rownames(gl) <- groupFile$sampleID[order(groupFile$group)]
colnames(gl) <- 'group'

dat1 <- dat[,rownames(gl)]
rownames(dat1) <- gsub('KEGG_','',rownames(dat1))

p <- pheatmap(dat1,show_colnames = F,cluster_cols =F,annotation_col = gl,cellheight = 12,
              color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
save_pheatmap_pdf <- function(x, filename, width=12, height=6) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p, "CRG/result/figs/18_GSVA_pheatcluster.pdf")







