###################################################LASSO COX MODEL###################################################
LASSO_COX_MODEL <- function(CRGs_surv_dat,plot){
  library(glmnet)
  library(survival)
  library(survminer)
  set.seed(10)
  X.train <- as.matrix(CRGs_surv_dat[,4:ncol(CRGs_surv_dat)])
  Y.train <- as.matrix(data.frame(time=CRGs_surv_dat$OS,status=CRGs_surv_dat$OS.E))
  fit.cv <- cv.glmnet(X.train,Y.train, family = "cox", nfolds = 10,alpha=1,keep=TRUE)
  print(paste0('lambda.min: ',fit.cv$lambda.min))
 
  idmin = match(fit.cv$lambda.min, fit.cv$lambda)
  lsocoef<-as.numeric(coef(fit.cv ,s="lambda.min"))
  lasso_genes <- colnames(X.train)[lsocoef!=0]
  
  if(plot=='T'){
    png("figures/02B_lassocoef.png",width = 1500,height = 1500,res=200)
    plot(fit.cv, lwd=3,cex.axis=1.5,ann = F)
    title(xlab="log(λ)", ylab="Partial likelihood Deviance",cex.lab=1.5,font.lab=2)
    dev.off()
    
    model <- glmnet(X.train,Y.train, family = "cox", alpha=1)
    

    png("figures/02A_coefficient.png",width = 1500,height = 1500,res=200)
    plot(model, lwd=3,cex.axis=1.5,ann = F)
    title(xlab="log(λ)", ylab="Coefficients",cex.lab=1.5,font.lab=2)
    dev.off()
    
    
    dmf <- read.xlsx(file = 'files/Fig2.xlsx',sheetName = 'Fig2D')
    dmf <- dmf[1:11,1:6]
    dmf <- dmf[order(dmf$coef,decreasing = T),]
    dmf$gene <- factor(dmf$gene,levels=dmf$gene)
    # Fig2C
    library(ggalt)
    dp <- ggplot(dmf,aes(y = gene, x=coef,xend=0))+
      geom_segment(aes(x=coef,xend=0,y=gene,yend=gene),color='grey',size=2)+
      geom_dumbbell(size_x=6, size_xend = 0,color="grey",colour_x = "black")+
      labs(x='Coefficient', y=NULL) +geom_vline(aes(xintercept = 0),lty=2,lwd=1)+
      theme_bw() + theme(axis.text.y = element_text(size=16,face='bold',color='black'),
                         axis.text.x = element_text(size=15,face='bold',color='black'),
                         text=element_text(size=16,face='bold'),axis.line = element_line(size=1),panel.grid = element_blank(),panel.border = element_rect(size=1.5))
    
    ggsave('figures/02C_dot_coef_bsr.png',dp,width = 6,height = 6,dpi=1000)
  }   
  return(lasso_genes)
}



lasso_genes <- LASSO_COX_MODEL(CRGs_surv_dat,plot='F')






