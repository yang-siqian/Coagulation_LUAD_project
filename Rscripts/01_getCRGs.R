library(dplyr)
library(writexl)
library(KEGGREST)
ccc <- keggGet('hsa04610') #Complement and coagulation cascades
pa <- keggGet('hsa04611')  #Platelet activation
ccc_gs <- unlist(lapply(ccc[[1]]$GENE,function(x) strsplit(x,';')))
ccc_gl <- ccc_gs[1:length(ccc_gs)%%3==2]
pa_gs <- unlist(lapply(pa[[1]]$GENE,function(x) strsplit(x,';')))
pa_gl <- pa_gs[1:length(pa_gs)%%3==2]
genelist <- c(ccc_gl,pa_gl) %>% as.data.frame()
write_xlsx(genelist,path = 'CRG/result/files/CRGs_209.xlsx')