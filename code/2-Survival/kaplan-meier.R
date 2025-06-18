library(WGCNA)
library(survival)
library(survminer)

load(file="data/RData/metabric.RData")
load(file="data/RData/meta.RData")
load(file="data/Complete_METABRIC_Clinical_Survival_Data__DSS.rbin")
load(file="results/2025/centrality_global.RData")
load(file="results/2025/centrality_basal.RData")

#############compute module eigengenes

MEs= moduleEigengenes(t(metabric), centrality_global$module)$eigengenes
colnames(MEs)<-gsub("^ME","", colnames(MEs))
MEs_basal= moduleEigengenes(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]), centrality_basal$module)$eigengenes
colnames(MEs_basal)<-gsub("^ME","", colnames(MEs_basal))

####DSS
for(i in 1:ncol(MEs)){
  a=data.frame(ID=colnames(metabric),OS=Complete_METABRIC_Clinical_Survival_Data__DSS[,1], IND=Complete_METABRIC_Clinical_Survival_Data__DSS[,2])
  b=data.frame(ID=rownames(MEs), eigen_gene=MEs[,i])
  module<-gsub("ME", "", colnames(MEs)[i])
  c=cbind(a,b)
  
  eigen_class<-cut(MEs[,i], breaks=quantile(MEs[,i],probs=c(0,1,0.5), na.rm=T))
  
  fit<-survfit(Surv(OS,IND)~eigen_class, data=c)
  ggsurv <- ggsurvplot(fit, data=c, risk.table=F, pval=T, palette=c("blue", "red"),legend.labs=c(paste("low", module, collapse=" "), paste("high",module, collapse=" ")), censor=F, tables.theme = clean_theme())
  png(paste("results/2025/survival_metabric_DSS_",module,".png", collapse=""), res=300, 1500, 1500)
  ggsurv$plot <- ggsurv$plot +
    theme(legend.text = element_text(size = 10))
  print(ggsurv)
  dev.off()
  
}



########5years surv
surv5years<-Complete_METABRIC_Clinical_Survival_Data__DSS
surv5years[Complete_METABRIC_Clinical_Survival_Data__DSS[,1]>1825,"status"]<-0
surv5years[Complete_METABRIC_Clinical_Survival_Data__DSS[,1]>1825,"time"]<-1825

for(i in 1:ncol(MEs)){
  a=data.frame(ID=colnames(metabric),OS=surv5years[,1], IND=surv5years[,2])
  b=data.frame(ID=rownames(MEs), eigen_gene=MEs[,i])
  module<-gsub("ME", "", colnames(MEs)[i])
  c=cbind(a,b)
  
  eigen_class<-cut(MEs[,i], breaks=quantile(MEs[,i],probs=c(0,1,0.5), na.rm=T))
  
  fit<-survfit(Surv(OS,IND)~eigen_class, data=c)
  ggsurv <- ggsurvplot(fit, data=c, risk.table=F, pval=T, palette=c("blue", "red"),legend.labs=c("low ME", "high ME"), censor=F, tables.theme = clean_theme(), ylim=c(0.6, 1))
  png(paste("results/2025/survival_metabric_DSS_",module,"_5years.png", collapse=""), res=300, 1000, 1000)
  ggsurv$plot <- ggsurv$plot +
    theme(legend.text = element_text(size = 10))
  print(ggsurv)
  dev.off()
  
}



TFblue<-read.csv("TF_blue_basal.csv", header=F, stringsAsFactors = F)
for(gene in TFblue[,1]){
  a=data.frame(ID=colnames(metabric),OS=surv5years[,1], IND=surv5years[,2])
  b=data.frame(eigen_gene=metabric[gene,])
  c=cbind(a,b)
  eigen_class<-cut(metabric[gene,], breaks=quantile(metabric[gene,],probs=c(0,1,0.5), na.rm=T))
  
  
  fit<-survfit(Surv(OS,IND)~eigen_class, data=c)
  ggsurv <- ggsurvplot(fit, data=c, risk.table=F, pval=T, legend.labs=c(paste("low", gene, collapse=" "), paste("high",gene, collapse=" ")), censor=F, tables.theme = clean_theme(), ylim=c(0.6,1), pval.coord=c(200,0.7), palette=c("blue", "red"))
  png(paste("survival_metabric_DSS_5years",gene,".png", collapse=""), 800, 800, res=200)
  ggsurv$plot <- ggsurv$plot + 
    theme(legend.text = element_text(size = 14))
  print(ggsurv)
  dev.off()
  
}


for(gene in TFblue[,1]){
  a=data.frame(ID=colnames(metabric),OS=surv5years[,1], IND=surv5years[,2])
  b=data.frame(eigen_gene=metabric[gene,])
  c=cbind(a,b)
  eigen_class<-cut(metabric[gene,], breaks=quantile(metabric[gene,],probs=c(0,1,0.5), na.rm=T))
  
  
  fit<-survfit(Surv(OS,IND)~eigen_class, data=c)
  pp<-survdiff(Surv(OS,IND)~eigen_class, data=c)
  pval<-1-pchisq(pp$chisq, length(pp$n)-1)
  print(paste(gene, signif(pval, digits = 2)))
  
}


for(i in 1:21){
  a=data.frame(ID=colnames(metabric),OS=surv5years[,1], IND=surv5years[,2])
  b=data.frame(ID=rownames(MEs_metabric), eigen_gene=MEs_metabric[,i])
  module<-gsub("ME", "", colnames(MEs_metabric)[i])
  c=cbind(a,b)
  fit=survfit(Surv(OS,IND)~1, data=c)
  plot(fit)
  
  #eigen_class<-cut(MEs_brain[type=="Basal",i],2)
  eigen_class<-cut(MEs_metabric[,i], breaks=quantile(MEs_metabric[,i],probs=c(0,1,0.5), na.rm=T))
  
  
  fit<-survfit(Surv(OS,IND)~eigen_class, data=c)
  pp<-survdiff(Surv(OS,IND)~eigen_class, data=c)
  pval<-1-pchisq(pp$chisq, length(pp$n)-1)
  print(paste(module, alt_names[match(module, labels2colors(1:20))], signif(pval, digits = 2)))
  
}
