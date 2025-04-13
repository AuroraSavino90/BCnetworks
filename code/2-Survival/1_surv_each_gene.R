##########survival each gene
###############
load("data/Complete_METABRIC_Clinical_Survival_Data__DSS.rbin")
load("data/RData/metabric.RData")
load("data/RData/meta.RData")

#load survival data, cut at 5years
surv<-Complete_METABRIC_Clinical_Survival_Data__DSS
surv5years<-surv
surv5years[surv[,1]>1825,"status"]<-0
surv5years[surv[,1]>1825,"time"]<-1825


#####survival######
#for each gene test the relationship with survival
#pvalue saved in two different vectors based on whether it correlates with poor or good survival

library(survival)
a=data.frame(ID=colnames(metabric),OS=surv[,1], IND=surv[,2])
pvalue_good<-c()
pvalue_poor<-c()

for(i in 1:nrow(metabric)){
  gene<-rownames(metabric)[i]
  b=data.frame(eigen_gene=metabric[gene,])
  c=cbind(a,b)
  
  eigen_class<-cut(metabric[gene,], breaks=quantile(metabric[gene,],probs=c(0,1,0.5), na.rm=T))
  
  pp<-survdiff(Surv(OS,IND)~eigen_class, data=c)
  pval<-1-pchisq(pp$chisq, length(pp$n)-1)
  
  if(pp$obs[1]>pp$exp[1]){
    pvalue_good<-c(pvalue_good, pval)
    pvalue_poor<-c(pvalue_poor, NA)
  } else {
    pvalue_good<-c(pvalue_good, NA)
    pvalue_poor<-c(pvalue_poor, pval)
  }
}


direction<-rep("None", nrow(metabric))
direction[!is.na(pvalue_good)]<-"Good"
direction[!is.na(pvalue_poor)]<-"Poor"


pval<-rep(NA, nrow(metabric))
pval[is.na(pvalue_good)==F]<-pvalue_good[is.na(pvalue_good)==F]
pval[is.na(pvalue_poor)==F]<-pvalue_poor[is.na(pvalue_poor)==F]

surv_data_global_all<-data.frame(gene=rownames(metabric), pval=pval, direction=direction)
rownames(surv_data_global_all)<-rownames(metabric)
save(surv_data_global_all, file="results/2025/surv_data_global_all.RData")



#####5years
a=data.frame(ID=colnames(metabric),OS=surv5years[,1], IND=surv5years[,2])
pvalue_good<-c()
pvalue_poor<-c()

for(i in 1:nrow(metabric)){
  gene<-rownames(metabric)[i]
  b=data.frame(eigen_gene=metabric[gene,])
  c=cbind(a,b)
  
  eigen_class<-cut(metabric[gene,], breaks=quantile(metabric[gene,],probs=c(0,1,0.5), na.rm=T))
  
  pp<-survdiff(Surv(OS,IND)~eigen_class, data=c)
  pval<-1-pchisq(pp$chisq, length(pp$n)-1)
  
  if(pp$obs[1]>pp$exp[1]){
    pvalue_good<-c(pvalue_good, pval)
    pvalue_poor<-c(pvalue_poor, NA)
  } else {
    pvalue_good<-c(pvalue_good, NA)
    pvalue_poor<-c(pvalue_poor, pval)
  }
}

direction<-rep("None", nrow(metabric))
direction[!is.na(pvalue_good)]<-"Good"
direction[!is.na(pvalue_poor)]<-"Poor"

pval<-rep(NA, nrow(metabric))
pval[is.na(pvalue_good)==F]<-pvalue_good[is.na(pvalue_good)==F]
pval[is.na(pvalue_poor)==F]<-pvalue_poor[is.na(pvalue_poor)==F]

surv_data_global_5years<-data.frame(gene=rownames(metabric), pval=pval, direction=direction)
rownames(surv_data_global_5years)<-rownames(metabric)
save(surv_data_global_5years, file="results/2025/surv_data_global_5years.RData")



#######Subtypes
#####survival######
for(subtype in unique(meta$NOT_IN_OSLOVAL_Pam50Subtype)){
a=data.frame(ID=colnames(metabric)[meta$NOT_IN_OSLOVAL_Pam50Subtype==subtype],OS=surv[meta$NOT_IN_OSLOVAL_Pam50Subtype==subtype,1], IND=surv[meta$NOT_IN_OSLOVAL_Pam50Subtype==subtype,2])
pvalue_good<-c()
pvalue_poor<-c()

for(i in 1:nrow(metabric)){
  gene<-rownames(metabric)[i]
  b=data.frame(eigen_gene=metabric[gene,meta$NOT_IN_OSLOVAL_Pam50Subtype==subtype])
  c=cbind(a,b)
  
  eigen_class<-cut(metabric[gene,meta$NOT_IN_OSLOVAL_Pam50Subtype==subtype], breaks=quantile(metabric[gene,meta$NOT_IN_OSLOVAL_Pam50Subtype==subtype],probs=c(0,1,0.5), na.rm=T))
  
  pp<-survdiff(Surv(OS,IND)~eigen_class, data=c)
  pval<-1-pchisq(pp$chisq, length(pp$n)-1)
  
  if(pp$obs[1]>pp$exp[1]){
    pvalue_good<-c(pvalue_good, pval)
    pvalue_poor<-c(pvalue_poor, NA)
  } else {
    pvalue_good<-c(pvalue_good, NA)
    pvalue_poor<-c(pvalue_poor, pval)
  }
}

direction<-rep("None", nrow(metabric))
direction[!is.na(pvalue_good)]<-"Good"
direction[!is.na(pvalue_poor)]<-"Poor"

pval<-rep(NA, nrow(metabric))
pval[is.na(pvalue_good)==F]<-pvalue_good[is.na(pvalue_good)==F]
pval[is.na(pvalue_poor)==F]<-pvalue_poor[is.na(pvalue_poor)==F]

surv_data_global_subtype<-data.frame(gene=rownames(metabric), pval=pval, direction=direction)
rownames(surv_data_global_subtype)<-rownames(metabric)
nam <- paste("surv_data_global_", subtype, sep = "")
assign(nam, surv_data_global_subtype)
save(nam, file=paste("results/2025/", nam,".RData", sep=""))
}


#####5years
for(subtype in unique(meta$NOT_IN_OSLOVAL_Pam50Subtype)){
a=data.frame(ID=colnames(metabric)[meta$NOT_IN_OSLOVAL_Pam50Subtype==subtype],OS=surv5years[meta$NOT_IN_OSLOVAL_Pam50Subtype==subtype,1], IND=surv5years[meta$NOT_IN_OSLOVAL_Pam50Subtype==subtype,2])
pvalue_good<-c()
pvalue_poor<-c()

for(i in 1:nrow(metabric)){
  gene<-rownames(metabric)[i]
  b=data.frame(eigen_gene=metabric[gene,meta$NOT_IN_OSLOVAL_Pam50Subtype==subtype])
  c=cbind(a,b)
  
  eigen_class<-cut(metabric[gene,meta$NOT_IN_OSLOVAL_Pam50Subtype==subtype], breaks=quantile(metabric[gene,meta$NOT_IN_OSLOVAL_Pam50Subtype==subtype],probs=c(0,1,0.5), na.rm=T))
  
  pp<-survdiff(Surv(OS,IND)~eigen_class, data=c)
  pval<-1-pchisq(pp$chisq, length(pp$n)-1)
  
  if(pp$obs[1]>pp$exp[1]){
    pvalue_good<-c(pvalue_good, pval)
    pvalue_poor<-c(pvalue_poor, NA)
  } else {
    pvalue_good<-c(pvalue_good, NA)
    pvalue_poor<-c(pvalue_poor, pval)
  }
}

direction<-rep("None", nrow(metabric))
direction[!is.na(pvalue_good)]<-"Good"
direction[!is.na(pvalue_poor)]<-"Poor"

pval<-rep(NA, nrow(metabric))
pval[is.na(pvalue_good)==F]<-pvalue_good[is.na(pvalue_good)==F]
pval[is.na(pvalue_poor)==F]<-pvalue_poor[is.na(pvalue_poor)==F]

surv_data_subtype_5years<-data.frame(gene=rownames(metabric), pval=pval, direction=direction)
rownames(surv_data_subtype_5years)<-rownames(metabric)

nam <- paste("surv_data_", subtype,"_5years", sep = "")
assign(nam, surv_data_subtype_5years)
save(nam, file=paste("results/2025/", nam,".RData", sep=""))

save(surv_data_global_5years, file="results/2025/surv_data_global_5years.RData")
}


save(surv_data_global_5years, 
     surv_data_Basal_5years,
     surv_data_LumA_5years,
     surv_data_LumB_5years,
     surv_data_Her2_5years,
     surv_data_global_5years,
     surv_data_global_all,
     surv_data_global_Basal,
     surv_data_global_LumA,
     surv_data_global_LumB,
     surv_data_global_Her2, file="results/2025/surv_data.RData")
