##########survival each gene
###############
load("data/Complete_METABRIC_Clinical_Survival_Data__DSS.rbin")
load("data/RData/metabric.RData")
load("data/RData/meta.RData")
#check that samples are in the same order in expression and metadata
stopifnot(identical(rownames(meta), colnames(metabric)))

#load survival data, cut at 5years
surv<-Complete_METABRIC_Clinical_Survival_Data__DSS
surv5years<-surv
surv5years[surv[,1]>1825,"status"]<-0
surv5years[surv[,1]>1825,"time"]<-1825


#####survival######
#for each gene test the relationship with survival (log-rank, median split)
#output: signed statistic z = sign(O-E in high group)*sqrt(chisq), z>0: high expression = poor prognosis
#p-values computed with lower.tail=F to avoid underflow to 0

library(survival)
dir.create("results/2025_fix", showWarnings=F)

surv_each_gene<-function(expr, time, status){
  res<-t(apply(expr, 1, function(x){
    grp<-factor(x > median(x, na.rm=T), levels=c(FALSE, TRUE), labels=c("low", "high"))
    if(length(unique(na.omit(grp)))<2) return(c(z=NA, pval=NA))
    pp<-survdiff(Surv(time, status)~grp)
    z<-sign(pp$obs[2]-pp$exp[2])*sqrt(pp$chisq)
    c(z=z, pval=pchisq(pp$chisq, length(pp$n)-1, lower.tail=F))
  }))
  direction<-rep("None", nrow(expr))
  direction[which(res[,"z"]<0)]<-"Good"
  direction[which(res[,"z"]>0)]<-"Poor"
  out<-data.frame(gene=rownames(expr), z=res[,"z"], pval=res[,"pval"], direction=direction)
  rownames(out)<-rownames(expr)
  out
}

#####all follow-up and 5years, all patients
surv_data_global_all<-surv_each_gene(metabric, surv[,1], surv[,2])
save(surv_data_global_all, file="results/2025_fix/surv_data_global_all.RData")

surv_data_global_5years<-surv_each_gene(metabric, surv5years[,1], surv5years[,2])
save(surv_data_global_5years, file="results/2025_fix/surv_data_global_5years.RData")


#######Subtypes
subtypes<-setdiff(unique(meta$NOT_IN_OSLOVAL_Pam50Subtype), NA)
for(subtype in subtypes){
  idx<-which(meta$NOT_IN_OSLOVAL_Pam50Subtype==subtype)
  
  nam<-paste("surv_data_global_", subtype, sep="")
  assign(nam, surv_each_gene(metabric[,idx], surv[idx,1], surv[idx,2]))
  save(list=nam, file=paste("results/2025_fix/", nam,".RData", sep=""))
  
  nam<-paste("surv_data_", subtype,"_5years", sep="")
  assign(nam, surv_each_gene(metabric[,idx], surv5years[idx,1], surv5years[idx,2]))
  save(list=nam, file=paste("results/2025_fix/", nam,".RData", sep=""))
}


save(surv_data_global_5years,
     surv_data_Basal_5years,
     surv_data_LumA_5years,
     surv_data_LumB_5years,
     surv_data_Her2_5years,
     surv_data_global_all,
     surv_data_global_Basal,
     surv_data_global_LumA,
     surv_data_global_LumB,
     surv_data_global_Her2, file="results/2025_fix/surv_data.RData")