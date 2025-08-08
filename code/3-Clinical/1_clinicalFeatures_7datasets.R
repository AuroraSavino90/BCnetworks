library(breastCancerNKI)
library(Biobase)
library(WGCNA)
library(genefu)
library(TCGAbiolinks)
library(pheatmap)
library(survcomp)
data(scmod2.robust)
data(pam50.robust)

#############load network objects and data

load("data/RData/metabric.RData")
load("data/RData/meta.RData")
load("data/RData/net_metabric_Basal.RData")
load("results/2025/centrality_basal.RData")
load("results/2025/centrality_global.RData")
load("data/Complete_METABRIC_Clinical_Survival_Data__DSS.rbin")

library("breastCancerMAINZ")
library("breastCancerTRANSBIG")
library("breastCancerUPP")
library("breastCancerUNT")
library("breastCancerNKI")
library("breastCancerVDX")


dn <- c("transbig", "unt", "upp", "mainz", "nki", "vdx")
dn.platform <- c("affy", "affy", "affy", "affy", "agilent", "affy")
feat<-c("grade", "age")

#############compute module eigengenes

MEs= moduleEigengenes(t(metabric), centrality_global$module)$eigengenes
colnames(MEs)<-gsub("^ME","", colnames(MEs))
MEs_basal= moduleEigengenes(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]), centrality_basal$module)$eigengenes
colnames(MEs_basal)<-gsub("^ME","", colnames(MEs_basal))

#######################
### FUNCTION: module vs clinical feature
######################

moduletrait<-function(features, metadata, original_dataset=metabric, new_dataset,genes_annot, centrality, MEs, new_dataset_name="TRANSBIG", basal=F){
  
  if(basal){
  PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=t(new_dataset),
                                  annot=genes_annot)
  PAM50<-PAM50Preds$subtype
  new_dataset<-new_dataset[,PAM50=="Basal"]
  metadata<-metadata[PAM50=="Basal",]
  }
  
  new_dataset<-new_dataset[which(rowSums(is.na(new_dataset))==0),]
  
  # MEs approximation for new_dataset
  # perform principal components analysis
  # project new data onto the PCA space
  
  pcaproj<-matrix(nrow=ncol(new_dataset), ncol=length(unique(centrality$module)))
  for(i in 1:length(unique(centrality$module))){
    pca <- prcomp(t(original_dataset[centrality$module==unique(centrality$module)[i],]))
    if(cor(pca$x[,1], MEs[,unique(centrality$module)[i]])<0){
      pca$rotation[,1]<-(-pca$rotation[,1])
    }
    commongenes<-rownames(new_dataset)[which(rownames(new_dataset) %in% rownames(original_dataset[centrality$module==unique(centrality$module)[i],]))]
    data_forpca<-new_dataset[commongenes,]
    pcaproj[,i]<-colSums(t(scale(t(data_forpca), pca$center[commongenes], pca$scale)) * c(pca$rotation[commongenes,1]))
  }

  colnames(pcaproj)<-unique(centrality$module)
  
  
  ##Correlation between the projected MEs and clinical features
  
  moduleTraitCor = list();
moduleTraitPvalue = list();
# Calculate the correlations
for (set in features)
{
  moduleTraitCor[[set]] = cor(pcaproj, as.numeric(as.character(metadata[,set])), use = "p");
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], nrow(pcaproj));
}

moduleTraitCor_tot<-moduleTraitCor[[1]]
moduleTraitPvalue_tot<-moduleTraitPvalue[[1]]

for (set in 2:length(features)){
  moduleTraitCor_tot<-cbind(moduleTraitCor_tot,moduleTraitCor[[set]])
  moduleTraitPvalue_tot<-cbind(moduleTraitPvalue_tot, moduleTraitPvalue[[set]])
}
colnames(moduleTraitCor_tot)<-features
colnames(moduleTraitPvalue_tot)<-features

if(basal){
  save(moduleTraitCor_tot, file=paste("results/2025/moduleTraitCor_",new_dataset_name,"_basal.RData", sep=""))
  save(moduleTraitPvalue_tot, file=paste("results/2025/moduleTraitPvalue_", new_dataset_name, "_basal.RData", sep=""))
  
} else {
save(moduleTraitCor_tot, file=paste("results/2025/moduleTraitCor_",new_dataset_name,".RData", sep=""))
save(moduleTraitPvalue_tot, file=paste("results/2025/moduleTraitPvalue_", new_dataset_name, ".RData", sep=""))
}
}

####################################
############ FUNCTION: MODULE vs SURVIVAL
#####################################

modulesurv<-function(features, surv, original_dataset=metabric, new_dataset,genes_annot, centrality, MEs, basal=F){
  
  if(basal){
    PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=t(new_dataset),
                                    annot=genes_annot)
    PAM50<-PAM50Preds$subtype
    new_dataset<-new_dataset[,PAM50=="Basal"]
    metadata<-metadata[PAM50=="Basal",]
    surv<-surv[PAM50=="Basal",]
  }
  
  new_dataset<-new_dataset[which(rowSums(is.na(new_dataset))==0),]
  
  # MEs approximation for new_dataset
  # perform principal components analysis
  # project new data onto the PCA space
  
  pcaproj<-matrix(nrow=ncol(new_dataset), ncol=length(unique(centrality$module)))
  for(i in 1:length(unique(centrality$module))){
    pca <- prcomp(t(original_dataset[centrality$module==unique(centrality$module)[i],]))
    if(cor(pca$x[,1], MEs[,unique(centrality$module)[i]])<0){
      pca$rotation[,1]<-(-pca$rotation[,1])
    }
    commongenes<-rownames(new_dataset)[which(rownames(new_dataset) %in% rownames(original_dataset[centrality$module==unique(centrality$module)[i],]))]
    data_forpca<-new_dataset[commongenes,]
    pcaproj[,i]<-colSums(t(scale(t(data_forpca), pca$center[commongenes], pca$scale)) * c(pca$rotation[commongenes,1]))
  }
  
  colnames(pcaproj)<-unique(centrality$module)
  
  ################# relationship between projected ME and survival
  
  
  surv[which(surv[,1]>1825),"status"]<-0
  surv[which(surv[,1]>1825),"time"]<-1825
  
  pvalue_good<-c()
  pvalue_poor<-c()
  
  for(module in unique(centrality$module)){
    a=data.frame(ID=colnames(new_dataset),OS=surv$time, IND=surv$status)
    b=data.frame(eigen_gene=pcaproj[,module])
    c=cbind(a,b)
    eigen_class<-cut(pcaproj[,module], breaks=quantile(pcaproj[,module],probs=c(0,0.5,1), na.rm=T))
    
    
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
  
  names(pvalue_poor)<-unique(centrality$module)
  names(pvalue_good)<-unique(centrality$module)
  
  return(list(pvalue_poor, pvalue_good))
  
}

#################################
####### All datasets
#################################

pvalue_poor_tot<-matrix(nrow=length(unique(centrality_global$module)), ncol=8)
pvalue_good_tot<-matrix(nrow=length(unique(centrality_global$module)), ncol=8)
pvalue_poor_tot_b<-matrix(nrow=length(unique(centrality_basal$module)), ncol=8)
pvalue_good_tot_b<-matrix(nrow=length(unique(centrality_basal$module)), ncol=8)

#######Correlation with grade

i<-1
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
data <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(data)<-dannot[rownames(data),3]

metadata<-pData(dd)
surv<-data.frame(status=metadata$e.rfs, time=metadata$t.rfs)

moduletrait(features=feat, metadata=metadata, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
            centrality=centrality_global, MEs=MEs, new_dataset_name="TRANSBIG", basal=F)
moduletrait(features=feat, metadata=metadata, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
            centrality=centrality_basal, MEs=MEs_basal, new_dataset_name="TRANSBIG", basal=T)


p<-modulesurv(features=feat, surv=surv, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
              centrality=centrality_global, MEs=MEs, basal=F)

pvalue_poor_tot[,2]<-p[[1]]
pvalue_good_tot[,2]<-p[[2]]

p<-modulesurv(features=feat, surv=surv, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
              centrality=centrality_basal, MEs=MEs_basal, basal=T)

pvalue_poor_tot_b[,2]<-p[[1]]
pvalue_good_tot_b[,2]<-p[[2]]


i<-2
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
data <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(data)<-dannot[rownames(data),3]

metadata<-pData(dd)
surv<-data.frame(status=metadata$e.rfs, time=metadata$t.rfs)


moduletrait(features=feat, metadata=metadata, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
            centrality=centrality_global, MEs=MEs, new_dataset_name="UNT", basal=F)
moduletrait(features=feat, metadata=metadata, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
            centrality=centrality_basal, MEs=MEs_basal, new_dataset_name="UNT", basal=T)
p<-modulesurv(features=feat, surv=surv, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
              centrality=centrality_global, MEs=MEs, basal=F)

pvalue_poor_tot[,3]<-p[[1]]
pvalue_good_tot[,3]<-p[[2]]

p<-modulesurv(features=feat, surv=surv, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
              centrality=centrality_basal, MEs=MEs_basal, basal=T)
pvalue_poor_tot_b[,3]<-p[[1]]
pvalue_good_tot_b[,3]<-p[[2]]

i<-3
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
data <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(data)<-dannot[rownames(data),3]

metadata<-pData(dd)
surv<-data.frame(status=metadata$e.rfs, time=metadata$t.rfs)


moduletrait(features=feat, metadata=metadata, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
            centrality=centrality_global, MEs=MEs, new_dataset_name="UPP", basal=F)
moduletrait(features=feat, metadata=metadata, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
            centrality=centrality_basal, MEs=MEs_basal, new_dataset_name="UPP", basal=T)
p<-modulesurv(features=feat, surv=surv, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
              centrality=centrality_global, MEs=MEs, basal=F)

pvalue_poor_tot[,4]<-p[[1]]
pvalue_good_tot[,4]<-p[[2]]

p<-modulesurv(features=feat, surv=surv, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
              centrality=centrality_basal, MEs=MEs_basal, basal=T)
pvalue_poor_tot_b[,4]<-p[[1]]
pvalue_good_tot_b[,4]<-p[[2]]


i<-4
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
data <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(data)<-dannot[rownames(data),3]

metadata<-pData(dd)
surv<-data.frame(status=metadata$e.dmfs, time=metadata$t.dmfs)

moduletrait(features=feat, metadata=metadata, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
            centrality=centrality_global, MEs=MEs, new_dataset_name="MAINZ", basal=F)
moduletrait(features=feat, metadata=metadata, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
            centrality=centrality_basal, MEs=MEs_basal, new_dataset_name="MAINZ", basal=T)



p<-modulesurv(features=feat, surv=surv, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
              centrality=centrality_global, MEs=MEs, basal=F)

pvalue_poor_tot[,5]<-p[[1]]
pvalue_good_tot[,5]<-p[[2]]

p<-modulesurv(features=feat, surv=surv, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
              centrality=centrality_basal, MEs=MEs_basal, basal=T)
pvalue_poor_tot_b[,5]<-p[[1]]
pvalue_good_tot_b[,5]<-p[[2]]


i<-5
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
data <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(data)<-dannot[rownames(data),6]

metadata<-pData(dd)
surv<-data.frame(status=metadata$e.rfs, time=metadata$t.rfs)

moduletrait(features=feat, metadata=metadata, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
            centrality=centrality_global, MEs=MEs, new_dataset_name="NKI", basal=F)
moduletrait(features=feat, metadata=metadata, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
            centrality=centrality_basal, MEs=MEs_basal, new_dataset_name="NKI", basal=T)

p<-modulesurv(features=feat, surv=surv, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
              centrality=centrality_global, MEs=MEs, basal=F)

pvalue_poor_tot[,1]<-p[[1]]
pvalue_good_tot[,1]<-p[[2]]

p<-modulesurv(features=feat, surv=surv, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
              centrality=centrality_basal, MEs=MEs_basal, basal=T)
pvalue_poor_tot_b[,1]<-p[[1]]
pvalue_good_tot_b[,1]<-p[[2]]



i<-6
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
data <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(data)<-dannot[rownames(data),3]

metadata<-pData(dd)
surv<-data.frame(status=metadata$e.dmfs, time=metadata$t.dmfs)

moduletrait(features=feat, metadata=metadata, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
            centrality=centrality_global, MEs=MEs, new_dataset_name="VDX", basal=F)
moduletrait(features=feat, metadata=metadata, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
            centrality=centrality_basal, MEs=MEs_basal, new_dataset_name="VDX", basal=T)


p<-modulesurv(features=feat, surv=surv, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
              centrality=centrality_global, MEs=MEs, basal=F)

pvalue_poor_tot[,6]<-p[[1]]
pvalue_good_tot[,6]<-p[[2]]

p<-modulesurv(features=feat, surv=surv, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
              centrality=centrality_basal, MEs=MEs_basal, basal=T)
pvalue_poor_tot_b[,6]<-p[[1]]
pvalue_good_tot_b[,6]<-p[[2]]

#TCGA

load("../../../Data/Data_downloaded/TCGA/TCGA-BRCA_Primary Tumor_dedupl.RData")
load("../../../Data/Data_downloaded/TCGA/TCGA-BRCA_clinical.RData")

library(TCGAbiolinks)
BRCA_subtype <- TCGAquery_subtype(tumor = "BRCA")

anno<-data.frame(barcode=colnames(TCGA_dedupl), 
                 subtype=BRCA_subtype$BRCA_Subtype_PAM50[match(substr(colnames(TCGA_dedupl),1,12), BRCA_subtype$patient)])
alldata_basal<-TCGA_dedupl[,which(anno$subtype=="Basal")]

surv<-data.frame(status=clinical.BCRtab.all$clinical_follow_up_v4.0_brca$vital_status[match(substr(colnames(TCGA_dedupl),1,12), clinical.BCRtab.all$clinical_follow_up_v4.0_brca$bcr_patient_barcode)],
                 timed=clinical.BCRtab.all$clinical_follow_up_v4.0_brca$death_days_to[match(substr(colnames(TCGA_dedupl),1,12), clinical.BCRtab.all$clinical_follow_up_v4.0_brca$bcr_patient_barcode)],
                 timef=clinical.BCRtab.all$clinical_follow_up_v4.0_brca$last_contact_days_to[match(substr(colnames(TCGA_dedupl),1,12), clinical.BCRtab.all$clinical_follow_up_v4.0_brca$bcr_patient_barcode)])
surv$status[surv$status=="Alive"]<-0
surv$status[surv$status=="Dead"]<-1
surv$time<-as.numeric(surv$timed)
surv$time[which(surv$status==0)]<-as.numeric(surv$timef[which(surv$status==0)])
surv$status<-as.numeric(surv$status)

p<-modulesurv(features=feat, surv=surv[,c("status", "time")], original_dataset=metabric, genes_annot = dannot, new_dataset=TCGA_dedupl, 
              centrality=centrality_global, MEs=MEs, basal=F)


pvalue_poor_tot[,7]<-p[[1]]
pvalue_good_tot[,7]<-p[[2]]

p<-modulesurv(features=feat, surv=surv[which(anno$subtype=="Basal"),c("status", "time")], original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = NA, new_dataset=TCGA_dedupl[,which(anno$subtype=="Basal")], 
              centrality=centrality_basal, MEs=MEs_basal, basal=F)


pvalue_poor_tot_b[,7]<-p[[1]]
pvalue_good_tot_b[,7]<-p[[2]]


#################################
#########METABRIC global
#################################

###moduletrait
feat<-c("grade", "age_at_diagnosis")

moduleTraitCor = list();
moduleTraitPvalue = list();
# Calculate the correlations
for (set in feat){
  moduleTraitCor[[set]] = cor(MEs, as.numeric(as.character(meta[,set])), use = "p");
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], nrow(MEs));
}

moduleTraitCor_tot<-moduleTraitCor[[1]]
moduleTraitPvalue_tot<-moduleTraitPvalue[[1]]

for (set in 2:length(feat)){
  moduleTraitCor_tot<-cbind(moduleTraitCor_tot,moduleTraitCor[[set]])
  moduleTraitPvalue_tot<-cbind(moduleTraitPvalue_tot, moduleTraitPvalue[[set]])
}
colnames(moduleTraitCor_tot)<-feat
colnames(moduleTraitPvalue_tot)<-feat

save(moduleTraitCor_tot, file="results/2025/moduleTraitCor_METABRIC.RData")
save(moduleTraitPvalue_tot, file="results/2025/moduleTraitPvalue_METABRIC.RData")

##########survival

pvalue_good<-c()
pvalue_poor<-c()

surv<-Complete_METABRIC_Clinical_Survival_Data__DSS
surv[surv[,1]>1825,"status"]<-0
surv[surv[,1]>1825,"time"]<-1825


for(module in unique(centrality_global$module)){
  a=data.frame(ID=colnames(metabric),OS=surv[,1], IND=surv[,2])
  b=data.frame(ID=rownames(MEs), eigen_gene=MEs[,module])
  c=cbind(a,b)
  
  eigen_class<-cut(MEs[,module], breaks=quantile(MEs[,module],probs=c(0,1,0.5), na.rm=T))
  
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

names(pvalue_poor)<-unique(centrality_global$module)
names(pvalue_good)<-unique(centrality_global$module)

pvalue_poor_tot[,8]<-pvalue_poor
pvalue_good_tot[,8]<-pvalue_good

#################################
#########METABRIC basal
#################################

###moduletrait
moduleTraitCor = list();
moduleTraitPvalue = list();
# Calculate the correlations
for (set in feat){
  moduleTraitCor[[set]] = cor(MEs_basal, as.numeric(as.character(meta[meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal",set])), use = "p");
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], nrow(MEs_basal));
}

moduleTraitCor_tot<-moduleTraitCor[[1]]
moduleTraitPvalue_tot<-moduleTraitPvalue[[1]]

for (set in 2:length(feat)){
  moduleTraitCor_tot<-cbind(moduleTraitCor_tot,moduleTraitCor[[set]])
  moduleTraitPvalue_tot<-cbind(moduleTraitPvalue_tot, moduleTraitPvalue[[set]])
}
colnames(moduleTraitCor_tot)<-feat
colnames(moduleTraitPvalue_tot)<-feat

save(moduleTraitCor_tot, file="results/2025/moduleTraitCor_METABRIC_basal.RData")
save(moduleTraitPvalue_tot, file="results/2025/moduleTraitPvalue_METABRIC_basal.RData")


############survival

pvalue_good<-c()
pvalue_poor<-c()

surv<-Complete_METABRIC_Clinical_Survival_Data__DSS
surv[surv[,1]>1825,"status"]<-0
surv[surv[,1]>1825,"time"]<-1825


for(module in unique(centrality_basal$module)){
  a=data.frame(ID=colnames(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]),OS=surv[meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal",1], IND=surv[meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal",2])
  b=data.frame(ID=rownames(MEs_basal), eigen_gene=MEs_basal[,module])
  c=cbind(a,b)
  
  eigen_class<-cut(MEs_basal[,module], breaks=quantile(MEs_basal[,module],probs=c(0,1,0.5), na.rm=T))
  
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

names(pvalue_poor)<-unique(centrality_basal$module)
names(pvalue_good)<-unique(centrality_basal$module)

pvalue_poor_tot_b[,8]<-pvalue_poor
pvalue_good_tot_b[,8]<-pvalue_good



###########################Overall representation of survival data
##########################Diversi tipi di survival nei diversi dataset -> DSS, DMFS, RFS

pvalue_good_tot[pvalue_good_tot<2.2*10^(-16)]<-2.2*10^(-16)
pvalue_poor_tot[pvalue_poor_tot<2.2*10^(-16)]<-2.2*10^(-16)

pvalue_merge<-(-log10(pvalue_good_tot))
pvalue_merge[is.na(pvalue_poor_tot)==F]<-(log10(pvalue_poor_tot)[is.na(pvalue_poor_tot)==F])

rownames(pvalue_merge)<-unique(centrality_global$module)
colnames(pvalue_merge)<-c("NKI", "TRANSBIG", "UNT", "UPP", "MAINZ", "VDX", "TCGA","METABRIC")

save(pvalue_merge, file="results/2025/pvalue_merge_surv_global.RData")

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths


png("results/2025/Suvr_7datasets_5years.png", res=500, 4000, 3000)
pheatmap(pvalue_merge[-which(rownames(pvalue_merge)=="Unconnected"),], cellwidth = 15, cellheight = 15, color = myColor)
dev.off()
pdf("results/2025/Suvr_7datasets_5years.pdf", 20, 15)
pheatmap(pvalue_merge[-which(rownames(pvalue_merge)=="Unconnected"),], cellwidth = 15, cellheight = 15, color = myColor)
dev.off()

pvalue_good_tot_b[pvalue_good_tot_b<2.2*10^(-16)]<-2.2*10^(-16)
pvalue_poor_tot_b[pvalue_poor_tot_b<2.2*10^(-16)]<-2.2*10^(-16)

pvalue_merge_b<-(-log10(pvalue_good_tot_b))
pvalue_merge_b[is.na(pvalue_poor_tot_b)==F]<-(log10(pvalue_poor_tot_b)[is.na(pvalue_poor_tot_b)==F])

rownames(pvalue_merge_b)<-unique(centrality_basal$module)
colnames(pvalue_merge_b)<-c("NKI", "TRANSBIG", "UNT", "UPP", "MAINZ", "VDX", "TCGA","METABRIC")

save(pvalue_merge_b, file="results/2025/pvalue_merge_surv_basal.RData")

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths


png("results/2025/Suvr_7datasets_5years_b.png", res=500, 4000, 3000)
pheatmap(pvalue_merge_b[-which(rownames(pvalue_merge_b)=="b_Unconnected"),], cellwidth = 15, cellheight = 15, color = myColor)
dev.off()
pdf("results/2025/Suvr_7datasets_5years_b.pdf",  20, 15)
pheatmap(pvalue_merge_b[-which(rownames(pvalue_merge_b)=="b_Unconnected"),], cellwidth = 15, cellheight = 15, color = myColor)
dev.off()



