library(WGCNA)
load(file="data/RData/metabric.RData")
load(file="data/RData/meta.RData")
load(file="results/2025/centrality_basal.RData")

###moduleColors_basal in other subtypes
moduleColors_basal<-centrality_basal$module

ADJ1=adjacency(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Normal"]), power = 6, type= "signed")
Alldegrees_basalmod_normal=intramodularConnectivity(ADJ1, moduleColors_basal)


ADJ1=adjacency(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="LumA"]), power = 6, type= "signed")
Alldegrees_basalmod_lumA=intramodularConnectivity(ADJ1, moduleColors_basal)


ADJ1=adjacency(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="LumB"]), power = 6, type= "signed")
Alldegrees_basalmod_lumB=intramodularConnectivity(ADJ1, moduleColors_basal)


ADJ1=adjacency(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Her2"]), power = 6, type= "signed")
Alldegrees_basalmod_her2=intramodularConnectivity(ADJ1, moduleColors_basal)

ADJ1=adjacency(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]), power = 6, type= "signed")
Alldegrees_basalmod_basal=intramodularConnectivity(ADJ1, moduleColors_basal)

library(ggpubr)
for(i in 1:length(unique(moduleColors_basal))){
  
  df<-data.frame(kwithin=c(Alldegrees_basalmod_basal[which(moduleColors_basal==unique(moduleColors_basal)[i]),2], 
                           Alldegrees_basalmod_lumA[which(moduleColors_basal==unique(moduleColors_basal)[i]),2], 
                           Alldegrees_basalmod_lumB[which(moduleColors_basal==unique(moduleColors_basal)[i]),2], 
                           Alldegrees_basalmod_her2[which(moduleColors_basal==unique(moduleColors_basal)[i]),2], 
                           Alldegrees_basalmod_normal[which(moduleColors_basal==unique(moduleColors_basal)[i]),2]),
                 subtype=rep(c("Basal", "LumA", "LumB", "Her2", "Normal"), each=length(which(moduleColors_basal==unique(moduleColors_basal)[i]))))
  
  my_comparisons <- list( c("Basal", "LumA"), c("Basal", "LumB"), c("Basal", "Her2"), c("Basal", "Normal") )
  p <- ggboxplot(df, x = "subtype", y = "kwithin", title=unique(moduleColors_basal)[i])
  #  Add p-value
  p <- p + stat_compare_means(comparisons = my_comparisons, ref.group = "Basal", 
                              method = "wilcox.test", method.args = list(alternative = "greater")) # Add pairwise comparisons p-value
  
  
  png(paste("results/2025/",unique(moduleColors_basal)[i],"moduleColorsbasal_inothersubtypes_kwithin.png",sep=""), res=300, 1500, 1500)
 print(p)
 dev.off()
}

####module eigengenes across subtypes

MEs_basal= moduleEigengenes(t(metabric), centrality_basal$module)$eigengenes
colnames(MEs_basal)<-gsub("^ME","", colnames(MEs_basal))


for(i in 1:length(unique(moduleColors_basal))){
  
  df<-data.frame(ME=c(MEs_basal[meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal",unique(moduleColors_basal)[i]], 
                      MEs_basal[meta$NOT_IN_OSLOVAL_Pam50Subtype=="Her2",unique(moduleColors_basal)[i]], 
                      MEs_basal[meta$NOT_IN_OSLOVAL_Pam50Subtype=="LumA",unique(moduleColors_basal)[i]], 
                      MEs_basal[meta$NOT_IN_OSLOVAL_Pam50Subtype=="LumB",unique(moduleColors_basal)[i]], 
                      MEs_basal[meta$NOT_IN_OSLOVAL_Pam50Subtype=="Normal",unique(moduleColors_basal)[i]]),
                 subtype=rep(c("Basal", "Her2", "LumA", "LumB", "Normal"), times=table(meta$NOT_IN_OSLOVAL_Pam50Subtype)[c("Basal", "Her2", "LumA", "LumB", "Normal")]))
  
  my_comparisons <- list( c("Basal", "Her2"), c("Basal", "LumA"), c("Basal", "LumB"),  c("Basal", "Normal") )
  p <- ggboxplot(df, x = "subtype", y = "ME", title=unique(moduleColors_basal)[i])
  #  Add p-value
  p <- p + stat_compare_means(comparisons = my_comparisons, ref.group = "Basal", 
                              method = "wilcox.test", method.args = list(alternative = "greater")) # Add pairwise comparisons p-value
  
  
  png(paste("results/2025/",unique(moduleColors_basal)[i],"moduleColorsbasal_inothersubtypes_ME.png",sep=""), res=300, 1500, 1500)
  print(p)
  dev.off()
}

#######################
### in TCGA
########################
load("../../../Data/Data_downloaded/TCGA/TCGA-BRCA_Primary Tumor_dedupl.RData")
load("../../../Data/Data_downloaded/TCGA/TCGA-BRCA_clinical.RData")

library(TCGAbiolinks)
BRCA_subtype <- TCGAquery_subtype(tumor = "BRCA")

anno<-data.frame(barcode=colnames(TCGA_dedupl), 
                 subtype=BRCA_subtype$BRCA_Subtype_PAM50[match(substr(colnames(TCGA_dedupl),1,12), BRCA_subtype$patient)])

incommon<-intersect(rownames(metabric), rownames(TCGA_dedupl))
TCGA_dedupl_ic<-TCGA_dedupl[incommon,]
names(moduleColors_basal)<-rownames(metabric)
moduleColors_basal_ic<-moduleColors_basal[incommon]

ADJ1=adjacency(t(TCGA_dedupl_ic[,anno$subtype=="Normal"]), power = 6, type= "signed")
Alldegrees_basalmod_normal=intramodularConnectivity(ADJ1, moduleColors_basal_ic)


ADJ1=adjacency(t(TCGA_dedupl_ic[,anno$subtype=="LumA"]), power = 6, type= "signed")
Alldegrees_basalmod_lumA=intramodularConnectivity(ADJ1, moduleColors_basal_ic)


ADJ1=adjacency(t(TCGA_dedupl_ic[,anno$subtype=="LumB"]), power = 6, type= "signed")
Alldegrees_basalmod_lumB=intramodularConnectivity(ADJ1, moduleColors_basal_ic)


ADJ1=adjacency(t(TCGA_dedupl_ic[,anno$subtype=="Her2"]), power = 6, type= "signed")
Alldegrees_basalmod_her2=intramodularConnectivity(ADJ1, moduleColors_basal_ic)

ADJ1=adjacency(t(TCGA_dedupl_ic[,anno$subtype=="Basal"]), power = 6, type= "signed")
Alldegrees_basalmod_basal=intramodularConnectivity(ADJ1, moduleColors_basal_ic)


for(i in 1:length(unique(moduleColors_basal_ic))){
  
  df<-data.frame(kwithin=c(Alldegrees_basalmod_basal[which(moduleColors_basal_ic==unique(moduleColors_basal_ic)[i]),2], 
                           Alldegrees_basalmod_lumA[which(moduleColors_basal_ic==unique(moduleColors_basal_ic)[i]),2], 
                           Alldegrees_basalmod_lumB[which(moduleColors_basal_ic==unique(moduleColors_basal_ic)[i]),2], 
                           Alldegrees_basalmod_her2[which(moduleColors_basal_ic==unique(moduleColors_basal_ic)[i]),2], 
                           Alldegrees_basalmod_normal[which(moduleColors_basal_ic==unique(moduleColors_basal_ic)[i]),2]),
                 subtype=rep(c("Basal", "LumA", "LumB", "Her2", "Normal"), each=length(which(moduleColors_basal_ic==unique(moduleColors_basal_ic)[i]))))
  
  my_comparisons <- list( c("Basal", "LumA"), c("Basal", "LumB"), c("Basal", "Her2"), c("Basal", "Normal") )
  p <- ggboxplot(df, x = "subtype", y = "kwithin", title=unique(moduleColors_basal_ic)[i])
  #  Add p-value
  p <- p + stat_compare_means(comparisons = my_comparisons, ref.group = "Basal", 
                              method = "wilcox.test", method.args = list(alternative = "greater")) # Add pairwise comparisons p-value
  
  
  png(paste("results/2025/",unique(moduleColors_basal_ic)[i],"moduleColorsbasal_inothersubtypes_kwithin_TCGA.png",sep=""), res=300, 1500, 1500)
  print(p)
  dev.off()
}

