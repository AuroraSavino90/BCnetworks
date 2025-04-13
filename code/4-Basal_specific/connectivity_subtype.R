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