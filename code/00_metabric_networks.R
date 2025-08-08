library(illuminaHumanv3.db)
library(WGCNA)

####################################
#####load metabric data
####################################

load(file="data/RData/metabric.RData")
load(file="data/RData/meta.RData")

####################################
## compute WGCNA networks
####################################


net_metabric = blockwiseModules(t(metabric),power= 6, corType="pearson", networkType = "signed", minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25,numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = F,
                                verbose = 3, maxBlockSize = 30000, nThreads=6)



save(net_metabric, "results/2025/net_metabric_oneblock.RData")

net_metabric_basal = blockwiseModules(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]),power= 6,TOMType ="signed", corType="pearson", networkType = "signed", minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25,numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = F,
                                      verbose = 3, maxBlockSize = 30000)

save(net_metabric_basal, "results/2025/net_metabric_Basal.RData")

####################################
##compute centrality
####################################
load(file="results/2025/net_metabric_oneblock.RData")
load(file="results/2025/net_metabric_Basal.RData")

moduleColors = labels2colors(net_metabric$colors)
moduleColors_basal = labels2colors(net_metabric_basal$colors)

ADJ=adjacency(t(metabric), power = 6, type= "signed")
Alldegrees=intramodularConnectivity(ADJ, moduleColors)

save(Alldegrees, file="results/2025/Alldegrees.RData")

ADJ_basal=adjacency(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]), power = 6, type= "signed")
Alldegrees_basal=intramodularConnectivity(ADJ_basal, moduleColors_basal)

save(Alldegrees_basal, file="results/2025/Alldegrees_basal.RData")


####################################

load(file="results/2025/Alldegrees.RData")
load(file="results/2025/Alldegrees_basal.RData")

centrality_global <- Alldegrees
centrality_basal <- Alldegrees_basal

################################################
#### Modules' enrichment
#################################################

library("msigdbr")
library(clusterProfiler)
library(openxlsx)

m_df = msigdbr(species = "Homo sapiens", category = "H")
m_t2g = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()


ids_background <- rownames(metabric)[moduleColors!="grey"]


ego<-list()
for(i in 1:length(unique(moduleColors))){
  ids <- rownames(metabric)[moduleColors==unique(moduleColors)[i]]
  ego[[i]]<-enricher(gene = ids, TERM2GENE = m_t2g, universe=ids_background, pvalueCutoff = 1)
  if(nrow(as.data.frame(ego[[i]]))>0){
    write.xlsx(as.data.frame(ego[[i]]), file=paste("results/2025/",unique(moduleColors)[i], "_msigdb_all.xlsx", sep=""))
  }
}

ids_background <- rownames(metabric)[moduleColors_basal!="grey"]


ego_basal<-list()
for(i in 1:length(unique(moduleColors_basal))){
  ids <- rownames(metabric)[moduleColors_basal==unique(moduleColors_basal)[i]]
  ego_basal[[i]]<-enricher(gene = ids, TERM2GENE = m_t2g, universe=ids_background, pvalueCutoff = 1)
  if(nrow(as.data.frame(ego_basal[[i]]))>0){
    write.xlsx(as.data.frame(ego_basal[[i]]), file=paste("results/",unique(moduleColors_basal)[i], "basal_msigdb_all.xlsx", sep=""))
  }
}

##########define alt_names based on ego and ego_basal

library(WGCNA)

alt_names<-c()
for(i in 1:length(unique(moduleColors))){
  alt_names<-c(alt_names, gsub("HALLMARK_", "", ego[[i]]$ID[1]))
}

names(alt_names)<-unique(moduleColors)
alt_names["grey"]<-"Unconnected"
alt_names[is.na(alt_names)]<-paste("NC", c(1:length(which(is.na(alt_names)))), sep="")
alt_names[which(duplicated(alt_names))]<-paste(alt_names[which(duplicated(alt_names))], "2", sep="")

alt_names_basal<-c()
for(i in 1:length(unique(moduleColors_basal))){
  alt_names_basal<-c(alt_names_basal, gsub("HALLMARK_", "", ego_basal[[i]]$ID[1]))
}

names(alt_names_basal)<-unique(moduleColors_basal)
alt_names_basal["grey"]<-"Unconnected"
alt_names_basal[which(duplicated(alt_names_basal))]<-paste(alt_names_basal[which(duplicated(alt_names_basal))], "2", sep="")
alt_names_basal<-paste("b_",alt_names_basal, sep="")
names(alt_names_basal)<-unique(moduleColors_basal)

#rename modules

moduleColors<-alt_names[moduleColors]
moduleColors_basal<-alt_names_basal[moduleColors_basal]

centrality_global$module<-moduleColors

rank<-rep(0, length(moduleColors))
for(i in 1:length(unique(moduleColors))){
  rank[moduleColors==unique(moduleColors)[i]]<-rank(-centrality_global$kWithin[moduleColors==unique(moduleColors)[i]])
  
}

centrality_global$rank<-rank

save(centrality_global, file="results/2025/centrality_global.RData")


centrality_basal$module<-moduleColors_basal

rank<-rep(0, length(moduleColors_basal))
for(i in 1:length(unique(moduleColors_basal))){
  rank[moduleColors_basal==unique(moduleColors_basal)[i]]<-rank(-centrality_basal$kWithin[moduleColors_basal==unique(moduleColors_basal)[i]])
  
}

centrality_basal$rank<-rank

centrality_basal[which(centrality_basal$module=="b_NA"),"module"]<-"b_NC1"
centrality_basal[which(centrality_basal$module=="b_NA2"),"module"]<-"b_NC2"


save(centrality_basal, file="results/2025/centrality_basal.RData")
