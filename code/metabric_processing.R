set.seed(575060)
library(illuminaHumanv3SYMBOL)

metabric<-read.csv("data/Complete_normalized_expression_data_METABRIC.txt", sep=",", dec=".")
load("data/Complete_METABRIC_Clinical_Features_Data.rbin")
meta<-Complete_METABRIC_Clinical_Features_Data

x <- illuminaHumanv3SYMBOL
# Get the probe identifiers that are mapped to a gene symbol
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
if(length(xx) > 0) {
  # Get the SYMBOL for the first five probes
  xx[1:5]
  # Get the first one
  xx[[1]]
}
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
adf <- data.frame(probe=names(vals), gene=vals)

# 5 Process and Export the data
# ned is our normalized expression mat
annotation_sel=adf[match( rownames(metabric), adf[,1]),2]

aned<-aned[which(is.na(annotation_sel)==F),]
annotation_sel<-na.omit(annotation_sel)
annotation_sel=as.character(annotation_sel)

for(i in 1:length(unique(annotation_sel))){
  if(length(which(annotation_sel==unique(annotation_sel)[i]))>1){
    m=which.max(rowMeans(metabric[which(annotation_sel==unique(annotation_sel)[i]),], na.rm=T))
    metabric=metabric[-which(annotation_sel==unique(annotation_sel)[i])[-m],]
    annotation_sel=annotation_sel[-which(annotation_sel==unique(annotation_sel)[i])[-m]]
  }
}

aned_conv=aned_conv[which(is.na(annotation_sel)==F),]
annotation_sel=na.omit(annotation_sel)
which(duplicated(annotation_sel))


rownames(metabric)=annotation_sel

net_metabric = blockwiseModules(t(metabric),power= 6,TOMType ="unsigned", corType="pearson", networkType = "signed", minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25,numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = T,
                               verbose = 3)

save(net_metabric, "results/2025/net_metabric.RData")

net_basal = blockwiseModules(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]),power= 6,TOMType ="unsigned", corType="pearson", networkType = "signed", minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25,numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = T,
                                verbose = 3)

save(net_basal, "results/2025/net_basal.RData")


moduleLabels_metabric = net_metabric$colors
moduleColors_metabric = labels2colors(net_metabric$colors)
MEs0_metabric = moduleEigengenes(t(metabric), moduleColors_metabric)$eigengenes
MEs_metabric = orderMEs(MEs0_metabric)
write.csv(data.frame(gene=rownames(metabric),module=moduleColors_NKI),"Modules_metabric.csv")


