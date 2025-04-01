#####Comparison with lines modules
library(WGCNA)
#change names in gene symbols
changenames<-function(data, anno){
  annotation_sel=anno[match( rownames(data), anno[,1]),2]
  
  if(length(which(annotation_sel==""))>0){
    data<-data[-which(annotation_sel==""),]
    annotation_sel<-annotation_sel[-which(annotation_sel=="")]
  }
  
  a<-which(duplicated(annotation_sel))
  while(length(a)>0){
    for(i in 1:length(unique(annotation_sel))){
      if(length(which(annotation_sel==unique(annotation_sel)[i]))>1){
        m=which.max(rowMeans(data[which(annotation_sel==unique(annotation_sel)[i]),], na.rm=T))
        data=data[-which(annotation_sel==unique(annotation_sel)[i])[-m],]
        annotation_sel=annotation_sel[-which(annotation_sel==unique(annotation_sel)[i])[-m]]
        print(i)
        print(dim(data))
      }
    }
    
    data=data[which(is.na(annotation_sel)==F),]
    annotation_sel=na.omit(annotation_sel)
    a<-which(duplicated(annotation_sel))
  }
  
  rownames(data)=annotation_sel
  return(data)
}

load(file="data/RData/metabric.RData")
load(file="data/RData/meta.RData")
load(file="data/RData/net_metabric_oneblock.RData")
load(file="data/RData/net_metabric_Basal.RData")
load(file="results/2025/centrality_global.RData")
load(file="results/2025/centrality_basal.RData")
meta_lines<-read.csv("data/GSE48213_series_matrix.csv", header=F)

setwd("data/GSE48213/")
files<-list.files()
x<-read.csv(files[1], row.names=1, header=T, sep="\t")
for(i in 2:length(files)){
  y<-read.csv(files[i], row.names=1, header=T, sep="\t")
  x<-cbind(x,y)
}

library(biomaRt)
ensembl <- useEnsembl(biomart = "ensembl", 
                      dataset = "hsapiens_gene_ensembl", 
                      mirror = "useast")

anno<-getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
      filters = 'ensembl_gene_id',
      values = rownames(x), 
      mart = ensembl)

sc_data<-changenames(x, anno)
basal_lines<-which(meta_lines[11,]=="subtype: Basal")
sc_data_basal<-sc_data[,basal_lines]

#gs<-goodSamplesGenes(t(x_basal))
#x_basal<-x_basal[gs$goodGenes,]
#gs<-goodSamplesGenes(metabric)

setLabels = c("BC", "lines");
multiExpr = list(BC = list(data = t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"])), lines = list(data = t(sc_data_basal)));
  multiColor = list(BC = centrality_basal$module )
system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3,
                          checkData=F)
} );
# Save the results
setwd("../../")
save(mp, file = "results/2025/modulePreservation_metabricVslines_basal.RData");

pdf("results/2025/preservation_metabricVslines_basal.pdf")
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

# Module labels and module sizes are also contained in the results
modColors = rep("grey", 20)
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("b_Unconnected", "gold"));
# Text labels for points
text = rownames(mp$preservation$observed$ref.BC$inColumnsAlsoPresentIn.lines);
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
p<-2
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  
}
# If plotting into a file, close it
dev.off()


###############All lines and global modules

setLabels = c("BC", "lines");
multiExpr = list(BC = list(data = t(metabric)), lines = list(data = t(sc_data)));
multiColor = list(BC = centrality_global$module )
system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3,
                          checkData=F)
} );
# Save the results
save(mp, file = "results/2025/modulePreservation_metabricVslines_global.RData");

pdf("results/2025/preservation_metabricVslines_basal.pdf")
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

# Module labels and module sizes are also contained in the results
modColors = rep("grey", 20)
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("Unconnected", "gold"));
# Text labels for points
text = rownames(mp$preservation$observed$ref.BC$inColumnsAlsoPresentIn.lines);
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
p<-2
min = min(plotData[, p], na.rm = TRUE);
max = max(plotData[, p], na.rm = TRUE);
# Adjust ploting ranges appropriately
if (p==2)
{
  if (min > -max/10) min = -max/10
  ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
} else
  ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
     main = mains[p],
     cex = 2.4,
     ylab = mains[p], xlab = "Module size", log = "x",
     ylim = ylim,
     xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
# For Zsummary, add threshold lines
if (p==2)
{
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
  
}
# If plotting into a file, close it
dev.off()



