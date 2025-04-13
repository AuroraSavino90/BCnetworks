library(WGCNA)
load("../../../Data_downloaded/TCGA/TCGA-BRCA_Primary Tumor_dedupl.RData")
load("../../../Data_downloaded/TCGA/TCGA-BRCA_clinical.RData")

load(file="data/RData/metabric.RData")
load(file="data/RData/meta.RData")
load(file="data/RData/net_metabric_oneblock.RData")
load(file="data/RData/net_metabric_Basal.RData")
load(file="results/2025/centrality_global.RData")
load(file="results/2025/centrality_basal.RData")

library(TCGAbiolinks)
BRCA_subtype <- TCGAquery_subtype(tumor = "BRCA")

anno<-data.frame(barcode=colnames(TCGA_dedupl), 
                 subtype=BRCA_subtype$BRCA_Subtype_PAM50[match(substr(colnames(TCGA_dedupl),1,12), BRCA_subtype$patient)])
alldata_basal<-TCGA_dedupl[,which(anno$subtype=="Basal")]

setLabels = c("BC", "TCGA");
multiExpr = list(BC = list(data = t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"])), NKI = list(data = t(alldata_basal)));
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
save(mp, file = "results/2025/modulePreservation_metabricVsTCGA_basal.RData");



pdf("results/2025/preservation_metabricVsTCGA_basal.pdf")
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
text = rownames(mp$preservation$observed$ref.BC$inColumnsAlsoPresentIn.NKI);
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


##############################################

setLabels = c("BC", "TCGA");
multiExpr = list(BC = list(data = t(metabric)), NKI = list(data = t(TCGA_dedupl)));
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
save(mp, file = "results/2025/modulePreservation_metabricVsTCGA.RData");



pdf("results/2025/preservation_metabricVsTCGA.pdf")
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
text = rownames(mp$preservation$observed$ref.BC$inColumnsAlsoPresentIn.NKI);
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

