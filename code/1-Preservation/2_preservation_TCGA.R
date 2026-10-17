library(WGCNA)
load("../../../Data/Data_downloaded/TCGA/TCGA-BRCA_Primary Tumor_dedupl.RData")
load("../../../Data/Data_downloaded/TCGA/TCGA-BRCA_clinical.RData")

load(file="data/RData/metabric.RData")
load(file="data/RData/meta.RData")
load(file="data/RData/net_metabric_oneblock.RData")
load(file="data/RData/net_metabric_Basal.RData")
load(file="results/2025/centrality_global.RData")
load(file="results/2025/centrality_basal.RData")

## separate output folder, so that the original results are not overwritten
outdir <- "results/2025_fix/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## PAM50 subtypes from TCGAquery_subtype, matched by patient barcode
## (TCGA_dedupl: log2(FPKM-UQ + 1), gene symbols, one sample per patient; see Methods)
library(TCGAbiolinks)
BRCA_subtype <- TCGAquery_subtype(tumor = "BRCA")

anno<-data.frame(barcode=colnames(TCGA_dedupl), 
                 subtype=BRCA_subtype$BRCA_Subtype_PAM50[match(substr(colnames(TCGA_dedupl),1,12), BRCA_subtype$patient)])
print(table(anno$subtype, useNA = "ifany"))
alldata_basal<-TCGA_dedupl[,which(anno$subtype=="Basal")]

################################
## Module preservation (same parameters as in 1_preservation_bulkBC.R and 3_preservation_cell_lines.R)
################################
## diagnostic plot of Zsummary vs module size (gold and Unconnected excluded)
plot_preservation <- function(mp, file){
  ref = 1
  test = 2
  statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
  statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])
  print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
               signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
  Z <- mp$preservation$Z[[ref]][[test]]
  plotMods <- !(rownames(Z) %in% c("gold", "Unconnected", "b_Unconnected"))
  moduleSizes <- Z[plotMods, "moduleSize"]
  zs <- Z[plotMods, "Zsummary.pres"]
  min = min(zs, na.rm = TRUE)
  max = max(zs, na.rm = TRUE)
  if (min > -max/10) min = -max/10
  ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  pdf(file)
  plot(moduleSizes, zs, col = 1, bg = "grey", pch = 21,
       main = "Preservation Zsummary", cex = 2.4,
       ylab = "Preservation Zsummary", xlab = "Module size", log = "x",
       ylim = ylim, cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.4)
  labelPoints(moduleSizes, zs, rownames(Z)[plotMods], cex = 1, offs = 0.08)
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
  dev.off()
}

run_preservation <- function(ref_expr, ref_modules, test_expr, test_label, file_stub){
  # remove genes/samples with too many missing values or zero variance in each set (modulePreservation does not remove them itself)
  gsg_ref  <- goodSamplesGenes(t(ref_expr),  verbose = 0)
  gsg_test <- goodSamplesGenes(t(test_expr), verbose = 0)
  message(file_stub, ": removed ", sum(!gsg_ref$goodGenes), " genes / ", sum(!gsg_ref$goodSamples), " samples (reference), ",
          sum(!gsg_test$goodGenes), " genes / ", sum(!gsg_test$goodSamples), " samples (test)")
  ref_expr    <- ref_expr[gsg_ref$goodGenes, gsg_ref$goodSamples]
  ref_modules <- ref_modules[gsg_ref$goodGenes]      # module colors kept aligned with the retained genes
  test_expr   <- test_expr[gsg_test$goodGenes, gsg_test$goodSamples]
  # no subsampling of modules: maxModuleSize = size of the largest module (Unconnected excluded)
  maxSize <- max(table(ref_modules[!ref_modules %in% c("Unconnected", "b_Unconnected")]))
  multiExpr = list(BC = list(data = t(ref_expr)), test = list(data = t(test_expr)))
  names(multiExpr)[2] <- test_label
  multiColor = list(BC = ref_modules)
  print(system.time( {
    mp = modulePreservation(multiExpr, multiColor,
                            referenceNetworks = 1,
                            networkType = "signed",     # as in the reference networks
                            maxModuleSize = maxSize,
                            nPermutations = 200,
                            randomSeed = 1,
                            quickCor = 0,
                            verbose = 3,
                            checkData = TRUE)           # check only: problematic genes/samples have already been removed
  } ))
  save(mp, file = paste0(outdir, "modulePreservation_metabricVs", file_stub, ".RData"))
  plot_preservation(mp, paste0(outdir, "preservation_metabricVs", file_stub, ".pdf"))
  invisible(mp)
}

## basal modules, PAM50 basal samples only
run_preservation(metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], centrality_basal$module,
                 alldata_basal, "TCGA", "TCGA_basal")

## global modules, all samples (file name as in the original, read by 4_preservation_pheatmap.R)
run_preservation(metabric, centrality_global$module, TCGA_dedupl, "TCGA", "TCGA")
