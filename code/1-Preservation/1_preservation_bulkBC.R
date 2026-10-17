library(breastCancerNKI)
library(Biobase)
library(WGCNA)
library(genefu)
load(file="data/RData/metabric.RData")
load(file="data/RData/meta.RData")
load(file="data/RData/net_metabric_oneblock.RData")
load(file="data/RData/net_metabric_Basal.RData")
load(file="results/2025/centrality_global.RData")
load(file="results/2025/centrality_basal.RData")

## separate output folder, so that the original results are not overwritten
outdir <- "results/2025_fix/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#change names in gene symbols  (unchanged)
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

################################
## Duplicated patients across TRANSBIG, UNT, UPP  (unchanged, from the genefu vignette)
################################
library("breastCancerMAINZ")
library("breastCancerTRANSBIG")
library("breastCancerUPP")
library("breastCancerUNT")
library("breastCancerVDX")
library(survcomp)
data(scmod2.robust)
data(pam50.robust)   # required by molecular.subtyping(sbt.model = "pam50"), as in the original script

data(breastCancerData)
cinfo <- colnames(pData(mainz7g))
data.all <- c("transbig7g"=transbig7g, "unt7g"=unt7g, "upp7g"=upp7g,
              "mainz7g"=mainz7g, "nki7g"=nki7g)
idtoremove.all <- NULL
duplres <- NULL
## No overlaps in the MainZ and NKI datasets.
## Focus on UNT vs UPP vs TRANSBIG
demo.all <- rbind(pData(transbig7g), pData(unt7g), pData(upp7g))
dn2 <- c("TRANSBIG", "UNT", "UPP")
## Karolinska
## Search for the VDXKIU, KIU, UPPU series
ds2 <- c("VDXKIU", "KIU", "UPPU")
demot <- demo.all[complete.cases(demo.all[ , c("series")]) & is.element(demo.all[ , "series"], ds2), ]
# Find the duplicated patients in that series
duplid <- sort(unique(demot[duplicated(demot[ , "id"]), "id"]))
duplrest <- NULL
for(i in 1:length(duplid)) {
  tt <- NULL
  for(k in 1:length(dn2)) {
    myx <- sort(row.names(demot)[complete.cases(demot[ , c("id", "dataset")]) &
                                   demot[ , "id"] == duplid[i] & demot[ , "dataset"] == dn2[k]])
    if(length(myx) > 0) { tt <- c(tt, myx) }
  }
  duplrest <- c(duplrest, list(tt))
}
names(duplrest) <- duplid
duplres <- c(duplres, duplrest)
## Oxford
## Search for the VVDXOXFU, OXFU series
ds2 <- c("VDXOXFU", "OXFU")
demot <- demo.all[complete.cases(demo.all[ , c("series")]) & is.element(demo.all[ , "series"], ds2), ]
# Find the duplicated patients in that series
duplid <- sort(unique(demot[duplicated(demot[ , "id"]), "id"]))
duplrest <- NULL
for(i in 1:length(duplid)) {
  tt <- NULL
  for(k in 1:length(dn2)) {
    myx <- sort(row.names(demot)[complete.cases(demot[ , c("id", "dataset")]) &
                                   demot[ , "id"] == duplid[i] & demot[ , "dataset"] == dn2[k]])
    if(length(myx) > 0) { tt <- c(tt, myx) }
  }
  duplrest <- c(duplrest, list(tt))
}
names(duplrest) <- duplid
duplres <- c(duplres, duplrest)
## Full set duplicated patients (the TRANSBIG copy is kept, the UNT/UPP copies are discarded)
duPL <- sort(unlist(lapply(duplres, function(x) { return(x[-1]) } )))


################################
## Same procedure for all datasets
################################
dn      <- c("transbig", "unt", "upp", "mainz", "nki", "vdx")
labels  <- c(transbig="TRANSBIG", unt="UNT", upp="UPP", mainz="MAINZ", nki="NKI", vdx="VDX")
sym_col <- c(transbig=3, unt=3, upp=3, mainz=3, nki=7, vdx=3)   # fData column with gene symbols (as in the original scripts)

## load a dataset, remove duplicated patients, assign PAM50, convert probes -> gene symbols
prepare_dataset <- function(name){
  dd <- get(data(list=name, envir=environment()))
  ddata  <- exprs(dd)
  dannot <- featureData(dd)@data
  ddata  <- ddata[, !colnames(ddata) %in% duPL]
  # PAM50 computed on the whole dataset (duplicates removed), at probe level
  PAM50Preds <- molecular.subtyping(sbt.model = "pam50", data=t(ddata),
                                    annot=dannot, do.mapping=TRUE)
  pam50 <- as.character(PAM50Preds$subtype[colnames(ddata)])
  message(name, ": gene symbol column = ", colnames(dannot)[sym_col[name]])
  print(table(pam50, useNA="ifany"))
  # probe -> gene symbol; for symbols with multiple probes the most expressed probe is kept (as in METABRIC)
  ddata <- changenames(data=ddata,
                       anno=data.frame(probe=rownames(ddata),
                                       symbol=as.character(dannot[rownames(ddata), sym_col[name]])))
  list(expr=ddata, pam50=pam50)
}

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

## module preservation with the same parameters for all datasets
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

metabric_basal <- metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]

for(name in dn){
  d <- prepare_dataset(name)
  # global modules, all samples (NKI file name as in the original, read by 4_preservation_pheatmap.R)
  stub_global <- ifelse(name == "nki", "NKI_global", labels[name])
  run_preservation(metabric, centrality_global$module, d$expr, labels[name], stub_global)
  # basal modules, PAM50 basal samples only
  run_preservation(metabric_basal, centrality_basal$module,
                   d$expr[, which(d$pam50 == "Basal")], labels[name], paste0(labels[name], "_basal"))
}
