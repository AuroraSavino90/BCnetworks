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

################################
##########NKI
###############################
data(nki)
NKI=exprs(nki)
annotation=fData(nki)
metadata=pData(nki)

NKI<-changenames(data=NKI, anno=cbind(annotation[,c(1,7)]))


########subtype prediction
library("breastCancerMAINZ")
library("breastCancerTRANSBIG")
library("breastCancerUPP")
library("breastCancerUNT")
library("breastCancerNKI")
library(survcomp)
data(scmod2.robust)
data(pam50.robust)

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
## Full set duplicated patients
duPL <- sort(unlist(lapply(duplres, function(x) { return(x[-1]) } )))




dn <- c("transbig", "unt", "upp", "mainz", "nki")
dn.platform <- c("affy", "affy", "affy", "affy", "agilent")
res <- ddemo.all <- ddemo.coln <- NULL
for(i in 1:length(dn)) {
  ## load dataset
  dd <- get(data(list=dn[i]))
  #Remove duplicates identified first
  message("obtained dataset!")
  #Extract expression set, pData, fData for each dataset
  ddata <- t(exprs(dd))
  ddemo <- phenoData(dd)@data
  if(length(intersect(rownames(ddata),duPL))>0)
  {
    ddata<-ddata[-which(rownames(ddata) %in% duPL),]
    ddemo<-ddemo[-which(rownames(ddemo) %in% duPL),]
  }
  dannot <- featureData(dd)@data
  # MOLECULAR SUBTYPING
  # Perform subtyping using scmod2.robust
  # scmod2.robust: List of parameters defining the subtype clustering model
  # (as defined by Wirapati et al)
  # OBSOLETE FUNCTION CALL - OLDER VERSIONS OF GENEFU
  # SubtypePredictions<-subtype.cluster.predict(sbt.model=scmod2.robust,data=ddata,
  # annot=dannot,do.mapping=TRUE,verbose=TRUE)
  # CURRENT FUNCTION CALL - NEWEST VERSION OF GENEFU
  SubtypePredictions<-molecular.subtyping(sbt.model = "scmod2",data = ddata,
                                          annot = dannot,do.mapping = TRUE)
  #Get sample counts pertaining to each subtype
  table(SubtypePredictions$subtype)
  #Select samples pertaining to Basal Subtype
  Basals<-names(which(SubtypePredictions$subtype == "ER-/HER2-"))
  #Select samples pertaining to HER2 Subtype
  HER2s<-names(which(SubtypePredictions$subtype == "HER2+"))
  #Select samples pertaining to Luminal Subtypes
  LuminalB<-names(which(SubtypePredictions$subtype == "ER+/HER2- High Prolif"))
  LuminalA<-names(which(SubtypePredictions$subtype == "ER+/HER2- Low Prolif"))
  #ASSIGN SUBTYPES TO EVERY SAMPLE, ADD TO THE EXISTING PHENODATA
  ddemo$SCMOD2<-SubtypePredictions$subtype
  ddemo[LuminalB,]$SCMOD2<-"LumB"
  ddemo[LuminalA,]$SCMOD2<-"LumA"
  ddemo[Basals,]$SCMOD2<-"Basal"
  ddemo[HER2s,]$SCMOD2<-"Her2"
  # Perform subtyping using PAM50
  # Matrix should have samples as ROWS, genes as COLUMNS
  # rownames(dannot)<-dannot$probe<-dannot$EntrezGene.ID
  # OLDER FUNCTION CALL
  # PAM50Preds<-intrinsic.cluster.predict(sbt.model=pam50,data=ddata,
  # annot=dannot,do.mapping=TRUE,verbose=TRUE)
  # NEWER FUNCTION CALL BASED ON MOST RECENT VERSION
  PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=ddata,
                                  annot=dannot,do.mapping=TRUE)
  table(PAM50Preds$subtype)
  ddemo$PAM50<-PAM50Preds$subtype
  LumA<-names(PAM50Preds$subtype)[which(PAM50Preds$subtype == "LumA")]
  LumB<-names(PAM50Preds$subtype)[which(PAM50Preds$subtype == "LumB")]
  ddemo[LumA,]$PAM50<-"LumA"
  ddemo[LumB,]$PAM50<-"LumB"
  ddemo.all <- rbind(ddemo, ddemo.all)
}
## obtained dat
##########

#gs<-goodSamplesGenes(t(NKI))
#NKI<-NKI[gs$goodGenes,]
#gs<-goodSamplesGenes(remaining_data)

setLabels = c("BC", "NKI");
multiExpr = list(BC = list(data = t(metabric)), NKI = list(data = t(NKI)));
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
save(mp, file = "results/2025/modulePreservation_metabricVsNKI_global.RData");



pdf("results/2025/preservation_metabricVsNKI.pdf")
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



#################preservation in basal
NKI_basal<-NKI[,ddemo.all[colnames(NKI),"PAM50"]=="Basal"]

setLabels = c("BC", "NKI");
multiExpr = list(BC = list(data = t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"])), NKI = list(data = t(NKI_basal)));
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
save(mp, file = "results/2025/modulePreservation_metabricVsNKI_basal.RData");



pdf("results/2025/preservation_metabricVsNKI_basal.pdf")
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


############preservation in other datasets
dn <- c("transbig", "unt", "upp", "mainz", "nki")
dn.platform <- c("affy", "affy", "affy", "affy", "agilent")

i<-1
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
ddata <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(ddata)<-dannot[rownames(ddata),3]

setLabels = c("BC", "TRANSBIG");
multiExpr = list(BC = list(data = t(metabric)), TRANSBIG = list(data = t(ddata)));
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
save(mp, file = "results/2025/modulePreservation_metabricVsTRANSBIG.RData");



i<-2
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
ddata <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(ddata)<-dannot[rownames(ddata),3]

setLabels = c("BC", "UNT");
multiExpr = list(BC = list(data = t(metabric)), UNT = list(data = t(ddata)));
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
save(mp, file = "results/2025/modulePreservation_metabricVsUNT.RData");





i<-3
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
ddata <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(ddata)<-dannot[rownames(ddata),3]

setLabels = c("BC", "UPP");
multiExpr = list(BC = list(data = t(metabric)), UPP = list(data = t(ddata)));
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
save(mp, file = "results/2025/modulePreservation_metabricVsUPP.RData");


i<-4
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
ddata <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(ddata)<-dannot[rownames(ddata),3]

setLabels = c("BC", "MAINZ");
multiExpr = list(BC = list(data = t(metabric)), MAINZ = list(data = t(ddata)));
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
save(mp, file = "results/2025/modulePreservation_metabricVsMAINZ.RData");


library("breastCancerVDX")


dn <- c("transbig", "unt", "upp", "mainz", "nki", "vdx")
dn.platform <- c("affy", "affy", "affy", "affy", "agilent", "affy")

i<-6

dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
ddata <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(ddata)<-dannot[rownames(ddata),3]

setLabels = c("BC", "VDX");
multiExpr = list(BC = list(data = t(metabric)), VDX = list(data = t(ddata)));
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
save(mp, file = "results/2025/modulePreservation_metabricVsVDX.RData");




########################################
######### Basal
######################################
metabric_basal<-metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]

##Basal
i<-1
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
ddata <- (exprs(dd))
dannot <- featureData(dd)@data
## load dataset
ddemo <- phenoData(dd)@data
# Perform subtyping using PAM50
# Matrix should have samples as ROWS, genes as COLUMNS
# rownames(dannot)<-dannot$probe<-dannot$EntrezGene.ID
# OLDER FUNCTION CALL
# PAM50Preds<-intrinsic.cluster.predict(sbt.model=pam50,data=ddata,
# annot=dannot,do.mapping=TRUE,verbose=TRUE)
# NEWER FUNCTION CALL BASED ON MOST RECENT VERSION
PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=t(ddata),
                                annot=dannot,do.mapping=TRUE)
table(PAM50Preds$subtype)
ddemo$PAM50<-PAM50Preds$subtype

rownames(ddata)<-dannot[rownames(ddata),3]

ddata<-ddata[,ddemo$PAM50=="Basal"]

setLabels = c("BC", "TRANSBIG");
multiExpr = list(BC = list(data = t(metabric_basal)), TRANSBIG = list(data = t(ddata)));
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
save(mp, file = "results/2025/modulePreservation_metabricVsTRANSBIG_basal.RData");


i<-2
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
ddata <- (exprs(dd))
dannot <- featureData(dd)@data
## load dataset
ddemo <- phenoData(dd)@data
# Perform subtyping using PAM50
# Matrix should have samples as ROWS, genes as COLUMNS
# rownames(dannot)<-dannot$probe<-dannot$EntrezGene.ID
# OLDER FUNCTION CALL
# PAM50Preds<-intrinsic.cluster.predict(sbt.model=pam50,data=ddata,
# annot=dannot,do.mapping=TRUE,verbose=TRUE)
# NEWER FUNCTION CALL BASED ON MOST RECENT VERSION
PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=t(ddata),
                                annot=dannot,do.mapping=TRUE)
table(PAM50Preds$subtype)
ddemo$PAM50<-PAM50Preds$subtype

rownames(ddata)<-dannot[rownames(ddata),3]

ddata<-ddata[,ddemo$PAM50=="Basal"]

setLabels = c("BC", "UNT");
multiExpr = list(BC = list(data = t(metabric_basal)), UNT = list(data = t(ddata)));
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
save(mp, file = "results/2025/modulePreservation_metabricVsUNT_basal.RData");



i<-3
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
ddata <- (exprs(dd))
dannot <- featureData(dd)@data
## load dataset
ddemo <- phenoData(dd)@data
# Perform subtyping using PAM50
# Matrix should have samples as ROWS, genes as COLUMNS
# rownames(dannot)<-dannot$probe<-dannot$EntrezGene.ID
# OLDER FUNCTION CALL
# PAM50Preds<-intrinsic.cluster.predict(sbt.model=pam50,data=ddata,
# annot=dannot,do.mapping=TRUE,verbose=TRUE)
# NEWER FUNCTION CALL BASED ON MOST RECENT VERSION
PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=t(ddata),
                                annot=dannot,do.mapping=TRUE)
table(PAM50Preds$subtype)
ddemo$PAM50<-PAM50Preds$subtype

rownames(ddata)<-dannot[rownames(ddata),3]

ddata<-ddata[,ddemo$PAM50=="Basal"]

setLabels = c("BC", "UPP");
multiExpr = list(BC = list(data = t(metabric_basal)), UPP = list(data = t(ddata)));
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
save(mp, file = "results/2025/modulePreservation_metabricVsUPP_basal.RData");


i<-4
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
ddata <- (exprs(dd))
dannot <- featureData(dd)@data
## load dataset
ddemo <- phenoData(dd)@data
# Perform subtyping using PAM50
# Matrix should have samples as ROWS, genes as COLUMNS
# rownames(dannot)<-dannot$probe<-dannot$EntrezGene.ID
# OLDER FUNCTION CALL
# PAM50Preds<-intrinsic.cluster.predict(sbt.model=pam50,data=ddata,
# annot=dannot,do.mapping=TRUE,verbose=TRUE)
# NEWER FUNCTION CALL BASED ON MOST RECENT VERSION
PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=t(ddata),
                                annot=dannot,do.mapping=TRUE)
table(PAM50Preds$subtype)
ddemo$PAM50<-PAM50Preds$subtype

rownames(ddata)<-dannot[rownames(ddata),3]

ddata<-ddata[,ddemo$PAM50=="Basal"]

setLabels = c("BC", "MAINZ");
multiExpr = list(BC = list(data = t(metabric_basal)), MAINZ = list(data = t(ddata)));
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
save(mp, file = "results/2025/modulePreservation_metabricVsMAINZ_basal.RData");


library("breastCancerVDX")


dn <- c("transbig", "unt", "upp", "mainz", "nki", "vdx")
dn.platform <- c("affy", "affy", "affy", "affy", "agilent", "affy")

i<-6

dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
ddata <- (exprs(dd))
dannot <- featureData(dd)@data
## load dataset
ddemo <- phenoData(dd)@data
# Perform subtyping using PAM50
# Matrix should have samples as ROWS, genes as COLUMNS
# rownames(dannot)<-dannot$probe<-dannot$EntrezGene.ID
# OLDER FUNCTION CALL
# PAM50Preds<-intrinsic.cluster.predict(sbt.model=pam50,data=ddata,
# annot=dannot,do.mapping=TRUE,verbose=TRUE)
# NEWER FUNCTION CALL BASED ON MOST RECENT VERSION
PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=t(ddata),
                                annot=dannot,do.mapping=TRUE)
table(PAM50Preds$subtype)
ddemo$PAM50<-PAM50Preds$subtype

rownames(ddata)<-dannot[rownames(ddata),3]

ddata<-ddata[,ddemo$PAM50=="Basal"]

setLabels = c("BC", "VDX");
multiExpr = list(BC = list(data = t(metabric_basal)), VDX = list(data = t(ddata)));
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
save(mp, file = "results/2025/modulePreservation_metabricVsVDX_basal.RData");



i<-5

dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
ddata <- (exprs(dd))
dannot <- featureData(dd)@data
## load dataset
ddemo <- phenoData(dd)@data
# Perform subtyping using PAM50
# Matrix should have samples as ROWS, genes as COLUMNS
# rownames(dannot)<-dannot$probe<-dannot$EntrezGene.ID
# OLDER FUNCTION CALL
# PAM50Preds<-intrinsic.cluster.predict(sbt.model=pam50,data=ddata,
# annot=dannot,do.mapping=TRUE,verbose=TRUE)
# NEWER FUNCTION CALL BASED ON MOST RECENT VERSION
PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=t(ddata),
                                annot=dannot,do.mapping=TRUE)
table(PAM50Preds$subtype)
ddemo$PAM50<-PAM50Preds$subtype

rownames(ddata)<-dannot[rownames(ddata),7]

ddata<-ddata[,ddemo$PAM50=="Basal"]

setLabels = c("BC", "NKI");
multiExpr = list(BC = list(data = t(metabric_basal)), NKI = list(data = t(ddata)));
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
save(mp, file = "results/2025/modulePreservation_metabricVsNKI_basal2.RData")

