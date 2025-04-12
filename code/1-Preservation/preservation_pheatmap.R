######################
## Load precomputed module preservations
#####################

load("results/2025/modulePreservation_metabricVsTCGA_basal.RData")
mp_TCGA_basal<-mp

load("results/2025/modulePreservation_metabricVsTCGA.RData")
mp_TCGA<-mp

load("results/2025/modulePreservation_metabricVsNKI_basal.RData")
mp_NKI_basal<-mp

load("results/2025/modulePreservation_metabricVsNKI_global.RData")
mp_NKI<-mp

load("results/2025/modulePreservation_metabricVslines_basal.RData")
mp_lines_basal<-mp

load("results/2025/modulePreservation_metabricVslines_global.RData")
mp_lines<-mp

load("results/2025/modulePreservation_metabricVsTRANSBIG.RData")
mp_transbig<-mp

load("results/2025/modulePreservation_metabricVsTRANSBIG_basal.RData")
mp_transbig_basal<-mp


load("results/2025/modulePreservation_metabricVsUNT.RData")
mp_unt<-mp

load("results/2025/modulePreservation_metabricVsUNT_basal.RData")
mp_unt_basal<-mp

load("results/2025/modulePreservation_metabricVsUPP.RData")
mp_upp<-mp

load("results/2025/modulePreservation_metabricVsUPP_basal.RData")
mp_upp_basal<-mp

load("results/2025/modulePreservation_metabricVsMAINZ.RData")
mp_mainz<-mp

load("results/2025/modulePreservation_metabricVsMAINZ_basal.RData")
mp_mainz_basal<-mp

load("results/2025/modulePreservation_metabricVsVDX.RData")
mp_vdx<-mp

load("results/2025/modulePreservation_metabricVsVDX_basal.RData")
mp_vdx_basal<-mp

ref = 1
test = 2

library(pheatmap)


##############################
#### Figure preservation global
############################

pheat<-cbind(mp_TCGA$preservation$Z[[ref]][[test]][, 2],
             mp_NKI$preservation$Z[[ref]][[test]][, 2],
             mp_transbig$preservation$Z[[ref]][[test]][, 2],
             mp_unt$preservation$Z[[ref]][[test]][, 2],
             mp_upp$preservation$Z[[ref]][[test]][, 2],
             mp_mainz$preservation$Z[[ref]][[test]][, 2],
             mp_vdx$preservation$Z[[ref]][[test]][, 2],
             mp_lines$preservation$Z[[ref]][[test]][, 2])
pheat[pheat>10]<-10
rownames(pheat)<-rownames(mp_NKI$preservation$Z[[ref]][[test]])
colnames(pheat)<-c("In TCGA", "In NKI", "In TRANSBIG", "In UNT", "In UPP", "In MAINZ", "In VDX","In BC cell lines")
#togliere moduli gold e Unconnected
pheat<-pheat[-which(rownames(pheat) %in% c("grey", "gold", "Unconnected")),]

##definire i colori per l'heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
myBreaks <- c(seq(min(unlist(log10(pheat)), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(unlist(log10(pheat)), na.rm=T)/paletteLength, max(unlist(log10(pheat)), na.rm=T), length.out=floor(paletteLength/2)))

png("results/2025/preservation_global_modules_alldatasets.png", res=300, 2500, 2000)
pheatmap(log10(pheat), cluster_cols = F,  cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor)
dev.off()

##############################
#### Figure preservation basal
############################

pheat<-cbind(mp_TCGA_basal$preservation$Z[[ref]][[test]][, 2],
             mp_NKI_basal$preservation$Z[[ref]][[test]][, 2],
             mp_transbig_basal$preservation$Z[[ref]][[test]][, 2],
             mp_unt_basal$preservation$Z[[ref]][[test]][, 2],
             mp_upp_basal$preservation$Z[[ref]][[test]][, 2],
             mp_mainz_basal$preservation$Z[[ref]][[test]][, 2],
             mp_vdx_basal$preservation$Z[[ref]][[test]][, 2],
             mp_lines_basal$preservation$Z[[ref]][[test]][, 2])
pheat[pheat>10]<-10
rownames(pheat)<-rownames(mp_NKI_basal$preservation$Z[[ref]][[test]])
colnames(pheat)<-c("In TCGA", "In NKI", "In TRANSBIG", "In UNT", "In UPP", "In MAINZ", "In VDX","In BC cell lines")
#togliere moduli gold e Unconnected
pheat<-pheat[-which(rownames(pheat) %in% c("grey", "gold", "b_Unconnected")),]

##definire i colori per l'heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
myBreaks <- c(seq(min(unlist(log10(pheat)), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(unlist(log10(pheat)), na.rm=T)/paletteLength, max(unlist(log10(pheat)), na.rm=T), length.out=floor(paletteLength/2)))

png("results/2025/preservation_basal_modules_alldatasets.png", res=300, 2500, 2000)
pheatmap(log10(pheat), cluster_cols = F,  cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor)
dev.off()
