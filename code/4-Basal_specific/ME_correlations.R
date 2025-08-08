#############load network objects and data
load("data/RData/metabric.RData")
load("data/RData/meta.RData")
load("results/2025/centrality_basal.RData")

#############compute module eigengenes
MEs_basal= moduleEigengenes(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]), centrality_basal$module)$eigengenes
colnames(MEs_basal)<-gsub("^ME","", colnames(MEs_basal))

library(pheatmap)
library(ggplot2)

MEs_basal<-MEs_basal[,-which(colnames(MEs_basal)=="b_Unconnected")]
toplot<-cor(MEs_basal)

paletteLength <- 50
# use floor and ceiling to deal with even/odd length pallettelengths
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

png("results/2025/Heatmap_b_cormod.png", res=300, 3000, 3000)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor)
dev.off()

pdf("results/2025/Heatmap_b_cormod.pdf", 10,10)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor)
dev.off()

png("results/2025/bE2F_vs_bEMT.png", res=300, 1300, 1300)
ggplot(MEs_basal, aes(x=b_E2F_TARGETS, y=b_EPITHELIAL_MESENCHYMAL_TRANSITION))+geom_point()+geom_smooth(method="lm", se=F)+theme_classic()
dev.off()

pdf("results/2025/bE2F_vs_bEMT.pdf", 4,4)
ggplot(MEs_basal, aes(x=b_E2F_TARGETS, y=b_EPITHELIAL_MESENCHYMAL_TRANSITION))+geom_point()+geom_smooth(method="lm", se=F)+theme_classic()
dev.off()


