library(WGCNA)
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
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
myBreaks <- seq(-1, 1, length.out=paletteLength + 1) # symmetric scale

png("results/2025_fix/Heatmap_b_cormod.png", res=300, 3000, 3000)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, clustering_distance_rows=as.dist(1-toplot), clustering_distance_cols=as.dist(1-toplot))
dev.off()

pdf("results/2025_fix/Heatmap_b_cormod.pdf", 10,10)
pheatmap(toplot, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, clustering_distance_rows=as.dist(1-toplot), clustering_distance_cols=as.dist(1-toplot))
dev.off()

ct<-cor.test(MEs_basal$b_E2F_TARGETS, MEs_basal$b_EPITHELIAL_MESENCHYMAL_TRANSITION)
write.csv(data.frame(r=ct$estimate, p=ct$p.value, n=nrow(MEs_basal)), "results/2025_fix/bE2F_vs_bEMT_cor.csv", row.names=FALSE)

png("results/2025_fix/bE2F_vs_bEMT.png", res=300, 1300, 1300)
print(ggplot(MEs_basal, aes(x=b_E2F_TARGETS, y=b_EPITHELIAL_MESENCHYMAL_TRANSITION))+geom_point()+geom_smooth(method="lm", se=F)+theme_classic())
dev.off()

pdf("results/2025_fix/bE2F_vs_bEMT.pdf", 4,4)
print(ggplot(MEs_basal, aes(x=b_E2F_TARGETS, y=b_EPITHELIAL_MESENCHYMAL_TRANSITION))+geom_point()+geom_smooth(method="lm", se=F)+theme_classic())
dev.off()


