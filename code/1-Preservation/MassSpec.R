#Testing within module correlations in proteomic data
load(file="results/2025/centrality_global.RData")
load(file="results/2025/centrality_basal.RData")


MS2619<-read.csv("data/Mass Spec BC/proteinGroups_PXD002619.txt", sep="\t")

HLnorm<-MS2619[,grep("Ratio.H.L.normalized", colnames(MS2619))]

ks_p<-matrix(NA, nrow=length(unique(centrality_global$module)), ncol=3)
rownames(ks_p)<-unique(centrality_global$module)

set.seed(54987656)
for(module in unique(centrality_global$module)){

HLnorm_mod<-HLnorm[MS2619[,7] %in% rownames(centrality_global)[centrality_global$module==module],]
corr_mod<-cor(t(HLnorm_mod), use="pairwise.complete.obs")
corr_mod<-corr_mod[corr_mod!=1]

HLnorm_rand<-HLnorm[sample(1:nrow(HLnorm), nrow(HLnorm_mod)),]
corr_rand<-cor(t(HLnorm_rand), use="pairwise.complete.obs")
corr_rand<-corr_rand[corr_rand!=1]


ks_p[module, 1]<-ks.test(corr_mod, corr_rand, alternative="less")[[2]]
}

PXD000815<-read.csv("data/Mass Spec BC/proteinGroups_PXD000815.txt", sep="\t")

set.seed(59465496)
for(module in unique(centrality_global$module)){

HLnormPXD000815<-PXD000815[,grep("Ratio.H.L.normalized", colnames(PXD000815))]
HLnormPXD000815_mod<-HLnormPXD000815[as.character(PXD000815[,7]) %in% rownames(centrality_global)[centrality_global$module==module],]
corrPXD000815_mod<-cor(t(HLnormPXD000815_mod), use="pairwise.complete.obs")
HLnormPXD000815_rand<-HLnormPXD000815[sample(1:nrow(HLnormPXD000815), nrow(HLnormPXD000815_mod)),]
corrPXD000815_rand<-cor(t(HLnormPXD000815_rand), use="pairwise.complete.obs")
corrPXD000815_mod<-corrPXD000815_mod[corrPXD000815_mod!=1]
corrPXD000815_rand<-corrPXD000815_rand[corrPXD000815_rand!=1]


ks_p[module, 2]<-ks.test(corrPXD000815_mod, corrPXD000815_rand, alternative="less")[[2]]
}

####RPPA TCGA
RPPA<-read.csv("data/Mass Spec BC/RPPA TCGA/TCGA-BRCA-L4.csv")
which(colnames(RPPA) %in% rownames(centrality_global)[centrality_global$module=="E2F_targets"])
#RPPA_E2F<-RPPA[,grep(paste(rownames(centrality_global)[centrality_global$module=="E2F_targets"], collapse="|"), colnames(RPPA))]
RPPA_E2F<-RPPA[,which(colnames(RPPA) %in% rownames(centrality_global)[centrality_global$module=="E2F_targets"])]
RPPA_rand<-RPPA[,sample(5:ncol(RPPA),11)]

set.seed(2648696)
for(module in unique(centrality_global$module)){
RPPA_module<-RPPA[,which(colnames(RPPA) %in% rownames(centrality_global)[centrality_global$module==module])]

RPPA_rand<-RPPA[,sample(5:ncol(RPPA),length(which(colnames(RPPA) %in% rownames(centrality_global)[centrality_global$module==module])))]
if(length(which(colnames(RPPA) %in% rownames(centrality_global)[centrality_global$module==module]))>1){
cor_module<-cor(RPPA_module)
cor_module<-cor_module[cor_module!=1]
cor_rand<-cor(RPPA_rand)
cor_rand<-cor_rand[cor_rand!=1]
cor_module<-unlist(cor_module)
cor_rand<-unlist(cor_rand)


ks_p[module, 3]<-ks.test(cor_module, cor_rand, alternative="less")[[2]]
}
}


#text<-matrix(paste("p =", signif(ks_p,2)),nrow=nrow(ks_p), byrow=F)
#text[text=="p = NA"]<-""
#text[text=="p = 2.2e-16"]<-"p < 2.2e-16"

ks_p[which(ks_p==0)]<-2.2*10^(-16)
ks_p<-ks_p[-which(rownames(ks_p)=="Unconnected"),]

library(pheatmap)
##definire i colori per l'heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
myBreaks <- c(seq(min(unlist(-log10(ks_p)), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(unlist(-log10(ks_p)), na.rm=T)/paletteLength, max(unlist(-log10(ks_p)), na.rm=T), length.out=floor(paletteLength/2)))

png("Protein_preservation_global_modules.png", res=300, 1700, 2000)
pheatmap(-log10(ks_p), cluster_cols = F,  cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor)
dev.off()

