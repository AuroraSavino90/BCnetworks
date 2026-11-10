library(WGCNA)
library(ggpubr)
load(file="data/RData/metabric.RData")
load(file="data/RData/meta.RData")
load(file="results/2025/centrality_basal.RData")
load(file="results/2025/centrality_global.RData")

###moduleColors_basal in other subtypes
moduleColors_global<-centrality_global$module
moduleColors_basal<-centrality_basal$module

ADJ1=adjacency(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]), power = 6, type= "signed")
Alldegrees_globalmod_basal=intramodularConnectivity(ADJ1, moduleColors_global)
Alldegrees_basalmod_basal=intramodularConnectivity(ADJ1, moduleColors_basal)
# kWithin divided by (module size - 1), i.e. mean adjacency: comparable between modules of different size
Alldegrees_globalmod_basal$kWithin<-Alldegrees_globalmod_basal$kWithin/(table(moduleColors_global)[moduleColors_global]-1)
Alldegrees_basalmod_basal$kWithin<-Alldegrees_basalmod_basal$kWithin/(table(moduleColors_basal)[moduleColors_basal]-1)


i<-6 #E2F_targets
stopifnot(unique(moduleColors_global)[i]=="E2F_TARGETS", unique(moduleColors_basal)[i]=="b_E2F_TARGETS")
df<-data.frame(kwithin=c(Alldegrees_globalmod_basal[which(moduleColors_global==unique(moduleColors_global)[i]),2],
                         Alldegrees_basalmod_basal[which(moduleColors_basal==unique(moduleColors_basal)[i]),2]),
               module=c(rep("E2F_TARGETS", length(which(moduleColors_global==unique(moduleColors_global)[i]))),
                         rep("b_E2F_TARGETS", length(which(moduleColors_basal==unique(moduleColors_basal)[i])))
                       ))


p <- ggboxplot(df, x = "module", y = "kwithin", ylab = "kWithin / (module size - 1)")
#  Add p-value
p <- p + stat_compare_means(method = "wilcox.test") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

png("results/2025_fix/E2F_globalvsbasal_kwithin.png", width = 1000, height = 2500, res = 300)
print(p)
dev.off()

pdf("results/2025_fix/E2F_globalvsbasal_kwithin.pdf", 2,6)
print(p)
dev.off()

#E2F and b_E2F have different size, check that this does not influence the result
#909 global, 1262 basal e2f
#random subsets of b_E2F_targets of the same size of E2F_targets, raw kWithin recomputed within each subset
g_b<-which(moduleColors_basal=="b_E2F_TARGETS")
n_g<-length(which(moduleColors_global=="E2F_TARGETS"))
k_glob<-mean(rowSums(ADJ1[which(moduleColors_global=="E2F_TARGETS"), which(moduleColors_global=="E2F_TARGETS")])-1)
set.seed(1)
k_red<-replicate(1000, {s<-sample(g_b, n_g); mean(rowSums(ADJ1[s, s])-1)})
write.csv(data.frame(mean_kwithin_E2F=k_glob, mean_kwithin_bE2F_sizematched=mean(k_red), q025=quantile(k_red, 0.025), q975=quantile(k_red, 0.975), p_emp=(1+sum(k_red<=k_glob))/1001),
          "results/2025_fix/E2F_globalvsbasal_kwithin_sizematched.csv", row.names=FALSE)

### same comparison in independent basal-like samples (TCGA)
load("../../../Data/Data_downloaded/TCGA/TCGA-BRCA_Primary Tumor_dedupl.RData")
library(TCGAbiolinks)
BRCA_subtype <- TCGAquery_subtype(tumor = "BRCA")
basal_tcga<-which(BRCA_subtype$BRCA_Subtype_PAM50[match(substr(colnames(TCGA_dedupl),1,12), BRCA_subtype$patient)]=="Basal")
mods_tcga<-list(E2F_TARGETS=intersect(rownames(centrality_global)[which(moduleColors_global=="E2F_TARGETS")], rownames(TCGA_dedupl)),
                b_E2F_TARGETS=intersect(rownames(centrality_basal)[which(moduleColors_basal=="b_E2F_TARGETS")], rownames(TCGA_dedupl)))
mods_tcga<-lapply(mods_tcga, function(g) g[apply(TCGA_dedupl[g, basal_tcga], 1, sd)>0]) # drop genes constant in TCGA basal (NA correlations)
df<-do.call(rbind, lapply(names(mods_tcga), function(m) data.frame(module=m,
       kwithin=(rowSums(adjacency(t(TCGA_dedupl[mods_tcga[[m]], basal_tcga]), power = 6, type= "signed"))-1)/(length(mods_tcga[[m]])-1))))

p <- ggboxplot(df, x = "module", y = "kwithin", ylab = "kWithin / (module size - 1)")
p <- p + stat_compare_means(method = "wilcox.test") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

png("results/2025_fix/E2F_globalvsbasal_kwithin_TCGA.png", width = 1000, height = 2500, res = 300)
print(p)
dev.off()

pdf("results/2025_fix/E2F_globalvsbasal_kwithin_TCGA.pdf", 2,6)
print(p)
dev.off()

library(ggVennDiagram)
# Combine into a named list
gene_lists <- list(
  E2F_TARGETS = rownames(centrality_global)[which(moduleColors_global==unique(moduleColors_global)[i])],
  b_E2F_TARGETS = rownames(centrality_basal)[which(moduleColors_basal==unique(moduleColors_basal)[i])]
)

png("results/2025_fix/E2F_globalvsbasal.png", width = 3000, height = 3500, res = 300)
print(ggVennDiagram(gene_lists,
              label_alpha = 0,
              edge_size = 0.5,
              category.names = names(gene_lists), label = "count") +
  scale_fill_gradient(low = "white", high = "#0073C2FF") +
  theme(text = element_text(size = 14),
        plot.margin = margin(2, 2, 2, 2, "cm")) + scale_x_continuous(expand = expansion(mult = .5))) # generous margin
dev.off()


pdf("results/2025_fix/E2F_globalvsbasal.pdf", 7,7)
print(ggVennDiagram(gene_lists,
              label_alpha = 0,
              edge_size = 0.5,
              category.names = names(gene_lists), label = "count") +
  scale_fill_gradient(low = "white", high = "#0073C2FF") +
  theme(text = element_text(size = 14),
        plot.margin = margin(2, 2, 2, 2, "cm")) + scale_x_continuous(expand = expansion(mult = .5))) # generous margin
dev.off()
