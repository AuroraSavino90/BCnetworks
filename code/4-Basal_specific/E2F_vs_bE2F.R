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


i<-6 #E2F_targets
df<-data.frame(kwithin=c(Alldegrees_globalmod_basal[which(moduleColors_global==unique(moduleColors_global)[i]),2],
                         Alldegrees_basalmod_basal[which(moduleColors_basal==unique(moduleColors_basal)[i]),2]),
               module=c(rep("E2F_TARGETS", length(which(moduleColors_global==unique(moduleColors_global)[i]))),
                         rep("b_E2F_TARGETS", length(which(moduleColors_basal==unique(moduleColors_basal)[i])))
                       ))


p <- ggboxplot(df, x = "module", y = "kwithin")
#  Add p-value
p <- p + stat_compare_means(method = "wilcox.test") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

png("results/2025/E2F_globalvsbasal_kwithin.png", width = 1000, height = 2500, res = 300)
print(p)
dev.off()

pdf("results/2025/E2F_globalvsbasal_kwithin.pdf", 2,6)
print(p)
dev.off()

#E2F and b_E2F have different size, check that this does not influence the result
#909 global, 1262 basal e2f
#reduce b_E2F_targets to the same size of E2F_targets
moduleColors_basalred<-moduleColors_basal
moduleColors_basalred[which(centrality_basal$module=="b_E2F_TARGETS" & centrality_basal$rank>909)]<-"removed"
Alldegrees_basalmodred_basal=intramodularConnectivity(ADJ1, moduleColors_basalrand)

df<-data.frame(kwithin=c(Alldegrees_globalmod_basal[which(moduleColors_global==unique(moduleColors_global)[i]),2],
                         Alldegrees_basalmodred_basal[which(moduleColors_basalrand==unique(moduleColors_basal)[i]),2]),
               module=c(rep("E2F_TARGETS", length(which(moduleColors_global==unique(moduleColors_global)[i]))),
                        rep("b_E2F_TARGETS_red", length(which(moduleColors_basalred==unique(moduleColors_basal)[i])))))


p <- ggboxplot(df, x = "module", y = "kwithin")
#  Add p-value
p <- p + stat_compare_means(method = "wilcox.test") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print(p)

library(VennDiagram)
library(ggVennDiagram)
# Combine into a named list
gene_lists <- list(
  E2F_TARGETS = rownames(centrality_global)[which(moduleColors_global==unique(moduleColors_global)[i])],
  b_E2F_TARGETS = rownames(centrality_basal)[which(moduleColors_basal==unique(moduleColors_basal)[i])]
)

png("results/2025/E2F_globalvsbasal.png", width = 3000, height = 3500, res = 300)
ggVennDiagram(gene_lists,
              label_alpha = 0,
              edge_size = 0.5,
              category.names = names(gene_lists), label = "count") +
  scale_fill_gradient(low = "white", high = "#0073C2FF") +
  theme(text = element_text(size = 14),
        plot.margin = margin(2, 2, 2, 2, "cm")) + scale_x_continuous(expand = expansion(mult = .5)) # generous margin
dev.off()


pdf("results/2025/E2F_globalvsbasal.pdf", 7,7)
ggVennDiagram(gene_lists,
              label_alpha = 0,
              edge_size = 0.5,
              category.names = names(gene_lists), label = "count") +
  scale_fill_gradient(low = "white", high = "#0073C2FF") +
  theme(text = element_text(size = 14),
        plot.margin = margin(2, 2, 2, 2, "cm")) + scale_x_continuous(expand = expansion(mult = .5)) # generous margin
dev.off()
