#############load network objects and data
load("data/RData/metabric.RData")
load("data/RData/meta.RData")
load("results/2025/centrality_basal.RData")
load("results/2025/centrality_global.RData")

#############compute module eigengenes

MEs= moduleEigengenes(t(metabric), centrality_global$module)$eigengenes
colnames(MEs)<-gsub("^ME","", colnames(MEs))
MEs_basal= moduleEigengenes(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]), centrality_basal$module)$eigengenes
colnames(MEs_basal)<-gsub("^ME","", colnames(MEs_basal))

png("results/2025/modulegrade_e2f.png", res = 300, width=1000, height = 1000)
boxplot(MEs[,"E2F_TARGETS"]~meta$grade, outline=F, xlab="grade", ylab="Module Eigengene")
dev.off()

png("results/2025/modulegrade_estrogen.png", res = 300, width=1000, height = 1000)
boxplot(MEs[,"ESTROGEN_RESPONSE_EARLY"]~meta$grade, outline=F, xlab="grade", ylab="Module Eigengene")
dev.off()