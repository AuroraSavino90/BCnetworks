library(limma)
library(openxlsx)

load(file="results/2025/centrality_basal.RData")
load(file="results/2025/surv_data.RData")
load(file="data/RData/metabric.RData")
load(file="data/RData/meta.RData")

TF<-read.csv("data/Transcription factor GO0003700.txt", header=T)


centrality_E2F<-centrality_basal[centrality_basal$module=="b_E2F_TARGETS",]
centrality_E2F_TF<-centrality_E2F[which(rownames(centrality_E2F) %in% TF[,3]),]

centrality_E2F_TF$surv_all_5y_dir<-surv_data_global_5years$direction[match(rownames(centrality_E2F_TF), surv_data_global_5years$gene)]
centrality_E2F_TF$surv_al_5y_p<-surv_data_global_5years$pval[match(rownames(centrality_E2F_TF), surv_data_global_5years$gene)]

##############################
### Compute genes differential expression between Basal-like and others
###############################

subtype<-meta$NOT_IN_OSLOVAL_Pam50Subtype
subtype<-ifelse(subtype=="Basal", "Basal", "nonBasal")
subtype<-factor(subtype, levels=c("nonBasal", "Basal"))

design <- model.matrix(~ 0+subtype)
colnames(design) <- c("nonBasal", "Basal")
cont_matrix <- makeContrasts(BasalvsNON = Basal-nonBasal, levels=design)

# Fit the expression matrix to a linear model
fit <- lmFit(metabric, design)
# Compute contrast
fit_contrast <- contrasts.fit(fit, cont_matrix)
# Bayes statistics of differential expression
# *There are several options to tweak!*
fit_contrast <- eBayes(fit_contrast)
# Generate a vocalno plot to visualize differential expression
volcanoplot(fit_contrast)
# Generate a list of top 100 differentially expressed genes
DEGs <- topTable(fit_contrast, number = nrow(fit_contrast), adjust = "BH")


centrality_E2F_TF$OE_log2FC<-DEGs$logFC[match(rownames(centrality_E2F_TF), rownames(DEGs))]
centrality_E2F_TF$OE_p<-DEGs$P.Value[match(rownames(centrality_E2F_TF), rownames(DEGs))]


write.xlsx(centrality_E2F_TF[order(centrality_E2F_TF$kWithin, decreasing=T),], "results/2025/TF_bE2F.xlsx", rowNames=T)

