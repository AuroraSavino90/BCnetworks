set.seed(123)
####setting paths
pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20250221/'

###loading input data
gene_annot <- read.csv("data/CRISPR/gene_identifiers_20241212.csv")
### gene_identifiers_20191101 downloaded from https://cog.sanger.ac.uk/cmp/download/gene_identifiers_20241212.csv on 20250221

CMP_annot <- read.csv("data/CRISPR/model_list_20241120.csv")
### model_list_20240110.csv downloaded from https://cog.sanger.ac.uk/cmp/download/model_list_20241120.csv on 20250129

#select only breast carcinoma cell lines
BC_lines<-CMP_annot$model_name[CMP_annot$cancer_type=="Breast Carcinoma"]


scaled_depFC<-read.csv('data/CRISPR//CRISPRGeneEffect.csv', row.names = 1)
colnames(scaled_depFC)<-gsub("\\..*","",colnames(scaled_depFC))
###scaled essentiality matrices downloaded from https://depmap.org/portal/data_page/?tab=allData on 20250129 (24Q4)

toremove<-which(is.na(CMP_annot$model_id[match(rownames(scaled_depFC),CMP_annot$BROAD_ID)]))
scaled_depFC<-scaled_depFC[-toremove,]
rownames(scaled_depFC)<-CMP_annot$model_id[match(rownames(scaled_depFC),CMP_annot$BROAD_ID)]

scaled_depFC<-t(scaled_depFC)

CMP_annot$model_name[CMP_annot$cancer_type=="Breast Carcinoma"][which(CMP_annot$model_id[CMP_annot$cancer_type=="Breast Carcinoma"] %in% colnames(scaled_depFC))]

BL<-c("HCC1143", "Hs-578-T", "BT-549", "MDA-MB-231", "MFM-223",
      "MDA-MB-157", "HCC38", "HCC70", "MDA-MB-468", "HCC1937", "HCC1187", "HCC1806",
      "CAL-120", "COLO-824", "DU-4475", "SUM-185PE", "SUM-159PT", "SUM-149PT",
      "SUM-229PE", "HMC-1-8", "VP229", "SUM-1315MO2", "MDA-MB-436", "HCC1395", 
      "SUM-102PT", "CAL-51")
BL_id<-CMP_annot$model_id[CMP_annot$model_name %in% BL]

###get CRISPR data for BL lines
scaled_depFC<-scaled_depFC[,BL_id]

#remove genes with missing values
scaled_depFC<-scaled_depFC[-which(rowSums(is.na(scaled_depFC))>0),]



###corelation between essentiality and centrality for each module
load(file="results/2025/centrality_basal.RData")

##keep only genes in common between centrality_basal and essentiality matrices
centrality_basal<-centrality_basal[which(rownames(centrality_basal) %in% rownames(scaled_depFC)),]

pval<-matrix(nrow=length(unique(centrality_basal$module)), ncol=ncol(scaled_depFC))
cc<-matrix(nrow=length(unique(centrality_basal$module)), ncol=ncol(scaled_depFC))
for(j in 1:length(unique(centrality_basal$module))){
  for(i in 1:ncol(scaled_depFC)){
    pval[j,i]<-cor.test(centrality_basal$kWithin[centrality_basal$module==unique(centrality_basal$module)[j]], scaled_depFC[rownames(centrality_basal)[centrality_basal$module==unique(centrality_basal$module)[j]],i])[[3]]
    cc[j,i]<-cor.test(centrality_basal$kWithin[centrality_basal$module==unique(centrality_basal$module)[j]], scaled_depFC[rownames(centrality_basal)[centrality_basal$module==unique(centrality_basal$module)[j]],i])[[4]]
  }
}
rownames(cc)<-unique(centrality_basal$module)
colnames(cc)<-colnames(scaled_depFC)
rownames(pval)<-unique(centrality_basal$module)
colnames(pval)<-colnames(scaled_depFC)

#remove unconnected
pval<-pval[-1,]
cc<-cc[-1,]

#compute average correlation
cc_mean<-rowMeans(cc)

#adjust pvalue
adj_pval<-matrix(p.adjust(pval, method = "BH"), nrow=nrow(pval))
colnames(adj_pval)<-colnames(pval)
rownames(adj_pval)<-rownames(pval)
#cut pvalue at 2.2*10^-16
adj_pval[adj_pval<2.2*10^(-16)]<-2.2*10^(-16)

#merge pvalues with the Fisher method
library(metap)
pmerged<-c()
for(r in 1:nrow(adj_pval)){
  istwo <- rep(T, ncol(adj_pval))
  toinvert <- ifelse(cc[r,]>0, T, F)
pmerged<-c(pmerged, sumlog(two2one(adj_pval[r,], two = istwo, invert = toinvert))$p)
}
pmerged[pmerged==0]<-10^(-299)

names(pmerged)<-rownames(adj_pval)

cc[adj_pval>0.05]<-NA
cc<-t(cc)
rownames(cc)<-CMP_annot$model_name[match(rownames(cc), CMP_annot$model_id)]


library(pheatmap)

paletteLength <- 50
myColor <- colorRampPalette(c("#4575B4", "white"))(paletteLength)
myBreaks <- c(seq(min(unlist(cc), na.rm=T),0, length.out=floor(paletteLength)))

anno_p<-data.frame(meanR= -cc_mean, mergedp= -log10(pmerged))

# Define a continuous color gradient
continuous_colors <- colorRampPalette(c("white", "mediumorchid4"))(100)

# Map the continuous annotation to colors
annotation_colors <- list(meanR = continuous_colors)


graphics.off()
png("results/2025/CRISPR_corr.png", res=300, 3000, 3000)
pheatmap((cc[,names(sort(colSums(cc, na.rm=T)))]), cluster_cols = F, cluster_rows = F, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, annotation_col = anno_p, annotation_colors = annotation_colors)
dev.off()


