library(DESeq2)
library(FactoMineR)
library(ggplot2)
library(openxlsx)
library(pheatmap)
library(WGCNA)
library(ggrepel)
library(ggpubr)

date<-"2025"

######network data to load
load("results/2025/centrality_basal.RData")
load("data/RData/metabric.RData")
load("data/RData/meta.RData")

########################################
##### FUNCTIONS #########################
#########################################


#####################################################################
##function to change gene names (e.g. from ENSEMBL to gene symbol)
###################################################################

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
      }
    }
    
    data=data[which(is.na(annotation_sel)==F),]
    annotation_sel=na.omit(annotation_sel)
    a<-which(duplicated(annotation_sel))
  }
  
  rownames(data)=annotation_sel
  return(data)
}


#####################################################################
##function to select DEGs with padj<x, up or downregulated
###################################################################

DEGsfilt<-function(DEGs, padj=0.05, FC="down"){
  DE_filt<-DEGs[which(DEGs$padj<=padj),]
  if(FC=="down"){
    DE_filt<-rownames(DE_filt)[DE_filt$log2FoldChange<0]
  } else if(FC=="up"){
    DE_filt<-rownames(DE_filt)[DE_filt$log2FoldChange>0]
  }
  return(DE_filt)
}

#####################################################################
##function to test the enrichment in basal modules
###################################################################

fishertest_basal<-function(genes, dataset){
  fisherp<-matrix(ncol=length(unique(centrality_basal$module)), nrow=1)
  for(i in 1:length(unique(centrality_basal$module))){
    metabric_tmp<-metabric[rownames(metabric) %in% rownames(dataset),]
    moduleColors_basal_tmp<-centrality_basal$module[rownames(metabric) %in% rownames(dataset)]
    counts<-matrix(c(length(intersect(rownames(metabric_tmp)[which(moduleColors_basal_tmp %in% unique(centrality_basal$module)[i])], genes )),
                     
                     
                     length(which(moduleColors_basal_tmp %in% unique(centrality_basal$module)[i]))-length(intersect(rownames(metabric_tmp)[which(moduleColors_basal_tmp %in% unique(centrality_basal$module)[i])], genes )),
                     length(intersect(genes, rownames(metabric_tmp)[-which(moduleColors_basal_tmp %in% c("grey", unique(centrality_basal$module)[i]))])),
                     length(setdiff(rownames(metabric_tmp), union(rownames(metabric_tmp)[which(moduleColors_basal_tmp %in% c("grey", unique(centrality_basal$module)[i]))], genes)))), nrow=2)
    fisherp[1,i]<-fisher.test(counts, alternative = "greater")[[1]]
  }
  
  colnames(fisherp)<-unique(centrality_basal$module)
  return(fisherp)
}

#####################################################################
##function to test the enrichment in basal modules for a list of datasets
###################################################################

fishertest_alldat<-function(alldat=c("dds468",
                                     "ddsHs",
                                     "dds231TFDP1", "dds231E2F3"),
                            names_alldat=c("MDAMB468 TFDP1", "Hs578 TFDP1", "MDAMB231 TFDP1", "MDAMB231 E2F3")){
  ft_up_tot<-matrix(nrow=20, ncol=length(alldat))
  ft_down_tot<-matrix(nrow=20, ncol=length(alldat))
  ft_all_tot<-matrix(nrow=20, ncol=length(alldat))
  pos<-0
  for(dat in alldat){
    i<-get(dat)
    pos<-pos+1
    i_down<-DEGsfilt(DEGs=i, padj=0.05, FC="down")
    i_up<-DEGsfilt(DEGs=i, padj=0.05, FC="up")
    ft_up<-fishertest_basal(i_up, counts_name)
    ft_down<-fishertest_basal(i_down, counts_name)
    ft_all<-fishertest_basal(c(i_up,i_down), counts_name)
    ft_up_tot[,pos]<-ft_up
    ft_down_tot[,pos]<-ft_down
    ft_all_tot[,pos]<-ft_all
  }
  
  rownames(ft_up_tot)<-colnames(ft_up)
  rownames(ft_down_tot)<-colnames(ft_up)
  rownames(ft_all_tot)<-colnames(ft_up)
  colnames(ft_up_tot)<-names_alldat
  colnames(ft_down_tot)<-names_alldat
  colnames(ft_all_tot)<-names_alldat
  
  
  ft_up_tot<-ft_up_tot[-1,]
  ft_down_tot<-ft_down_tot[-1,]
  ft_all_tot<-ft_all_tot[-1,]
  
  for(row in 1:nrow(ft_up_tot)){
    ft_up_tot[row,]<-p.adjust(ft_up_tot[row,])
  }
  for(row in 1:nrow(ft_down_tot)){
    ft_down_tot[row,]<-p.adjust(ft_down_tot[row,])
  }
  for(row in 1:nrow(ft_all_tot)){
    ft_all_tot[row,]<-p.adjust(ft_all_tot[row,])
  }
  
  ft_up_tot[ft_up_tot<2.2*10^(-16)]<-2.2*10^(-16)
  ft_down_tot[ft_down_tot<2.2*10^(-16)]<-2.2*10^(-16)
  ft_all_tot[ft_all_tot<2.2*10^(-16)]<-2.2*10^(-16)
  
  return(list(ft_up_tot, ft_down_tot, ft_all_tot))
}

#####################################################################
###function to project MEs computed on a dataset (original_data) on new transcriptomic data (newdata)
#####################################################################

pcaproject<-function(newdata, original_data, modules, ME=eigengenes){
  pcaproj<-matrix(nrow=ncol(newdata), ncol=length(unique(modules)))
  for(i in 1:length(unique(modules))){
    pca <- prcomp(t(original_data[modules==unique(modules)[i],]))
    if(cor(pca$x[,1], ME[,unique(modules)[i]])<0){
      pca$rotation[,1]<-(-pca$rotation[,1])
    }
    commongenes<-rownames(newdata)[which(rownames(newdata) %in% rownames(original_data[modules==unique(modules)[i],]))]
    newdata_forpca<-newdata[commongenes,]
    pcaproj[,i]<-colSums(t(scale(t(newdata_forpca), pca$center[commongenes], pca$scale)) * c(pca$rotation[commongenes,1]), na.rm=T)
  }
  return(pcaproj)
}



########################################
##### ANALYSES #########################
#########################################


#######################################################
### normalization, filtering, log transformation
###################################################

counts<-read.xlsx("data/RNAseq/gene_count.xlsx")
rownames(counts)<-counts[,1]
anno<-counts[,c("gene_id", "gene_name")]
counts_name<-changenames(counts[,c(2:34)], anno = anno)

counts2<-read.table("data/RNAseq/Poli_E2F3_KO-RNAseq-v1-run241011/RNAseq/dataset/v1-run241011/GEP.count", header = T, row.names = 1)

genes<-intersect(rownames(counts_name), rownames(counts2))

counts_tot<-cbind(counts_name[genes,], counts2[genes,])

RPM<-t(t(counts_tot)/colSums(counts_tot))*1000000

counts_tot<-counts_tot[rowSums(counts_tot>=10)>2,]
RPM<-RPM[rownames(counts_tot),]
RPMlog<-log2(RPM+1)

###########################
### quality checks: PCA
##########################

metadata<-read.xlsx("data/RNAseq/Novogene_metadata.xlsx", rowNames = T)
metadata$Clone<-factor(metadata$Clone)
metadata$Seq<-"Novogene"

metadata2<-data.frame(Cell.line=rep("MDAMB468", 9), KO.gene=rep(c("EV", "E2F3", "E2F3"), each=3),
                     Clone=rep(c(100, 1, 3), each=3), Replicate=rep(1:3, 3), Seq=rep("Oliviero", 9))

metadata<-rbind.data.frame(metadata, metadata2)
metadata$Clone<-factor(metadata$Clone, levels=c(levels(metadata$Clone), 1, 3))

metadata$Clone[37:39]<-1
metadata$Clone[40:42]<-3

rownames(metadata)<-colnames(RPMlog)

pca<-PCA(t(RPMlog))
df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2], PC3=pca$ind$coord[,3],
               PC4=pca$ind$coord[,4], PC5=pca$ind$coord[,5],
               metadata)

png(paste("results/",date, "/PCA.png", sep=""), res=300, 1300, 1000)
ggplot(df, aes(x=PC1, y=PC2, colour=Cell.line))+geom_point()+theme_classic()
dev.off()


################## PCA separating cell lines

pca231<-PCA(t(RPMlog[,metadata$Cell.line=="MDAMB231"]))
df<-data.frame(PC1=pca231$ind$coord[,1], PC2=pca231$ind$coord[,2], PC3=pca231$ind$coord[,3],
               PC4=pca231$ind$coord[,4], PC5=pca231$ind$coord[,5],
               metadata[metadata$Cell.line=="MDAMB231",])

png(paste("results/",date, "/PCA231.png", sep=""), res=300, 1300, 1000)
ggplot(df, aes(x=PC1, y=PC2, colour=KO.gene, shape=as.factor(Clone)))+geom_point()+theme_classic()
dev.off()

pca468<-PCA(t(RPMlog[,metadata$Cell.line=="MDAMB468"]))
df<-data.frame(PC1=pca468$ind$coord[,1], PC2=pca468$ind$coord[,2], PC3=pca468$ind$coord[,3],
               PC4=pca468$ind$coord[,4], PC5=pca468$ind$coord[,5],
               metadata[metadata$Cell.line=="MDAMB468",])

png(paste("results/",date, "/PCA468.png", sep=""), res=300, 1300, 1000)
ggplot(df, aes(x=PC1, y=PC2, colour=KO.gene, shape=as.factor(Clone)))+geom_point()+theme_classic()
dev.off()

pcaHs<-PCA(t(RPMlog[,metadata$Cell.line=="Hs578"]))
df<-data.frame(PC1=pcaHs$ind$coord[,1], PC2=pcaHs$ind$coord[,2], PC3=pcaHs$ind$coord[,3],
               PC4=pcaHs$ind$coord[,4], PC5=pcaHs$ind$coord[,5],
               metadata[metadata$Cell.line=="Hs578",])

png(paste("results/",date, "/PCAHs.png", sep=""), res=300, 1300, 1000)
ggplot(df, aes(x=PC1, y=PC2, colour=KO.gene, shape=as.factor(Clone)))+geom_point()+theme_classic()
dev.off()


#######################################################
########### DEGs
#########################################################

#clone can't be used as a covariate as we have only one clone for controls, leading to nested variables

############## pooling clones

dds468 <- DESeqDataSetFromMatrix(countData = counts_tot[,metadata$Cell.line=="MDAMB468"],
                                 colData = metadata[metadata$Cell.line=="MDAMB468",],
                                 design= ~ KO.gene)
dds468 <- DESeq(dds468)
dds468TFDP1 <- results(dds468, contrast = c("KO.gene","TFDP1", "EV"))
dds468E2F3 <- results(dds468, contrast = c("KO.gene","E2F3", "EV"))

ddsHs <- DESeqDataSetFromMatrix(countData = counts_tot[,metadata$Cell.line=="Hs578"],
                                   colData = metadata[metadata$Cell.line=="Hs578",],
                                   design= ~ KO.gene)
ddsHs <- DESeq(ddsHs)
ddsHs <- results(ddsHs, contrast = c("KO.gene","TFDP1", "EV"))

dds231 <- DESeqDataSetFromMatrix(countData = counts_tot[,metadata$Cell.line=="MDAMB231"],
                                    colData = metadata[metadata$Cell.line=="MDAMB231",],
                                    design= ~ KO.gene)
dds231 <- DESeq(dds231)
dds231TFDP1 <- results(dds231, contrast = c("KO.gene","TFDP1", "EV"))
dds231E2F3 <- results(dds231, contrast = c("KO.gene","E2F3", "EV"))


#################################
### DEGs overlap
##################################
DEGs_shared<-cbind(ddsHs$log2FoldChange, dds468TFDP1$log2FoldChange,dds468E2F3$log2FoldChange, dds231TFDP1$log2FoldChange, dds231E2F3$log2FoldChange)
DEGs_shared[which(ddsHs$padj>0.05),1]<-0
DEGs_shared[which(dds468TFDP1$padj>0.05),2]<-0
DEGs_shared[which(dds468E2F3$padj>0.05),3]<-0
DEGs_shared[which(dds231TFDP1$padj>0.05),4]<-0
DEGs_shared[which(dds231E2F3$padj>0.05),5]<-0
rownames(DEGs_shared)<-rownames(ddsHs)
colnames(DEGs_shared)<-c("Hs578 TFDP1", "MDAMB468 TFDP1","MDAMB468 E2F3", "MDAMB231 TFDP1", "MDAMB231 E2F3")

toplot<-DEGs_shared[DEGs_shared[,1]!=0,]
toplot<-toplot[-which(rowSums(is.na(toplot))==5),]
toplot<-toplot[which(rownames(toplot) %in% rownames(centrality_basal)[which(centrality_basal$module=="b_E2F_TARGETS")]),]
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

graphics.off()
png(paste("results/",date, "/shared_DEGs_Novogene.png", sep=""),res=300, 1500, 50000)
pheatmap(toplot,   cellwidth=15, cellheight=5,  keep.dendro=T, color =myColor, breaks = myBreaks)
dev.off()

##########################
### Enrichment test
###########################


############## pooling clones

ft<-fishertest_alldat(alldat=c("ddsHs", "dds231TFDP1", "dds231E2F3","dds468TFDP1","dds468E2F3"),
                  names_alldat=c("Hs578 TFDP1","MDAMB231 TFDP1", "MDAMB231 E2F3", "MDAMB468 TFDP1", "MDAMB468 E2F3"))

library(metap)
istwo <- rep(F, nrow(ft[[1]]))
toinvert <- rep(F, nrow(ft[[1]]))
pmerged<-c()
for(row in 1:nrow(ft[[1]])){
  pmerged<-c(pmerged, sumlog(two2one(ft[[1]][row,], two = istwo, invert = toinvert))$p)
}
names(pmerged)<-rownames(ft[[1]])
anno_p<-data.frame(mergedp= -log10(pmerged), npositive=rowSums(ft[[1]]<0.05))
rownames(anno_p)<-rownames(ft[[1]])

graphics.off()
png(paste("results/",date, "/Enrich_up_Novogene.png", sep=""), res=300, 2500, 2500)
pheatmap(-log10(ft[[1]]), cellwidth=15, cellheight=15, keep.dendro=T, annotation_row = anno_p)
dev.off()

istwo <- rep(F, nrow(ft[[2]]))
toinvert <- rep(F, nrow(ft[[2]]))
pmerged<-c()
for(row in 1:nrow(ft[[2]])){
  pmerged<-c(pmerged, sumlog(two2one(ft[[2]][row,], two = istwo, invert = toinvert))$p)
}
names(pmerged)<-rownames(ft[[2]])
anno_p<-data.frame(mergedp= -log10(pmerged), npositive=rowSums(ft[[2]]<0.05))
rownames(anno_p)<-rownames(ft[[2]])

png(paste("results/",date, "/Enrich_down_Novogene.png", sep=""), res=300, 2500, 2500)
pheatmap(-log10(ft[[2]]), cellwidth=15, cellheight=15, keep.dendro=T, annotation_row = anno_p)
dev.off()

istwo <- rep(F, nrow(ft[[3]]))
toinvert <- rep(F, nrow(ft[[3]]))
pmerged<-c()
for(row in 1:nrow(ft[[3]])){
  pmerged<-c(pmerged, sumlog(two2one(ft[[3]][row,], two = istwo, invert = toinvert))$p)
}
names(pmerged)<-rownames(ft[[3]])
anno_p<-data.frame(mergedp= -log10(pmerged), npositive=rowSums(ft[[3]]<0.05))
rownames(anno_p)<-rownames(ft[[3]])

png(paste("results/",date, "/Enrich_all_Novogene.png", sep=""), res=300, 2500, 2500)
pheatmap(-log10(ft[[3]]), cellwidth=15, cellheight=15, keep.dendro=T, annotation_row = anno_p)
dev.off()

png(paste("results/",date, "/Enrich_all_Novogene_scaled.png", sep=""), res=300, 2500, 2500)
pheatmap(-log10(ft[[3]]),cellwidth=15, cellheight=15, scale="column")
dev.off()

#########################################################
### projection of b_E2F_targets MEs on new data  
#########################################################
MEs<-moduleEigengenes(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]), centrality_basal$module)$eigengenes
colnames(MEs)<-gsub("ME", "", colnames(MEs))

#####project MEs on Novogene data

pcaproj_hubs<-pcaproject(newdata=RPMlog, original_data=metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], modules=colnames(MEs), ME=MEs)
colnames(pcaproj_hubs)<-colnames(MEs)
rownames(pcaproj_hubs)<-colnames(RPMlog)

pcaproj_hubs<-cbind.data.frame(pcaproj_hubs, metadata)
pcaproj_hubs$KO.gene<-factor(pcaproj_hubs$KO.gene, levels=c("EV", "TFDP1", "E2F3"))

library(rstatix)
stat.test <- pcaproj_hubs %>%
  group_by(Cell.line) %>%
  wilcox_test(b_E2F_TARGETS ~ KO.gene, ref.group = "EV")

stat.test <- stat.test %>% add_y_position()


p <- ggboxplot(pcaproj_hubs, x = "KO.gene", y = "b_E2F_TARGETS",
               color = "KO.gene", palette = "jco",
               facet.by = "Cell.line", 
               add = "jitter")
#  Add p-value
p<-p +
  stat_pvalue_manual(stat.test, label = "p = {p.adj}", tip.length = 0.01)


####plot b_E2F_targets across KOs
png(paste("results/",date, "/bE2F_project_Novogene.png", sep=""), res=300, 2000, 2000)
print(p)
dev.off()


stat.test <- pcaproj_hubs %>%
  group_by(Cell.line) %>%
  wilcox_test(b_E2F_TARGETS ~ KO.gene, ref.group = "EV")

stat.test <- stat.test %>% add_y_position(scales="free")

p <- ggboxplot(pcaproj_hubs, x = "KO.gene", y = "b_E2F_TARGETS",
               color = "KO.gene", palette = "jco",
               add = "jitter")+ scale_y_continuous(expand = c(.1, .1)) 

p<-facet(p, facet.by = "Cell.line", ncol = 1, scales="free")

p<-p+
  stat_pvalue_manual(stat.test, label = "p = {p.adj}", tip.length = 0.01) + 
  theme(strip.text.x = element_text(size = 12))

####plot b_E2F_targets across KOs
png(paste("results/",date, "/bE2F_project_Novogene_long.png", sep=""), res=300, 1000, 3000)
print(p)
dev.off()


#######compute the cohens'd

cd<-matrix(nrow=3, ncol=19)
ind<-0
for(l in unique(pcaproj_hubs$Cell.line)){
  ind<-ind+1
  
  ctrl<-subset(pcaproj_hubs, Cell.line==l &  KO.gene=="EV")
  trt<-subset(pcaproj_hubs, Cell.line==l &  KO.gene=="TFDP1")
  
  for(i in 1:19){
    cd[ind,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
  }
  
}

colnames(cd)<-colnames(pcaproj_hubs)[1:19]
rownames(cd)<-paste(unique(pcaproj_hubs$Cell.line), "TFDP1")


cd2<-matrix(nrow=2, ncol=19)
ind<-1

ctrl<-subset(pcaproj_hubs, Cell.line=="MDAMB231" &  KO.gene=="EV")
trt<-subset(pcaproj_hubs, Cell.line=="MDAMB231" &  KO.gene=="E2F3")

for(i in 1:19){
  cd2[ind,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}

ind<-2

ctrl<-subset(pcaproj_hubs, Cell.line=="MDAMB468" &  KO.gene=="EV")
trt<-subset(pcaproj_hubs, Cell.line=="MDAMB468" &  KO.gene=="E2F3")

for(i in 1:19){
  cd2[ind,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}


colnames(cd2)<-colnames(pcaproj_hubs)[1:19]
rownames(cd2)<-c("MDAMB231 E2F3", "MDAMB468 E2F3")

cd_all<-rbind(cd, cd2)

########## plot changes in MEs (Cohen's d) for each KO
## Hs578 behave oppositely to other cell lines (which behave as expected)

toplot<-cd_all
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

graphics.off()
png(paste("results/",date, "/Cohen_bE2F_Novogene.png", sep=""),res=300, 1500, 2500)
pheatmap(t(toplot),   cellwidth=15, cellheight=15,  keep.dendro=T, color =myColor, breaks = myBreaks)
dev.off()



######################
## correlation with functional impact
#####################
load(file=paste("results/",date, "/d.RData", sep=""))

df<-data.frame(cd=abs(toplot[c(2,3,4,1,5),"b_E2F_TARGETS"]), fi=scale(d[3,])+scale(d[1,]), condition=rownames(toplot)[c(2,3,4,1,5)])
png(paste("results/",date, "/Cohen_bE2F_vs_functional.png", sep=""),res=300, 1000, 1000)
ggplot(df, aes(x=cd, y=fi, label=condition))+geom_point(size=2)+geom_smooth(method="lm", se=F)+geom_label_repel()+theme_classic()+
  xlab("Cohen's d (absolute value)") + ylab("Scaled functional impact")
dev.off()
cor.test(abs(toplot[c(2,3,4,1,5),"b_E2F_TARGETS"]), scale(d[3,])+scale(d[1,]))

########## plot changes in MEs (Cohen's d) vs modules' correlation with b_E2F_targets
## for each KO, the most affected modules are either the most highly or lowly correlated with b_E2F_targets

df<-data.frame(cd=c(t(cd_all)), condition=rep(c("MDAMB468 TFDP1", "Hs578 TFDP1", "MDAMB231 TFDP1", "MDAMB231 E2F3", "MDAMB468 E2F3"), each= 19), corr=rep(cc["b_E2F_TARGETS",colnames(cd_all)],5),
               module=rep(colnames(cd_all),5))
df$condition<-factor(df$condition, levels=c("Hs578 TFDP1", "MDAMB231 TFDP1", "MDAMB231 E2F3","MDAMB468 TFDP1", "MDAMB468 E2F3"))

png(paste("results/",date, "/CohenVScorr_bE2F_Novogene.png", sep=""),res=300, 4500, 2000)
ggplot(df, aes(x=corr, y=cd, label=module))+geom_point(size=2)+facet_grid(~condition)+geom_smooth(method = lm)+stat_cor(label.x=-0.5, label.y = 17)+geom_text_repel(max.overlaps = 5)+theme_bw()+theme(strip.text=element_text(size = 12, face = "bold"))
dev.off()



######################################
######## GO enrichment
######################################

library(clusterProfiler)
library(org.Hs.eg.db)

ego_up<-list()
ego_dn<-list()
for(c in c(
             "ddsHs",
             "dds231TFDP1", "dds231E2F3",
             "dds468TFDP1","dds468E2F3")){
  
  i<-get(c)
 i_down<-DEGsfilt(DEGs=i, padj=0.05, FC="down")
  i_up<-DEGsfilt(DEGs=i, padj=0.05, FC="up")
  
  ego_up[[c]] <- enrichGO(gene          = i_up,
                          universe      = rownames(i),
                          OrgDb         = org.Hs.eg.db,
                          keyType = "SYMBOL",
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)
  ego_dn[[c]] <- enrichGO(gene          = i_down,
                          universe      = rownames(i),
                          OrgDb         = org.Hs.eg.db,
                          keyType = "SYMBOL",
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)
  
}


####GO up
allpaths<-c()
for(i in 1:length(ego_dn)){
  allpaths<-union(allpaths, ego_dn[[i]]$Description[ego_dn[[i]]$p.adjust<0.05])
}
allpaths_mat<-matrix(0,nrow=length(allpaths), ncol=length(ego_dn))
rownames(allpaths_mat)<-allpaths
for(i in 1:length(ego_dn)){
  allpaths_mat[ego_dn[[i]]$Description[which(ego_dn[[i]]$p.adjust<0.05)],i]<-1
}

colnames(allpaths_mat)<-names(ego_dn)

shared_paths<-allpaths_mat[rowSums(allpaths_mat)>3, ]

toplot<-shared_paths
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

#graphics.off()
#pdf(paste("results/",date, "/GO_up_shared.pdf", sep=""), 10, 30)
#pheatmap(toplot,   cellwidth=15, cellheight=15,  keep.dendro=T)
#dev.off()

###################################
############ GSEA
###################################

library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)


fgsea_MsigdbC2CP<-list()
for(c in c( "ddsHs",
  "dds231TFDP1", "dds231E2F3",
  "dds468TFDP1","dds468E2F3")){
  
  i<-get(c)
  
  forgesea<-unlist(i$log2FoldChange)
  names(forgesea)<-rownames(i)
  forgesea<-forgesea[!is.na(forgesea)]
  fgsea_MsigdbC2CP[[c]]<-GSEA(sort(forgesea, decreasing=T), TERM2GENE=m_df, pvalueCutoff = 1, maxGSSize = 10000)
}


allpaths<-c()
for(i in 1:length(fgsea_MsigdbC2CP)){
  allpaths<-union(allpaths, fgsea_MsigdbC2CP[[i]]$Description[fgsea_MsigdbC2CP[[i]]$p.adjust<0.05])
}
allpaths_mat<-matrix(0,nrow=length(allpaths), ncol=length(fgsea_MsigdbC2CP))
rownames(allpaths_mat)<-allpaths
for(i in 1:length(fgsea_MsigdbC2CP)){
  allpaths_mat[fgsea_MsigdbC2CP[[i]]$Description[which(fgsea_MsigdbC2CP[[i]]$p.adjust<0.05)],i]<-fgsea_MsigdbC2CP[[i]]$NES[which(fgsea_MsigdbC2CP[[i]]$p.adjust<0.05)]
}


colnames(allpaths_mat)<-names(fgsea_MsigdbC2CP)

shared_paths<-allpaths_mat[rowSums(allpaths_mat!=0)>0, ]
colnames(shared_paths)<-c("Hs578 TFDP1", "MDAMB231 TFDP1", "MDAMB231 E2F3","MDAMB468 TFDP1", "MDAMB468 E2F3")
  
toplot<-shared_paths
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

graphics.off()
png(paste("results/",date, "/GSEA_MSigDB_Hallmarks_shared.png", sep=""), res=300, 2000, 2000)
pheatmap(toplot,   cellwidth=15, cellheight=15,  keep.dendro=T, breaks = myBreaks, color = myColor)
dev.off()


####################################
####### enrichr
###################################

library(enrichR)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
if (websiteLive) dbs <- listEnrichrDbs()
if (websiteLive) head(dbs)

dbs <- c("GO_Biological_Process_2023","WikiPathways_2024_Human", "Reactome_2022", "TF_Perturbations_Followed_by_Expression", "ENCODE_TF_ChIP-seq_2015")

enriched_up<-list()
enriched_down<-list()
for(c in c(c("ddsHs",
             "dds231TFDP1", "dds231E2F3",
             "dds468TFDP1","dds468E2F3"))){
  
  i<-get(c)
  i_down<-DEGsfilt(DEGs=i, padj=0.05, FC="down")
  i_up<-DEGsfilt(DEGs=i, padj=0.05, FC="up")
  enriched_up[[c]] <- enrichr(i_up, dbs)
  enriched_down[[c]] <- enrichr(i_down, dbs)
}

for(j in 1:length(dbs)){
allpaths<-c()
for(i in 1:length(enriched_down)){
  allpaths<-union(allpaths, enriched_down[[i]][[j]]$Term[enriched_down[[i]][[j]]$Adjusted.P.value<0.05])
}
allpaths_mat<-matrix(0,nrow=length(allpaths), ncol=length(enriched_down))
rownames(allpaths_mat)<-allpaths
for(i in 1:length(enriched_down)){
  allpaths_mat[enriched_down[[i]][[j]]$Term[which(enriched_down[[i]][[j]]$Adjusted.P.value<0.05)],i]<-enriched_down[[i]][[j]]$Combined.Score[which(enriched_down[[i]][[j]]$Adjusted.P.value<0.05)]
}


colnames(allpaths_mat)<-names(enriched_down)

shared_paths<-allpaths_mat[rowSums(allpaths_mat!=0)>2, ]
colnames(shared_paths)<-c("Hs578 TFDP1", "MDAMB231 TFDP1", "MDAMB231 E2F3","MDAMB468 TFDP1", "MDAMB468 E2F3")

toplot<-shared_paths
paletteLength <- 50
myColor <- colorRampPalette(c("white", "red"))(paletteLength)
myBreaks <- c(seq( 0, length.out=ceiling(paletteLength)+1))
length(myBreaks) == length(paletteLength) + 1

graphics.off()
png(paste("results/",date, "/", dbs[j],"_Novogene.png", sep=""), res=300, 4000, nrow(toplot)*100)
pheatmap(toplot,   cellwidth=15, cellheight=15,  keep.dendro=T, breaks = myBreaks, color = myColor)
dev.off()
}


######################################
##### Gene Essentiality
#########################################

###load input data
load("data/Sanger_Broad_higQ_scaled_depFC.RData")
CMP_annot <- read.csv("data/model_list_20210611.csv") # from https://cog.sanger.ac.uk/cmp/download/model_list_20210611.csv

ess<-scaled_depFC[c("E2F3", "TFDP1", "TEAD4", "CEBPG", "PTTG1"),c("MDA-MB-468","MDA-MB-231","Hs-578-T")]
write.xlsx(data.frame(ess), file=paste("results/",date, "/essentiality.xlsx", sep=""), rowNames=T)



