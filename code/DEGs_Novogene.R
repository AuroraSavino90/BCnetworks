library(DESeq2)
library(FactoMineR)
library(ggplot2)
library(openxlsx)
library(pheatmap)
library(WGCNA)
library(ggrepel)
library(ggpubr)

setwd("D:/MBC Dropbox/Lab Poli PhD/Aurora/Projects_wd/BC networks")

date<-"20241008"


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

counts<-read.xlsx("data/gene_count.xlsx")
rownames(counts)<-counts[,1]
anno<-counts[,c("gene_id", "gene_name")]

counts_name<-changenames(counts[,c(2:34)], anno = anno)
RPM<-t(t(counts_name)/colSums(counts_name))*1000000

counts_name<-counts_name[rowSums(counts_name>=10)>2,]
RPM<-RPM[rownames(counts_name),]

RPMlog<-log2(RPM+1)

###########################
### quality checks: PCA
##########################

metadata<-read.xlsx("data/Novogene_metadata.xlsx", rowNames = T)
metadata$Clone<-factor(metadata$Clone)

pca<-PCA(t(RPMlog))
df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2], PC3=pca$ind$coord[,3],
               PC4=pca$ind$coord[,4], PC5=pca$ind$coord[,5],
               metadata)

png(paste("results/",date, "/PCA.png", sep=""), res=300, 1500, 1500)
ggplot(df, aes(x=PC1, y=PC2, colour=Cell.line))+geom_point()
dev.off()

################## PCA separating cell lines

pca231<-PCA(t(RPMlog[,metadata$Cell.line=="MDAMB231"]))
df<-data.frame(PC1=pca231$ind$coord[,1], PC2=pca231$ind$coord[,2], PC3=pca231$ind$coord[,3],
               PC4=pca231$ind$coord[,4], PC5=pca231$ind$coord[,5],
               metadata[metadata$Cell.line=="MDAMB231",])

png(paste("results/",date, "/PCA231.png", sep=""), res=300, 1500, 1500)
ggplot(df, aes(x=PC1, y=PC2, colour=KO.gene, shape=as.factor(Clone)))+geom_point()
dev.off()

pca468<-PCA(t(RPMlog[,metadata$Cell.line=="MDAMB468"]))
df<-data.frame(PC1=pca468$ind$coord[,1], PC2=pca468$ind$coord[,2], PC3=pca468$ind$coord[,3],
               PC4=pca468$ind$coord[,4], PC5=pca468$ind$coord[,5],
               metadata[metadata$Cell.line=="MDAMB468",])

png(paste("results/",date, "/PCA468.png", sep=""), res=300, 1500, 1500)
ggplot(df, aes(x=PC1, y=PC2, colour=KO.gene, shape=as.factor(Clone)))+geom_point()
dev.off()

pcaHs<-PCA(t(RPMlog[,metadata$Cell.line=="Hs578"]))
df<-data.frame(PC1=pcaHs$ind$coord[,1], PC2=pcaHs$ind$coord[,2], PC3=pcaHs$ind$coord[,3],
               PC4=pcaHs$ind$coord[,4], PC5=pcaHs$ind$coord[,5],
               metadata[metadata$Cell.line=="Hs578",])

png(paste("results/",date, "/PCAHs.png", sep=""), res=300, 1500, 1500)
ggplot(df, aes(x=PC1, y=PC2, colour=KO.gene, shape=as.factor(Clone)))+geom_point()
dev.off()


#######################################################
########### DEGs
#########################################################

#clone can't be used as a covariate as we have only one clone for controls, leading to nested variables

############## pooling clones

dds468 <- DESeqDataSetFromMatrix(countData = counts_name[,metadata$Cell.line=="MDAMB468"],
                                   colData = metadata[metadata$Cell.line=="MDAMB468",],
                                   design= ~ KO.gene)
dds468 <- DESeq(dds468)
dds468 <- results(dds468, contrast = c("KO.gene","TFDP1", "EV"))

ddsHs <- DESeqDataSetFromMatrix(countData = counts_name[,metadata$Cell.line=="Hs578"],
                                   colData = metadata[metadata$Cell.line=="Hs578",],
                                   design= ~ KO.gene)
ddsHs <- DESeq(ddsHs)
ddsHs <- results(ddsHs, contrast = c("KO.gene","TFDP1", "EV"))

dds231 <- DESeqDataSetFromMatrix(countData = counts_name[,metadata$Cell.line=="MDAMB231"],
                                    colData = metadata[metadata$Cell.line=="MDAMB231",],
                                    design= ~ KO.gene)
dds231 <- DESeq(dds231)
dds231TFDP1 <- results(dds231, contrast = c("KO.gene","TFDP1", "EV"))
dds231E2F3 <- results(dds231, contrast = c("KO.gene","E2F3", "EV"))


################# each clone separately

dds468_4 <- DESeqDataSetFromMatrix(countData = counts_name[,metadata$Cell.line=="MDAMB468" & metadata$Clone %in% c("100", "4")],
                                 colData = metadata[metadata$Cell.line=="MDAMB468" & metadata$Clone %in% c("100", "4"),],
                                 design= ~ KO.gene)
dds468_4 <- DESeq(dds468_4)
dds468_4 <- results(dds468_4, contrast = c("KO.gene","TFDP1", "EV"))

dds468_9 <- DESeqDataSetFromMatrix(countData = counts_name[,metadata$Cell.line=="MDAMB468" & metadata$Clone %in% c("100", "9")],
                                   colData = metadata[metadata$Cell.line=="MDAMB468" & metadata$Clone %in% c("100", "9"),],
                                   design= ~ KO.gene)
dds468_9 <- DESeq(dds468_9)
dds468_9 <- results(dds468_9, contrast = c("KO.gene","TFDP1", "EV"))


ddsHs_14 <- DESeqDataSetFromMatrix(countData = counts_name[,metadata$Cell.line=="Hs578" & metadata$Clone %in% c("100", "14")],
                                colData = metadata[metadata$Cell.line=="Hs578" & metadata$Clone %in% c("100", "14"),],
                                design= ~ KO.gene)
ddsHs_14 <- DESeq(ddsHs_14)
ddsHs_14 <- results(ddsHs_14, contrast = c("KO.gene","TFDP1", "EV"))

ddsHs_18 <- DESeqDataSetFromMatrix(countData = counts_name[,metadata$Cell.line=="Hs578" & metadata$Clone %in% c("100", "18")],
                                   colData = metadata[metadata$Cell.line=="Hs578" & metadata$Clone %in% c("100", "18"),],
                                   design= ~ KO.gene)
ddsHs_18 <- DESeq(ddsHs_18)
ddsHs_18 <- results(ddsHs_18, contrast = c("KO.gene","TFDP1", "EV"))

dds231_18 <- DESeqDataSetFromMatrix(countData = counts_name[,metadata$Cell.line=="MDAMB231" & metadata$Clone %in% c("100", "18")],
                                 colData = metadata[metadata$Cell.line=="MDAMB231" & metadata$Clone %in% c("100", "18"),],
                                 design= ~ KO.gene)
dds231_18 <- DESeq(dds231_18)
dds231_18 <- results(dds231_18, contrast = c("KO.gene","E2F3", "EV"))

dds231_20 <- DESeqDataSetFromMatrix(countData = counts_name[,metadata$Cell.line=="MDAMB231" & metadata$Clone %in% c("100", "20")],
                                    colData = metadata[metadata$Cell.line=="MDAMB231" & metadata$Clone %in% c("100", "20"),],
                                    design= ~ KO.gene)
dds231_20 <- DESeq(dds231_20)
dds231_20 <- results(dds231_20, contrast = c("KO.gene","E2F3", "EV"))

dds231_6 <- DESeqDataSetFromMatrix(countData = counts_name[,metadata$Cell.line=="MDAMB231" & metadata$Clone %in% c("100", "6")],
                                    colData = metadata[metadata$Cell.line=="MDAMB231" & metadata$Clone %in% c("100", "6"),],
                                    design= ~ KO.gene)
dds231_6 <- DESeq(dds231_6)
dds231_6 <- results(dds231_6, contrast = c("KO.gene","TFDP1", "EV"))


dds231_21 <- DESeqDataSetFromMatrix(countData = counts_name[,metadata$Cell.line=="MDAMB231" & metadata$Clone %in% c("100", "21")],
                                    colData = metadata[metadata$Cell.line=="MDAMB231" & metadata$Clone %in% c("100", "21"),],
                                    design= ~ KO.gene)
dds231_21 <- DESeq(dds231_21)
dds231_21 <- results(dds231_21, contrast = c("KO.gene","TFDP1", "EV"))

#################################
### DEGs overlap
##################################
DEGs_shared<-cbind(ddsHs$log2FoldChange, dds468$log2FoldChange, dds231TFDP1$log2FoldChange, dds231E2F3$log2FoldChange)
DEGs_shared[which(ddsHs$padj>0.05),1]<-0
DEGs_shared[which(dds468$padj>0.05),2]<-0
DEGs_shared[which(dds231TFDP1$padj>0.05),3]<-0
DEGs_shared[which(dds231E2F3$padj>0.05),4]<-0
rownames(DEGs_shared)<-rownames(ddsHs)
colnames(DEGs_shared)<-c("Hs578 TFDP1", "MDAMB468 TFDP1", "MDAMB231 TFDP1", "MDAMB231 E2F3")

toplot<-DEGs_shared[DEGs_shared[,1]!=0,]
toplot<-toplot[-which(rowSums(is.na(toplot))==4),]
toplot<-toplot[which(rownames(toplot) %in% rownames(centrality_basal)[which(centrality_basal$module=="b_E2F_targets")]),]
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

######network data to load
load("../../R analyses/METABRIC networks/centrality_basal.RData")
load("../../R analyses/METABRIC networks/metabric.RData")
load("../../R analyses/METABRIC networks/meta.RData")

############## pooling clones

ft<-fishertest_alldat(alldat=c("dds468", "ddsHs","dds231TFDP1", "dds231E2F3"),
                  names_alldat=c("MDAMB468 TFDP1", "Hs578 TFDP1", "MDAMB231 TFDP1", "MDAMB231 E2F3"))

graphics.off()
png(paste("results/",date, "/Enrich_up_Novogene.png", sep=""), res=300, 1500, 2500)
pheatmap(-log10(ft[[1]]),cellwidth=15, cellheight=15)
dev.off()

png(paste("results/",date, "/Enrich_down_Novogene.png", sep=""), res=300, 1500, 2500)
pheatmap(-log10(ft[[2]]),cellwidth=15, cellheight=15)
dev.off()

png(paste("results/",date, "/Enrich_all_Novogene.png", sep=""), res=300, 1500, 2500)
pheatmap(-log10(ft[[3]]),cellwidth=15, cellheight=15)
dev.off()

############ each clone separately

ft2<-fishertest_alldat(alldat=c("dds468_4","dds468_9",
                                "ddsHs_14","ddsHs_18",
                                "dds231_6","dds231_21",
                                "dds231_18", "dds231_20"),
                      names_alldat=c("MDAMB468 TFDP1 4", "MDAMB468 TFDP1 9", 
                                     "Hs578 TFDP1 14","Hs578 TFDP1 18",
                                     "MDAMB231 TFDP1 6","MDAMB231 TFDP1 21",
                                     "MDAMB231 E2F3 18", "MDAMB231 E2F3 20"))


graphics.off()
png(paste("results/",date, "/Enrich_up_Novogene_clones.png", sep=""), res=300, 2500, 2500)
pheatmap(-log10(ft2[[1]]),cellwidth=15, cellheight=15)
dev.off()

png(paste("results/",date, "/Enrich_down_Novogene_clones.png", sep=""), res=300, 2500, 2500)
pheatmap(-log10(ft2[[2]]),cellwidth=15, cellheight=15)
dev.off()

png(paste("results/",date, "/Enrich_all_Novogene_clones.png", sep=""), res=300, 2500, 2500)
pheatmap(-log10(ft2[[3]]),cellwidth=15, cellheight=15)
dev.off()



###########################
##correlation between modules
##########################

MEs<-moduleEigengenes(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]), centrality_basal$module)$eigengenes
colnames(MEs)<-gsub("ME", "", colnames(MEs))
cc<-cor(MEs)
cc<-cc[,-20]

#top altered (in significance) are the most positively and negatively correlated modules with E2F_targets
plot(cc["b_E2F_targets",], rowSums(-log10(ft[[1]][names(cc["b_E2F_targets",]),]))+rowSums(-log10(ft[[2]][names(cc["b_E2F_targets",]),])))
plot(cc["b_E2F_targets",], rowSums(-log10(ft_o[[1]][names(cc["b_E2F_targets",]),]))+rowSums(-log10(ft_o[[2]][names(cc["b_E2F_targets",]),])))


#########################################################
### projection of b_E2F_targets MEs on new data  
#########################################################

#####project MEs on Novogene data

pcaproj_hubs<-pcaproject(newdata=RPMlog, original_data=metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], modules=colnames(MEs), ME=MEs)
colnames(pcaproj_hubs)<-colnames(MEs)
rownames(pcaproj_hubs)<-colnames(RPMlog)

pcaproj_hubs<-cbind.data.frame(pcaproj_hubs, metadata)
pcaproj_hubs$KO.gene<-factor(pcaproj_hubs$KO.gene, levels=c("EV", "TFDP1", "E2F3"))

####plot b_E2F_targets across KOs
png(paste("results/",date, "/bE2F_project_Novogene.png", sep=""), res=300, 1500, 1500)
ggplot(pcaproj_hubs, aes(x=Cell.line, y=b_E2F_targets, fill=KO.gene))+geom_boxplot()+theme_classic()
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

ctrl<-subset(pcaproj_hubs, Cell.line=="MDAMB231" &  KO.gene=="EV")
trt<-subset(pcaproj_hubs, Cell.line=="MDAMB231" &  KO.gene=="E2F3")

cd2<-matrix(nrow=1, ncol=19)
ind<-1

for(i in 1:19){
  cd2[ind,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}
colnames(cd2)<-colnames(pcaproj_hubs)[1:19]
rownames(cd2)<-"MDAMB231 E2F3"

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

########## plot changes in MEs (Cohen's d) vs modules' correlation with b_E2F_targets
## for each KO, the most affected modules are either the most highly or lowly correlated with b_E2F_targets

df<-data.frame(cd=c(t(cd_all)), condition=rep(c("MDAMB468 TFDP1", "Hs578 TFDP1", "MDAMB231 TFDP1", "MDAMB231 E2F3"), each= 19), corr=rep(cc["b_E2F_targets",colnames(cd_all)],4),
               module=rep(colnames(cd_all),4))
df$condition<-factor(df$condition, levels=c("Hs578 TFDP1","MDAMB468 TFDP1", "MDAMB231 TFDP1", "MDAMB231 E2F3"))

png(paste("results/",date, "/CohenVScorr_bE2F_Novogene.png", sep=""),res=300, 4500, 2000)
ggplot(df, aes(x=corr, y=cd, label=module))+geom_point(size=2)+facet_grid(~condition)+geom_smooth(method = lm)+stat_cor(label.x=-0.5, label.y = 17)+geom_text_repel(max.overlaps = 5)+theme_bw()+theme(strip.text=element_text(size = 12, face = "bold"))
dev.off()


######for each clone separately
pcaproj_hubs$condition<-factor(paste(paste(pcaproj_hubs$Cell.line, pcaproj_hubs$Clone, sep=" "), pcaproj_hubs$Clone, sep=" "),
                               levels=c("Hs578 EV 100", "Hs578 TFDP1 14", "Hs578 TFDP1 18",
                                        "MDAMB231 EV 100",      "MDAMB231 TFDP1 6",  "MDAMB231 TFDP1 21",
                                        "MDAMB231 E2F3 18","MDAMB231 E2F3 20",
                                        "MDAMB468 EV 100",   "MDAMB468 TFDP1 4",  "MDAMB468 TFDP1 9"))

####plot b_E2F_targets across KOs
png(paste("results/",date, "/bE2F_project_Novogene_clones.png", sep=""), res=300, 2000, 3000)
ggplot(pcaproj_hubs, aes(x=Cell.line, y=b_E2F_targets, fill=condition))+geom_boxplot()+theme_classic()+
  scale_fill_manual(values=c("green", "red", "orange",
                        "green", "red", "orange", "blue", "lightblue",
                        "green", "red", "orange"))
dev.off()

cd_sep<-matrix(nrow=8, ncol=19)

ctrl<-subset(pcaproj_hubs, Cell.line=="MDAMB468" &  KO.gene=="EV")
trt<-subset(pcaproj_hubs, Cell.line=="MDAMB468" &  KO.gene=="TFDP1" & Clone==4)

for(i in 1:19){
  cd_sep[1,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}

ctrl<-subset(pcaproj_hubs, Cell.line=="MDAMB468" &  KO.gene=="EV")
trt<-subset(pcaproj_hubs, Cell.line=="MDAMB468" &  KO.gene=="TFDP1" & Clone==9)

for(i in 1:19){
  cd_sep[2,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}

ctrl<-subset(pcaproj_hubs, Cell.line=="Hs578" &  KO.gene=="EV")
trt<-subset(pcaproj_hubs, Cell.line=="Hs578" &  KO.gene=="TFDP1" & Clone==14)

for(i in 1:19){
  cd_sep[3,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}

ctrl<-subset(pcaproj_hubs, Cell.line=="Hs578" &  KO.gene=="EV")
trt<-subset(pcaproj_hubs, Cell.line=="Hs578" &  KO.gene=="TFDP1" & Clone==18)

for(i in 1:19){
  cd_sep[4,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}

ctrl<-subset(pcaproj_hubs, Cell.line=="MDAMB231" &  KO.gene=="EV")
trt<-subset(pcaproj_hubs, Cell.line=="MDAMB231" &  KO.gene=="E2F3" & Clone==18)

for(i in 1:19){
  cd_sep[5,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}

ctrl<-subset(pcaproj_hubs, Cell.line=="MDAMB231" &  KO.gene=="EV")
trt<-subset(pcaproj_hubs, Cell.line=="MDAMB231" &  KO.gene=="E2F3" & Clone==20)

for(i in 1:19){
  cd_sep[6,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}

ctrl<-subset(pcaproj_hubs, Cell.line=="MDAMB231" &  KO.gene=="EV")
trt<-subset(pcaproj_hubs, Cell.line=="MDAMB231" &  KO.gene=="TFDP1" & Clone==6)

for(i in 1:19){
  cd_sep[7,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}

ctrl<-subset(pcaproj_hubs, Cell.line=="MDAMB231" &  KO.gene=="EV")
trt<-subset(pcaproj_hubs, Cell.line=="MDAMB231" &  KO.gene=="TFDP1" & Clone==21)

for(i in 1:19){
  cd_sep[8,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}

colnames(cd_sep)<-colnames(pcaproj_hubs)[1:19]
rownames(cd_sep)<-c("MDAMB468 TFDP1 4", "MDAMB468 TFDP1 9", 
  "Hs578 TFDP1 14","Hs578 TFDP1 18",
  "MDAMB231 E2F3 18", "MDAMB231 E2F3 20",
  "MDAMB231 TFDP1 6","MDAMB231 TFDP1 21")

########## plot changes in MEs (Cohen's d) for each clone
## Clones are coherent (with some outliers)

toplot<-cd_sep
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

graphics.off()
png(paste("results/",date, "/Cohen_bE2F_Novogene_clones.png", sep=""),res=300, 1500, 2500)
pheatmap(t(toplot),   cellwidth=15, cellheight=15,  keep.dendro=T, color =myColor, breaks = myBreaks)
dev.off()

########## plot changes in MEs (Cohen's d) vs modules' correlation with b_E2F_targets
## for each KO, the most affected modules are either the most highly or lowly correlated with b_E2F_targets

df<-data.frame(cd=c(t(cd_sep)), condition=rep(c("MDAMB468 TFDP1 4", "MDAMB468 TFDP1 9", 
                                                "Hs578 TFDP1 14","Hs578 TFDP1 18",
                                                "MDAMB231 E2F3 18", "MDAMB231 E2F3 20",
                                                "MDAMB231 TFDP1 6","MDAMB231 TFDP1 21"), each= 19), corr=rep(cc["b_E2F_targets",colnames(cd_sep)],8),
               module=rep(colnames(cd_sep),8))
df$condition<-factor(df$condition, levels=c("Hs578 TFDP1 14","Hs578 TFDP1 18",
                                            "MDAMB468 TFDP1 4", "MDAMB468 TFDP1 9",
                                            "MDAMB231 TFDP1 6","MDAMB231 TFDP1 21",
                                            "MDAMB231 E2F3 18", "MDAMB231 E2F3 20"))

png(paste("results/",date, "/CohenVScorr_bE2F_Novogene_clones.png", sep=""),res=300, 6000, 2000)
ggplot(df, aes(x=corr, y=cd, label=module))+geom_point(size=2)+facet_grid(~condition)+geom_smooth(method = lm)+stat_cor(label.x=-0.5, label.y = 30)+geom_text_repel(max.overlaps = 5)+
  scale_y_continuous(limits = c(-30, 30))+
  theme_bw()+theme(strip.text=element_text(size = 12, face = "bold"))
dev.off()



######################################
######## GO enrichment
######################################

library(clusterProfiler)
library(org.Hs.eg.db)

ego_up<-list()
ego_dn<-list()
for(c in c("dds468",
             "ddsHs",
             "dds231TFDP1", "dds231E2F3")){
  
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

shared_paths<-allpaths_mat[rowSums(allpaths_mat)>2, ]

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
for(c in c("dds468",
           "ddsHs",
           "dds231TFDP1", "dds231E2F3")){
  
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
colnames(shared_paths)<-c("MDAMB468 TFDP1", "Hs578 TFDP1", "MDAMB231 TFDP1", "MDAMB231 E2F3")
  
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

############# each clone separately

fgsea_MsigdbC2CP_cl<-list()
for(c in c("dds468_4","dds468_9",
           "ddsHs_14","ddsHs_18",
           "dds231_6","dds231_21",
           "dds231_18", "dds231_20")){
  
  i<-get(c)
  
  forgesea<-unlist(i$log2FoldChange)
  names(forgesea)<-rownames(i)
  forgesea<-forgesea[!is.na(forgesea)]
  fgsea_MsigdbC2CP_cl[[c]]<-GSEA(sort(forgesea, decreasing=T), TERM2GENE=m_df, pvalueCutoff = 1, maxGSSize = 10000)
}


allpaths<-c()
for(i in 1:length(fgsea_MsigdbC2CP_cl)){
  allpaths<-union(allpaths, fgsea_MsigdbC2CP_cl[[i]]$Description[fgsea_MsigdbC2CP_cl[[i]]$p.adjust<0.05])
}
allpaths_mat<-matrix(0,nrow=length(allpaths), ncol=length(fgsea_MsigdbC2CP_cl))
rownames(allpaths_mat)<-allpaths
for(i in 1:length(fgsea_MsigdbC2CP_cl)){
  allpaths_mat[fgsea_MsigdbC2CP_cl[[i]]$Description[which(fgsea_MsigdbC2CP_cl[[i]]$p.adjust<0.05)],i]<-fgsea_MsigdbC2CP_cl[[i]]$NES[which(fgsea_MsigdbC2CP_cl[[i]]$p.adjust<0.05)]
}


colnames(allpaths_mat)<-names(fgsea_MsigdbC2CP_cl)

shared_paths<-allpaths_mat[rowSums(allpaths_mat!=0)>0, ]
colnames(shared_paths)<-c("MDAMB468 TFDP1 4", "MDAMB468 TFDP1 9", 
                          "Hs578 TFDP1 14","Hs578 TFDP1 18",
                          "MDAMB231 TFDP1 6","MDAMB231 TFDP1 21",
                          "MDAMB231 E2F3 18", "MDAMB231 E2F3 20")

toplot<-shared_paths
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

graphics.off()
png(paste("results/",date, "/GSEA_MSigDB_Hallmarks_shared_clones.png", sep=""), res=300, 2500, 2500)
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
enriched <- enrichr(DE_common, dbs)
plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

enriched_up<-list()
enriched_down<-list()
for(c in c(c("dds468",
             "ddsHs",
             "dds231TFDP1", "dds231E2F3"))){
  
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
colnames(shared_paths)<-c("MDAMB468 TFDP1", "Hs578 TFDP1", "MDAMB231 TFDP1", "MDAMB231 E2F3")

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


enriched_up_cl<-list()
enriched_down_cl<-list()
for(c in c("dds468_4","dds468_9",
           "ddsHs_14","ddsHs_18",
           "dds231_6","dds231_21",
           "dds231_18", "dds231_20")){
  
  i<-get(c)
  i_down<-DEGsfilt(DEGs=i, padj=0.05, FC="down")
  i_up<-DEGsfilt(DEGs=i, padj=0.05, FC="up")
  enriched_up_cl[[c]] <- enrichr(i_up, dbs)
  enriched_down_cl[[c]] <- enrichr(i_down, dbs)
}


#############################################################################################
###############################################################################################

###########################################
### Old sequencing data (shRNA)
############################################

setwd("D:/MBC Dropbox/Lab Poli PhD/Aurora/R analyses/METABRIC networks/Hubs functional validations/RNASeq hubs/counts")
files<-list.files()

data<-read.csv(files[1], sep="\t", header=F)
for(i in 2:length(files)){
  data<-cbind(data,read.csv(files[i], sep="\t", header=F)[,2])
}

rownames(data)<-data[,1]
data<-data[,-1]
colnames(data)<-files

data_info<-data[c(60677:60681),]
data<-data[-c(60677:60681),]

library("biomaRt")
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
geni_coord<-getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                  values = rownames(data), 
                  mart = ensembl)

anno<-geni_coord[match(rownames(data), geni_coord[,1]),]

data_name<-changenames(data, anno = anno)

RPM_old<-t(t(data_name)/colSums(data_name))*1000000

data_name<-data_name[rowSums(data_name>=3)>2,]
RPM_old<-RPM_old[rownames(data_name),]

RPM_oldlog<-log2(RPM_old+1)



coldata<-data.frame(treatment=c("C2","C2","C2","C5","C5","C5", "E5","E5","E5",
                                "E6", "E6", "E6", "P1","P1","P1","P4","P4","P4",
                                "PSC","PSC", "PSC", "PSE", "PSE", "PSE",
                                "SCR2","SCR2","SCR2", "SCR3", "SCR3", "SCR3",
                                "SCR4","SCR4","SCR4", "SCR5", "SCR5", "SCR5",
                                "T1","T1","T1", "T4", "T4", "T4", "TF4", "TF4", "TF4", "TF5", "TF5", "TF5"),
                    treatment2=c("CEBPG","CEBPG","CEBPG","CEBPG","CEBPG","CEBPG",
                                 "E2F3","E2F3","E2F3","E2F3","E2F3","E2F3",
                                 "PTTG1","PTTG1","PTTG1","PTTG1","PTTG1","PTTG1",
                                 "PSAT1","PSAT1","PSAT1","PSAT1","PSAT1","PSAT1",
                                 "CTRL","CTRL","CTRL","CTRL","CTRL","CTRL",
                                 "CTRL","CTRL","CTRL","CTRL","CTRL","CTRL",
                                 "TEAD4","TEAD4","TEAD4", "TEAD4", "TEAD4", "TEAD4",
                                 "TFDP1","TFDP1","TFDP1","TFDP1","TFDP1","TFDP1" ),
                    replicate=as.factor(c(1,2,3, 1,2,3, 1,2,3, 1,2,3, 1,2,3, 1,2,3, 1,2,3, 1,2,3, 1,2,3, 1,2,3, 1,2,3, 1,2,3, 1,2,3, 1,2,3, 1,2,3, 1,2,3 )))


#pcaproj_hubs<-pcaproject(newdata=data_filt_norm, original_data=metabric, modules=modules, ME=ME)
#colnames(pcaproj_hubs)<-unique(modules)
#rownames(pcaproj_hubs)<-colnames(data_filt_norm)

#incommon<-intersect(rownames(data_filt_norm), rownames(centrality_basal))
#centrality_basal_hubs<-centrality_basal[incommon,]
#data_filt_norm_sel<-data_filt_norm[incommon,]

#ssmod_hubs<-ssMod(centrality_basal=centrality_basal_hubs, data=data_filt_norm_sel)

##PCA on normalized data
setwd("D:/MBC Dropbox/Lab Poli PhD/Aurora/Projects_wd/BC networks")
pca<-PCA(t(RPM_oldlog))
df<-data.frame(PC1=pca$ind$coord[,1], PC2=pca$ind$coord[,2], PC3=pca$ind$coord[,3],
               PC4=pca$ind$coord[,4], PC5=pca$ind$coord[,5],
               coldata)

png(paste("results/",date, "/PCA_sh.png", sep=""), res=300, 1500, 1500)
ggplot(df, aes(x=PC1, y=PC2, colour=treatment2))+geom_point()
dev.off()


##control KD genes' expression
df<-data.frame(E2F3=RPM_oldlog["E2F3",],TFDP1=RPM_oldlog["TFDP1",],PSAT1=RPM_oldlog["PSAT1",],
               PTTG1=RPM_oldlog["PTTG1",],CEBPG=RPM_oldlog["CEBPG",],TEAD4=RPM_oldlog["TEAD4",],
               shRNA=coldata$treatment, gene=coldata$treatment2)

png(paste("results/",date, "/E2F3_sh.png", sep=""), res=300, 2500, 1000)
ggplot(df, aes(x=shRNA, y=E2F3, fill=gene))+geom_boxplot()+theme_bw()
dev.off()
png(paste("results/",date, "/PTTG1_sh.png", sep=""), res=300, 2500, 1000)
ggplot(df, aes(x=shRNA, y=PTTG1, fill=gene))+geom_boxplot()+theme_bw()
dev.off()
png(paste("results/",date, "/TFDP1_sh.png", sep=""), res=300, 2500, 1000)
ggplot(df, aes(x=shRNA, y=TFDP1, fill=gene))+geom_boxplot()+theme_bw()
dev.off()
png(paste("results/",date, "/TEAD4_sh.png", sep=""), res=300, 2500, 1000)
ggplot(df, aes(x=shRNA, y=TEAD4, fill=gene))+geom_boxplot()+theme_bw()
dev.off()
png(paste("results/",date, "/CEBPG_sh.png", sep=""), res=300, 2500, 1000)
ggplot(df, aes(x=shRNA, y=CEBPG, fill=gene))+geom_boxplot()+theme_bw()
dev.off()


###########
## DEGs
###########

coldata$treatment2<-factor(as.factor(coldata$treatment2), levels=c("CTRL", "CEBPG", "TEAD4", "E2F3", "PTTG1", "TFDP1", "PSAT1"))
coldata$treatment<-factor(as.factor(coldata$treatment), levels=c("SCR2", "SCR3","SCR4","SCR5", "C2", "C5", "E5","E6","P1","P4","T1", "T4", "TF4", "TF5", "PSC", "PSE"))
dds <- DESeqDataSetFromMatrix(countData = data_name,
                              colData = coldata,
                              design= ~ treatment2)

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

DE_CEBPG<-results(dds, contrast = c("treatment2", "CEBPG", "CTRL"))
DE_TEAD4<-results(dds, contrast = c("treatment2", "TEAD4", "CTRL"))
DE_PTTG1<-results(dds, contrast = c("treatment2", "PTTG1",  "CTRL"))
DE_E2F3<-results(dds, contrast = c("treatment2",  "E2F3",  "CTRL"))
DE_PSAT1<-results(dds, contrast = c("treatment2", "PSAT1",  "CTRL"))
DE_TFDP1<-results(dds, contrast = c("treatment2", "TFDP1",  "CTRL"))

ft_o<-fishertest_alldat(alldat=c("DE_CEBPG",
                                 "DE_TEAD4",
                                 "DE_PTTG1", "DE_E2F3", "DE_TFDP1", "DE_PSAT1"),
                        names_alldat=c("CEBPG", "TEAD4", "PTTG1", "E2F3", "TFDP1", "PSAT1"))


graphics.off()
png(paste("results/",date, "/Enrich_up_sh.png", sep=""), res=300, 1500, 2500)
pheatmap(-log10(ft_o[[1]]),cellwidth=15, cellheight=15)
dev.off()

png(paste("results/",date, "/Enrich_down_sh.png", sep=""), res=300, 1500, 2500)
pheatmap(-log10(ft_o[[2]]),cellwidth=15, cellheight=15)
dev.off()

png(paste("results/",date, "/Enrich_all_sh.png", sep=""), res=300, 1500, 2500)
pheatmap(-log10(ft_o[[3]]),cellwidth=15, cellheight=15)
dev.off()


##facendo ogni sh singolarmente con tutti i controlli

dds <- DESeqDataSetFromMatrix(countData = data_name,
                              colData = coldata,
                              design= ~ treatment)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds)


DE_CEBPG_2<-results(dds, contrast = c("treatment","C2", "SCR2"))
DE_CEBPG_5<-results(dds, contrast = c("treatment", "C5", "SCR2"))

DE_TEAD4_1<-results(dds, contrast = c("treatment",  "T1", "SCR2"))
DE_TEAD4_4<-results(dds, contrast = c("treatment",  "T4", "SCR2"))

DE_PTTG1_1<-results(dds, contrast = c("treatment",  "P1", "SCR4"))
DE_PTTG1_4<-results(dds, contrast = c("treatment",  "P4", "SCR4"))

DE_E2F3_5<-results(dds, contrast = c("treatment","E5", "SCR4"))
DE_E2F3_6<-results(dds, contrast = c("treatment", "E6", "SCR5"))

DE_PSAT1_C<-results(dds, contrast = c("treatment",  "PSC", "SCR5"))
DE_PSAT1_E<-results(dds, contrast = c("treatment",  "PSE", "SCR5"))

DE_TFDP1_4<-results(dds, contrast = c("treatment",  "TF4", "SCR3"))
DE_TFDP1_5<-results(dds, contrast = c("treatment",  "TF5", "SCR3"))


ft_o_sep<-fishertest_alldat(alldat=c("DE_CEBPG_2","DE_CEBPG_5",
                                     "DE_TEAD4_1", "DE_TEAD4_4",
                                     "DE_PTTG1_1","DE_PTTG1_4",
                                     "DE_E2F3_5","DE_E2F3_6", 
                                     "DE_TFDP1_4", "DE_TFDP1_5"),
                            names_alldat=c("CEBPG 2","CEBPG 5", "TEAD4 1","TEAD4 4",
                                           "PTTG1 1","PTTG1 4", 
                                           "E2F3 5", "E2F3 6", 
                                           "TFDP1 4", "TFDP1 5"))


graphics.off()
png(paste("results/",date, "/Enrich_up_sh_sep.png", sep=""), res=300, 1500, 2500)
pheatmap(-log10(ft_o_sep[[1]]),cellwidth=15, cellheight=15)
dev.off()

png(paste("results/",date, "/Enrich_down_sh_sep.png", sep=""), res=300, 1500, 2500)
pheatmap(-log10(ft_o_sep[[2]]),cellwidth=15, cellheight=15)
dev.off()

png(paste("results/",date, "/Enrich_all_sh_sep.png", sep=""), res=300, 1500, 2500)
pheatmap(-log10(ft_o_sep[[3]]),cellwidth=15, cellheight=15)
dev.off()

plot(c(DE_CEBPG_2["CEBPG","log2FoldChange"], DE_CEBPG_5["CEBPG","log2FoldChange"],
       DE_TEAD4_1["TEAD4", "log2FoldChange"], DE_TEAD4_4["TEAD4", "log2FoldChange"],
       DE_PTTG1_1["PTTG1","log2FoldChange"], DE_PTTG1_4["PTTG1","log2FoldChange"], 
       DE_E2F3_5["E2F3","log2FoldChange"], DE_E2F3_6["E2F3","log2FoldChange"],
       DE_TFDP1_4["TFDP1","log2FoldChange"], DE_TFDP1_5["TFDP1","log2FoldChange"]), -log10(ft_o_sep[[2]])["b_E2F_targets",])

df<-cbind.data.frame(-log10(t(ft_o_sep[[2]])), TF_log2FC=c(DE_CEBPG_2["CEBPG","log2FoldChange"], DE_CEBPG_5["CEBPG","log2FoldChange"],
                                                           DE_TEAD4_1["TEAD4", "log2FoldChange"], DE_TEAD4_4["TEAD4", "log2FoldChange"],
                                                           DE_PTTG1_1["PTTG1","log2FoldChange"], DE_PTTG1_4["PTTG1","log2FoldChange"], 
                                                           DE_E2F3_5["E2F3","log2FoldChange"], DE_E2F3_6["E2F3","log2FoldChange"],
                                                           DE_TFDP1_4["TFDP1","log2FoldChange"], DE_TFDP1_5["TFDP1","log2FoldChange"]),
                     TF=c("CEBPG 2", "CEBPG 5", "TEAD4 1", "TEAD4 4", "PTTG1 1", "PTTG1 4", "E2F3 5", "E2F3 6", "TFDP1 4", "TFDP1 5"))

png(paste("results/",date, "/EnrichVSlog2FC_sh_sep.png", sep=""),res=300, 1500, 1500)
ggplot(df, aes(x=TF_log2FC, y=b_E2F_targets, label=TF))+geom_point()+geom_label_repel()+theme_bw()
dev.off()


df<-cbind.data.frame(-log10(t(ft_o_sep[[1]])), TF_log2FC=c(DE_CEBPG_2["CEBPG","log2FoldChange"], DE_CEBPG_5["CEBPG","log2FoldChange"],
                                                           DE_TEAD4_1["TEAD4", "log2FoldChange"], DE_TEAD4_4["TEAD4", "log2FoldChange"],
                                                           DE_PTTG1_1["PTTG1","log2FoldChange"], DE_PTTG1_4["PTTG1","log2FoldChange"], 
                                                           DE_E2F3_5["E2F3","log2FoldChange"], DE_E2F3_6["E2F3","log2FoldChange"],
                                                           DE_TFDP1_4["TFDP1","log2FoldChange"], DE_TFDP1_5["TFDP1","log2FoldChange"]),
                     TF=c("CEBPG 2", "CEBPG 5", "TEAD4 1", "TEAD4 4", "PTTG1 1", "PTTG1 4", "E2F3 5", "E2F3 6", "TFDP1 4", "TFDP1 5"))

ggplot(df, aes(x=TF_log2FC, y=b_E2F_targets, label=TF))+geom_point()+geom_label_repel()+theme_bw()

#######################
## putting all data together
#######################
anno<-data.frame(batch=c(rep("KO", ncol(ft[[1]])), rep("sh", ncol(ft_o[[1]]))))
rownames(anno)<-colnames(cbind(-log10(ft[[1]]), -log10(ft_o[[1]])))

graphics.off()
png(paste("results/",date, "/Enrich_up_alldata.png", sep=""), res=300, 2500, 2500)
pheatmap(cbind(-log10(ft[[1]]), -log10(ft_o[[1]])),cellwidth=15, cellheight=15, annotation_col = anno)
dev.off()

png(paste("results/",date, "/Enrich_down_alldata.png", sep=""), res=300, 2500, 2500)
pheatmap(cbind(-log10(ft[[2]]), -log10(ft_o[[2]])),cellwidth=15, cellheight=15, annotation_col = anno)
dev.off()

png(paste("results/",date, "/Enrich_all_alldata.png", sep=""), res=300, 2500, 2500)
pheatmap(cbind(-log10(ft[[3]]), -log10(ft_o[[3]])),cellwidth=15, cellheight=15, annotation_col = anno)
dev.off()


##########################################
#####project MEs on sh data
###########################################

pcaproj_hubs_o<-pcaproject(newdata=RPM_oldlog, original_data=metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], modules=colnames(MEs), ME=MEs)
colnames(pcaproj_hubs_o)<-colnames(MEs)
rownames(pcaproj_hubs_o)<-colnames(RPM_oldlog)
pheatmap(t(pcaproj_hubs_o))

pcaproj_hubs_o<-cbind.data.frame(pcaproj_hubs_o, coldata, t(RPM_oldlog))

png(paste("results/",date, "/bE2F_project_sh.png", sep=""), res=300, 2000, 1500)
ggplot(pcaproj_hubs_o, aes(x=treatment2, y=b_E2F_targets))+geom_boxplot()+theme_classic()
dev.off()



ggplot(pcaproj_hubs_o, aes(x=treatment2, y=TFDP1))+geom_boxplot()
ggplot(pcaproj_hubs_o, aes(x=E2F3, y=b_E2F_targets, colour=treatment2))+geom_point()
ggplot(pcaproj_hubs_o, aes(x=CEBPG, y=b_E2F_targets, colour=treatment2))+geom_point()
ggplot(pcaproj_hubs_o, aes(x=TFDP1, y=b_E2F_targets, colour=treatment2))+geom_point()
ggplot(pcaproj_hubs_o, aes(x=PTTG1, y=b_E2F_targets, colour=treatment2))+geom_point()
ggplot(pcaproj_hubs_o, aes(x=TEAD4, y=b_E2F_targets, colour=treatment2))+geom_point()


cd_o<-matrix(nrow=7, ncol=19)
ind<-0
for(l in unique(pcaproj_hubs_o$treatment2)){
  ind<-ind+1
  
  ctrl<-subset(pcaproj_hubs_o, treatment2=="CTRL")
  trt<-subset(pcaproj_hubs_o, treatment2==l)
  
  for(i in 1:19){
    cd_o[ind,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
  }
  
}

colnames(cd_o)<-colnames(pcaproj_hubs_o)[1:19]
rownames(cd_o)<-unique(pcaproj_hubs_o$treatment2)

df_o<-data.frame(cd=c(t(cd_o)), trt=rep(unique(pcaproj_hubs_o$treatment2), each= 19), corr=rep(cc["b_E2F_targets",colnames(cd_o)],7))

ggplot(df_o, aes(x=corr, y=cd))+geom_point()+facet_grid(~trt)

########## plot changes in MEs (Cohen's d) for each clone
## Clones are coherent (with some outliers)

toplot<-cd_o[-which(rownames(cd_o) %in% c("CTRL", "PSAT1")),]
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

graphics.off()
png(paste("results/",date, "/Cohen_bE2F_sh.png", sep=""),res=300, 1500, 2500)
pheatmap(t(toplot),   cellwidth=15, cellheight=15,  keep.dendro=T, color =myColor, breaks = myBreaks)
dev.off()

########## plot changes in MEs (Cohen's d) vs modules' correlation with b_E2F_targets
## for each KO, the most affected modules are either the most highly or lowly correlated with b_E2F_targets

df<-data.frame(cd=c(t(cd_o)), trt=rep(unique(pcaproj_hubs_o$treatment2), each= 19), corr=rep(cc["b_E2F_targets",colnames(cd_o)],7),
               module=rep(colnames(cd_o),7))

png(paste("results/",date, "/CohenVScorr_bE2F_sh.png", sep=""),res=300, 3000, 1500)
ggplot(subset(df, trt %in% c("CEBPG", "TEAD4", "E2F3", "PTTG1", "TFDP1")), aes(x=corr, y=cd, label=module))+geom_point(size=2)+facet_grid(~trt)+geom_smooth(method = lm)+stat_cor(label.x=-0.5, label.y = 7)+geom_text_repel(max.overlaps = 5)+
  theme_bw()+theme(strip.text=element_text(size = 12, face = "bold"))
dev.off()

################################
## for each sh separately


png(paste("results/",date, "/bE2F_project_sh_sep.png", sep=""), res=300, 2500, 1500)
ggplot(pcaproj_hubs_o, aes(x=treatment, y=b_E2F_targets))+geom_boxplot()+theme_classic()
dev.off()


cd_o_sep<-matrix(nrow=10, ncol=19)

ctrl<-subset(pcaproj_hubs_o, treatment=="SCR2")
trt<-subset(pcaproj_hubs_o, treatment=="C2")
for(i in 1:19){
  cd_o_sep[1,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}
ctrl<-subset(pcaproj_hubs_o, treatment=="SCR2")
trt<-subset(pcaproj_hubs_o, treatment=="C5")
for(i in 1:19){
  cd_o_sep[2,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}  
ctrl<-subset(pcaproj_hubs_o, treatment=="SCR2")
trt<-subset(pcaproj_hubs_o, treatment=="T1")
for(i in 1:19){
  cd_o_sep[3,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}
ctrl<-subset(pcaproj_hubs_o, treatment=="SCR2")
trt<-subset(pcaproj_hubs_o, treatment=="T4")
for(i in 1:19){
  cd_o_sep[4,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}
ctrl<-subset(pcaproj_hubs_o, treatment=="SCR4")
trt<-subset(pcaproj_hubs_o, treatment=="P1")
for(i in 1:19){
  cd_o_sep[5,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}
ctrl<-subset(pcaproj_hubs_o, treatment=="SCR4")
trt<-subset(pcaproj_hubs_o, treatment=="P4")
for(i in 1:19){
  cd_o_sep[6,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}
ctrl<-subset(pcaproj_hubs_o, treatment=="SCR4")
trt<-subset(pcaproj_hubs_o, treatment=="E5")
for(i in 1:19){
  cd_o_sep[7,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}
ctrl<-subset(pcaproj_hubs_o, treatment=="SCR5")
trt<-subset(pcaproj_hubs_o, treatment=="E6")
for(i in 1:19){
  cd_o_sep[8,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}
ctrl<-subset(pcaproj_hubs_o, treatment=="SCR3")
trt<-subset(pcaproj_hubs_o, treatment=="TF4")
for(i in 1:19){
  cd_o_sep[9,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}
ctrl<-subset(pcaproj_hubs_o, treatment=="SCR3")
trt<-subset(pcaproj_hubs_o, treatment=="TF5")
for(i in 1:19){
  cd_o_sep[10,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
}



colnames(cd_o_sep)<-colnames(pcaproj_hubs_o)[1:19]
rownames(cd_o_sep)<-c("CEBPG 2", "CEBPG 5", "TEAD4 1", "TEAD4 4", "PTTG1 1", "PTTG1 4", "E2F3 5", "E2F3 6", "TFDP1 4", "TFDP1 5")


########## plot changes in MEs (Cohen's d) for each sh

toplot<-cd_o_sep
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

graphics.off()
png(paste("results/",date, "/Cohen_bE2F_sh_sep.png", sep=""),res=300, 2000, 2500)
pheatmap(t(toplot),   cellwidth=15, cellheight=15,  keep.dendro=T, color =myColor, breaks = myBreaks)
dev.off()

########## plot changes in MEs (Cohen's d) vs modules' correlation with b_E2F_targets

df<-data.frame(cd=c(t(cd_o_sep)), trt=rep(rownames(cd_o_sep), each= 19), corr=rep(cc["b_E2F_targets",colnames(cd_o_sep)],10),
               module=rep(colnames(cd_o_sep),10))

png(paste("results/",date, "/CohenVScorr_bE2F_sh_sep.png", sep=""),res=300, 6000, 1500)
ggplot(df, aes(x=corr, y=cd, label=module))+geom_point(size=2)+facet_grid(~trt)+geom_smooth(method = lm)+stat_cor(label.x=-0.5, label.y = 25)+geom_text_repel(max.overlaps = 5)+
  scale_y_continuous(limits = c(-25, 25))+
  theme_bw()+theme(strip.text=element_text(size = 12, face = "bold"))
dev.off()


##########relationship with sh impact on the gene's expression
df<-cbind.data.frame(cd_o_sep, TF_log2FC=c(DE_CEBPG_2["CEBPG","log2FoldChange"], DE_CEBPG_5["CEBPG","log2FoldChange"],
                                           DE_TEAD4_1["TEAD4", "log2FoldChange"], DE_TEAD4_4["TEAD4", "log2FoldChange"],
                                           DE_PTTG1_1["PTTG1","log2FoldChange"], DE_PTTG1_4["PTTG1","log2FoldChange"], 
                                           DE_E2F3_5["E2F3","log2FoldChange"], DE_E2F3_6["E2F3","log2FoldChange"],
                                           DE_TFDP1_4["TFDP1","log2FoldChange"], DE_TFDP1_5["TFDP1","log2FoldChange"]),
                     TF=c("CEBPG 2", "CEBPG 5", "TEAD4 1", "TEAD4 4", "PTTG1 1", "PTTG1 4", "E2F3 5", "E2F3 6", "TFDP1 4", "TFDP1 5"))

png(paste("results/",date, "/CohenVSlog2FC_sh_sep.png", sep=""),res=300, 1500, 1500)
ggplot(df, aes(x=TF_log2FC, y=b_E2F_targets, label=TF))+geom_point()+geom_label_repel()+theme_bw()
dev.off()


##################################
############ GSEA
##################################

fgsea_MsigdbC2CP_sh<-list()
for(c in c("DE_CEBPG",
           "DE_TEAD4",
           "DE_PTTG1", "DE_E2F3", "DE_TFDP1", "DE_PSAT1")){
  
  i<-get(c)
  
  forgesea<-unlist(i$log2FoldChange)
  names(forgesea)<-rownames(i)
  forgesea<-forgesea[!is.na(forgesea)]
  fgsea_MsigdbC2CP_sh[[c]]<-GSEA(sort(forgesea, decreasing=T), TERM2GENE=m_df, pvalueCutoff = 1, maxGSSize = 10000)
}


allpaths<-c()
for(i in 1:length(fgsea_MsigdbC2CP_sh)){
  allpaths<-union(allpaths, fgsea_MsigdbC2CP_sh[[i]]$Description[fgsea_MsigdbC2CP_sh[[i]]$p.adjust<0.05])
}
allpaths_mat<-matrix(0,nrow=length(allpaths), ncol=length(fgsea_MsigdbC2CP_sh))
rownames(allpaths_mat)<-allpaths
for(i in 1:length(fgsea_MsigdbC2CP_sh)){
  allpaths_mat[fgsea_MsigdbC2CP_sh[[i]]$Description[which(fgsea_MsigdbC2CP_sh[[i]]$p.adjust<0.05)],i]<-fgsea_MsigdbC2CP_sh[[i]]$NES[which(fgsea_MsigdbC2CP_sh[[i]]$p.adjust<0.05)]
}


colnames(allpaths_mat)<-names(fgsea_MsigdbC2CP_sh)
allpaths_mat<-allpaths_mat[,-6]

shared_paths<-allpaths_mat[rowSums(allpaths_mat!=0)>0, ]
colnames(shared_paths)<-c("CEBPG",
                          "TEAD4",
                          "PTTG1", "E2F3", "TFDP1")

toplot<-shared_paths
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

graphics.off()
png(paste("results/",date, "/GSEA_MSigDB_Hallmarks_shared_sh.png", sep=""), res=300, 2000, 2500)
pheatmap(toplot,   cellwidth=15, cellheight=15,  keep.dendro=T, breaks = myBreaks, color = myColor)
dev.off()



fgsea_MsigdbC2CP_sh_sep<-list()
for(c in c("DE_CEBPG_2","DE_CEBPG_5",
           "DE_TEAD4_1", "DE_TEAD4_4",
           "DE_PTTG1_1","DE_PTTG1_4",
           "DE_E2F3_5","DE_E2F3_6", 
           "DE_TFDP1_4", "DE_TFDP1_5")){
  
  i<-get(c)
  
  forgesea<-unlist(i$log2FoldChange)
  names(forgesea)<-rownames(i)
  forgesea<-forgesea[!is.na(forgesea)]
  fgsea_MsigdbC2CP_sh_sep[[c]]<-GSEA(sort(forgesea, decreasing=T), TERM2GENE=m_df, pvalueCutoff = 1, maxGSSize = 10000)
}



allpaths<-c()
for(i in 1:length(fgsea_MsigdbC2CP_sh_sep)){
  allpaths<-union(allpaths, fgsea_MsigdbC2CP_sh_sep[[i]]$Description[fgsea_MsigdbC2CP_sh_sep[[i]]$p.adjust<0.05])
}
allpaths_mat<-matrix(0,nrow=length(allpaths), ncol=length(fgsea_MsigdbC2CP_sh_sep))
rownames(allpaths_mat)<-allpaths
for(i in 1:length(fgsea_MsigdbC2CP_sh_sep)){
  allpaths_mat[fgsea_MsigdbC2CP_sh_sep[[i]]$Description[which(fgsea_MsigdbC2CP_sh_sep[[i]]$p.adjust<0.05)],i]<-fgsea_MsigdbC2CP_sh_sep[[i]]$NES[which(fgsea_MsigdbC2CP_sh_sep[[i]]$p.adjust<0.05)]
}


colnames(allpaths_mat)<-names(fgsea_MsigdbC2CP_sh_sep)


shared_paths<-allpaths_mat[rowSums(allpaths_mat!=0)>0, ]
colnames(shared_paths)<-c("CEBPG_2","CEBPG_5",
                          "TEAD4_1", "TEAD4_4",
                          "PTTG1_1","PTTG1_4",
                          "E2F3_5","E2F3_6", 
                          "TFDP1_4", "TFDP1_5")

toplot<-shared_paths
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(unlist(toplot), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1),
              seq(max(unlist(toplot), na.rm=T)/paletteLength, max(unlist(toplot), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

graphics.off()
png(paste("results/",date, "/GSEA_MSigDB_Hallmarks_shared_sh_sep.png", sep=""), res=300, 3000, 4000)
pheatmap(toplot,   cellwidth=15, cellheight=15,  keep.dendro=T, breaks = myBreaks, color = myColor)
dev.off()



######################################
##### Gene Essentiality
#########################################

###load input data
load("data/Sanger_Broad_higQ_scaled_depFC.RData")
CMP_annot <- read.csv("data/model_list_20210611.csv") # from https://cog.sanger.ac.uk/cmp/download/model_list_20210611.csv

ess<-scaled_depFC[c("E2F3", "TFDP1", "TEAD4", "CEBPG", "PTTG1"),c("MDA-MB-468","MDA-MB-231","Hs-578-T")]
write.xlsx(data.frame(ess), file=paste("results/",date, "/essentiality.xlsx", sep=""), rowNames=T)



