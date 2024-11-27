library(DESeq2)
library(FactoMineR)
library(ggplot2)
library(openxlsx)
library(pheatmap)
library(WGCNA)
library(ggrepel)
library(ggpubr)

date<-"20241121"

######network data to load
load("data/Networks/centrality_basal.RData")
load("data/Networks/metabric.RData")
load("data/Networks/meta.RData")

########################################
##### FUNCTIONS #########################
#########################################

##single sample network function
source("code/ssMod.R")

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



########################################
##### ANALYSES #########################
#########################################


#######################################################
### normalization, filtering, log transformation
###################################################

counts<-read.xlsx("data/gene_count.xlsx")
rownames(counts)<-counts[,1]
anno<-counts[,c("gene_id", "gene_name")]

metadata<-read.xlsx("data/Novogene_metadata.xlsx", rowNames = T)
metadata$Clone<-factor(metadata$Clone)

counts_name<-changenames(counts[,c(2:34)], anno = anno)
RPM<-t(t(counts_name)/colSums(counts_name))*1000000

counts_name<-counts_name[rowSums(counts_name>=10)>2,]
RPM<-RPM[rownames(counts_name),]

RPMlog<-log2(RPM+1)

counts<-read.table("data/Poli_E2F3_KO-RNAseq-v1-run241011/RNAseq/dataset/v1-run241011/GEP.count", header = T, row.names = 1)
RPM<-t(t(counts)/colSums(counts))*1000000

counts<-counts[rowSums(counts>=10)>2,]
RPM<-RPM[rownames(counts),]
RPMlog2<-log2(RPM+1)

RPMtot<-cbind(RPMlog[intersect(rownames(RPMlog), rownames(RPMlog2)), ], RPMlog2[intersect(rownames(RPMlog), rownames(RPMlog2)), ])


metadata<-read.xlsx("data/Novogene_metadata.xlsx", rowNames = T)
metadata$Clone<-factor(metadata$Clone)
metadata$Seq<-"Novogene"

metadata2<-data.frame(Cell.line=rep("MDAMB468", 9), KO.gene=rep(c("EV", "E2F3", "E2F3"), each=3),
                      Clone=rep(c(100, 1, 3), each=3), Replicate=rep(1:3, 3), Seq=rep("Oliviero", 9))

metadata<-rbind.data.frame(metadata, metadata2)
metadata$Clone<-factor(metadata$Clone, levels=c(levels(metadata$Clone), 1, 3))

metadata$Clone[37:39]<-1
metadata$Clone[40:42]<-3

rownames(metadata)<-colnames(RPMtot)

###################################
ssmod<-ssMod(data_orig=metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], centrality=centrality_basal, data=RPMtot, modules="b_E2F_targets")

ssmod_E2F<-metadata
ssmod_E2F$ssMod<-colMeans(abs(ssmod[[1]]), na.rm=T)
ssmod_E2F$KO.gene<-factor(ssmod_E2F$KO.gene, levels=c("EV", "TFDP1", "E2F3"))
ssmod_E2F$treatment2<-paste(ssmod_E2F$KO.gene, ssmod_E2F$Clone)

####plot b_E2F_targets across KOs
ggplot(ssmod_E2F, aes(x=Cell.line, y=ssMod, fill=KO.gene))+geom_boxplot()+theme_classic()
ggplot(ssmod_E2F, aes(x=Cell.line, y=ssMod, fill=treatment2))+geom_boxplot()+theme_classic()


####relate ssmod changes and functional impact of the KO
load(file=paste("results/",date, "/dc.RData", sep=""))
load(file=paste("results/",date, "/d.RData", sep=""))
save(p, file=paste("results/",date, "/p.RData", sep=""))

#individual clones
tt<-array2DF(by(ssmod_E2F$ssMod, list(ssmod_E2F$Cell.line, ssmod_E2F$treatment2), mean))
tt[tt$Var1=="Hs578",]
tt[tt$Var1=="MDAMB231",]
tt[tt$Var1=="MDAMB468",]

ssmod_val<-tt$Value

ssmod_FC<-c(ssmod_val[c(13,16)]/ssmod_val[c(10)] , ssmod_val[c(20,26,2,5)]/ssmod_val[c(11)], ssmod_val[c(24,30,9,9)]/ssmod_val[c(12)])
plot(ssmod_FC, scale(dc[3,])+scale(dc[1,]))
cor.test(ssmod_FC, scale(dc[3,])+scale(dc[1,]), method="s")

##pooling clones
tt<-array2DF(by(ssmod_E2F$ssMod, list(ssmod_E2F$Cell.line, ssmod_E2F$KO.gene), mean))
tt[tt$Var1=="Hs578",]
tt[tt$Var1=="MDAMB231",]
tt[tt$Var1=="MDAMB468",]

ssmod_val<-tt$Value

ssmod_FC<-c(ssmod_val[c(4)]/ssmod_val[c(1)] , ssmod_val[c(5,8)]/ssmod_val[c(2)], ssmod_val[c(6,9)]/ssmod_val[c(3)])
plot(ssmod_FC, scale(d[3,])+scale(d[1,]))

cor.test(ssmod_FC, scale(d[3,])+scale(d[1,]), method="s")