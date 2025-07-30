library(WGCNA)
library(limma)
load("data/RData/meta.RData")
load("data/RData/metabric.RData")
load("results/2025/centrality_basal.RData")

fishertest_basal<-function(genes, dataset){
  fisherp<-matrix(ncol=length(unique(centrality_basal$module)), nrow=1)
 for(i in 1:length(unique(centrality_basal$module))){
  metabric_tmp<-metabric[rownames(metabric) %in% rownames(dataset),]
  moduleColors_basal_tmp<-centrality_basal[rownames(metabric) %in% rownames(dataset), "module"]
  counts<-matrix(c(length(intersect(rownames(metabric_tmp)[which(moduleColors_basal_tmp %in% unique(centrality_basal$module)[i])], genes )),
                   length(which(moduleColors_basal_tmp %in% unique(centrality_basal$module)[i]))-length(intersect(rownames(metabric_tmp)[which(moduleColors_basal_tmp %in% unique(centrality_basal$module)[i])], genes )),
                   length(intersect(genes, rownames(metabric_tmp)[-which(moduleColors_basal_tmp %in% c("grey", unique(centrality_basal$module)[i]))])),
                   length(setdiff(rownames(metabric_tmp), union(rownames(metabric_tmp)[which(moduleColors_basal_tmp %in% c("grey", unique(centrality_basal$module)[i]))], genes)))), nrow=2)
  fisherp[1,i]<-fisher.test(counts, alternative = "greater")[[1]]
 }

colnames(fisherp)<-unique(centrality_basal$module)
return(fisherp)
}


#change names in gene symbols
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

#returns DE genes, positive logFC when the genes are higher in condition 2
DE<-function(cond1, cond2, data){
  gs <- factor(c(rep("Cond1",ncol(cond1)),rep("Cond2",ncol(cond2))))
  design <- model.matrix(~gs + 0, cbind.data.frame(cond1, cond2))
  fit <- lmFit(cbind.data.frame(cond1, cond2), design)  # fit linear model
  # set up contrasts of interest and recalculate model coefficients
  #cts <- paste(levels(gs)[1], levels(gs)[2], sep="-")
  cts<-"gsCond2-gsCond1"
  cont.matrix <- makeContrasts(contrasts=cts, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  
  # compute statistics and table of top significant genes
  fit2 <- eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(data))
  return(tT)
}



#####FOXM1
FOXM1_GSE2222=read.csv("data/datasets/GSE2222_series_matrix.txt", sep="\t", row.names = 1)
library(hgu133a.db)
x <- hgu133aSYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
adf <- data.frame(probe=names(vals), gene=vals)

FOXM1_GSE2222<-changenames(data=FOXM1_GSE2222, anno=cbind(rownames(FOXM1_GSE2222), adf[match( rownames(FOXM1_GSE2222), adf[,1]),2]))

#cond1 is the ctrl, cond2 is the trt (FOXM1 siRNA)
DE_FOXM1_GSE2222<-DE(cond1=FOXM1_GSE2222[,4:6], cond2=FOXM1_GSE2222[,7:9], FOXM1_GSE2222)

genes_up_FOXM1_GSE2222<-rownames(DE_FOXM1_GSE2222)[which(DE_FOXM1_GSE2222$P.Value<0.05 & DE_FOXM1_GSE2222$logFC>0)]
genes_down_FOXM1_GSE2222<-rownames(DE_FOXM1_GSE2222)[which(DE_FOXM1_GSE2222$P.Value<0.05 & DE_FOXM1_GSE2222$logFC<0)]

ETup_1<-fishertest_basal(genes_up_FOXM1_GSE2222, FOXM1_GSE2222)
ETdn_1<-fishertest_basal(genes_down_FOXM1_GSE2222, FOXM1_GSE2222)
ETall_1<-fishertest_basal(c(genes_down_FOXM1_GSE2222,genes_up_FOXM1_GSE2222), FOXM1_GSE2222)

#####FOXM1
FOXM1_GSE55204=read.csv("data/datasets/GSE55204_series_matrix.txt", sep="\t", row.names = 1)
library(hgu133a2.db)
x <- hgu133a2SYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
adf <- data.frame(probe=names(vals), gene=vals)

FOXM1_GSE55204<-changenames(data=FOXM1_GSE55204, anno=cbind(rownames(FOXM1_GSE55204), adf[match( rownames(FOXM1_GSE55204), adf[,1]),2]))

#cond1 is the ctrl, cond2 is the trt (FOXM1 siRNA)
DE_FOXM1_GSE55204<-DE(cond1=FOXM1_GSE55204[,1:2], cond2=FOXM1_GSE55204[,5:6], FOXM1_GSE55204)

genes_up_FOXM1_GSE55204<-rownames(DE_FOXM1_GSE55204)[which(DE_FOXM1_GSE55204$P.Value<0.05 & DE_FOXM1_GSE55204$logFC>0)]
genes_down_FOXM1_GSE55204<-rownames(DE_FOXM1_GSE55204)[which(DE_FOXM1_GSE55204$P.Value<0.05 & DE_FOXM1_GSE55204$logFC<0)]

ETup_2<-fishertest_basal(genes_up_FOXM1_GSE55204, FOXM1_GSE55204)
ETdn_2<-fishertest_basal(genes_down_FOXM1_GSE55204, FOXM1_GSE55204)
ETall_2<-fishertest_basal(c(genes_down_FOXM1_GSE55204,genes_up_FOXM1_GSE55204 ), FOXM1_GSE55204)

###############
FOXM1_GSE25741=read.csv("data/datasets/GSE25741_series_matrix.txt", sep="\t", row.names = 1)
library("illuminaHumanv3.db")
x <- illuminaHumanv3SYMBOL
# Get the probe identifiers that are mapped to a gene name
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
adf <- data.frame(probe=names(vals), gene=vals)

FOXM1_GSE25741<-changenames(data=FOXM1_GSE25741, anno=cbind(rownames(FOXM1_GSE25741), adf[match( rownames(FOXM1_GSE25741), adf[,1]),2]))

#cond1 is the ctrl, cond2 is the trt (FOXM1 siRNA)
DE_FOXM1_GSE25741<-DE(cond1=FOXM1_GSE25741[,1:3], cond2=FOXM1_GSE25741[,4:6], FOXM1_GSE25741)

genes_up_FOXM1_GSE25741<-rownames(DE_FOXM1_GSE25741)[which(DE_FOXM1_GSE25741$P.Value<0.05 & DE_FOXM1_GSE25741$logFC>0)]
genes_down_FOXM1_GSE25741<-rownames(DE_FOXM1_GSE25741)[which(DE_FOXM1_GSE25741$P.Value<0.05 & DE_FOXM1_GSE25741$logFC<0)]

ETup_3<-fishertest_basal(genes_up_FOXM1_GSE25741, FOXM1_GSE25741)
ETdn_3<-fishertest_basal(genes_down_FOXM1_GSE25741, FOXM1_GSE25741)
ETall_3<-fishertest_basal(c(genes_down_FOXM1_GSE25741, genes_up_FOXM1_GSE25741), FOXM1_GSE25741)

#####PTTG1
PTTG1=read.csv("data/datasets/GSE48928_series_matrix.txt", sep="\t", row.names = 1)
library(illuminaHumanv4.db)
x <- illuminaHumanv4SYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
adf <- data.frame(probe=names(vals), gene=vals)

PTTG1<-changenames(data=PTTG1, anno=cbind(rownames(PTTG1), adf[match( rownames(PTTG1), adf[,1]),2]))

#cond1 is the ctrl, cond2 is the trt (PTTG1 siRNA)
DE_PTTG1<-DE(cond1=PTTG1[,1:6], cond2=PTTG1[,7:12], PTTG1)

genes_up_PTTG1<-rownames(DE_PTTG1)[which(DE_PTTG1$P.Value<0.05 & DE_PTTG1$logFC>0)]
genes_down_PTTG1<-rownames(DE_PTTG1)[which(DE_PTTG1$P.Value<0.05 & DE_PTTG1$logFC<0)]

ETup_4<-fishertest_basal(genes_up_PTTG1, PTTG1)
ETdn_4<-fishertest_basal(genes_down_PTTG1, PTTG1)
ETall_4<-fishertest_basal(c(genes_down_PTTG1,genes_up_PTTG1), PTTG1)


#####EZH2 GSE103242
EZH2=read.csv("data/datasets/MCF7_overEZH2_readCounts.txt", sep="\t", row.names = 1)

EZH2_RPKM<-(t(t(EZH2)/colSums(EZH2))*10^6)/EZH2[,1]
EZH2_RPKM<-EZH2_RPKM[,-1]

#cond1 is the ctrl, cond2 is the trt (EZH2 OE)
DE_EZH2<-DE(cond1=EZH2_RPKM[,5:6], cond2=EZH2_RPKM[,11:12], EZH2_RPKM)

genes_up_EZH2<-rownames(DE_EZH2)[which(DE_EZH2$P.Value<0.05 & DE_EZH2$logFC>0)]
genes_down_EZH2<-rownames(DE_EZH2)[which(DE_EZH2$P.Value<0.05 & DE_EZH2$logFC<0)]

ETup_5<-fishertest_basal(genes_up_EZH2, EZH2_RPKM)
ETdn_5<-fishertest_basal(genes_down_EZH2, EZH2_RPKM)
ETall_5<-fishertest_basal(c(genes_down_EZH2,genes_up_EZH2), EZH2_RPKM)


#####EZH2
EZH2_GSE48979=read.csv("data/datasets/GSE48979_series_matrix.txt", sep="\t", row.names = 1)
library(illuminaHumanv3.db)
x <- illuminaHumanv3SYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
adf <- data.frame(probe=names(vals), gene=vals)

EZH2_GSE48979<-changenames(data=EZH2_GSE48979, anno=cbind(rownames(EZH2_GSE48979), adf[match( rownames(EZH2_GSE48979), adf[,1]),2]))

#cond1 is the ctrl, cond2 is the trt (EZH2 siRNA)
DE_EZH2_GSE48979<-DE(cond1=EZH2_GSE48979[,c(1:3, 7:9)], cond2=EZH2_GSE48979[,10:12], EZH2_GSE48979)

genes_up_EZH2_GSE48979<-rownames(DE_EZH2_GSE48979)[which(DE_EZH2_GSE48979$P.Value<0.05 & DE_EZH2_GSE48979$logFC>0)]
genes_down_EZH2_GSE48979<-rownames(DE_EZH2_GSE48979)[which(DE_EZH2_GSE48979$P.Value<0.05 & DE_EZH2_GSE48979$logFC<0)]

ETup_6<-fishertest_basal(genes_up_EZH2_GSE48979, EZH2_GSE48979)
ETdn_6<-fishertest_basal(genes_down_EZH2_GSE48979, EZH2_GSE48979)
ETall_6<-fishertest_basal(c(genes_down_EZH2_GSE48979, genes_up_EZH2_GSE48979), EZH2_GSE48979)



#####FOXC1
FOXC1_GSE73234=read.csv("data/datasets/GSE73234_series_matrix.txt", sep="\t", row.names = 1, header=F)
FOXC1_GSE73234<-FOXC1_GSE73234[-1,]
colnames(FOXC1_GSE73234)<-FOXC1_GSE73234[1,]
FOXC1_GSE73234<-FOXC1_GSE73234[-1,]

FOXC1mat<-as.data.frame(lapply(FOXC1_GSE73234,function(x) as.numeric(as.character(x))))
colnames(FOXC1mat)<-colnames(FOXC1_GSE73234)
rownames(FOXC1mat)<-rownames(FOXC1_GSE73234)
FOXC1_GSE73234<-FOXC1mat

Affy<-read.csv("data/datasets/AffyHG1st.txt")
adf<-data.frame(probe=as.character(Affy$AFFY.HuGene.1.0.st.v1.probe),
                gene=as.character(Affy$HGNC.symbol))

FOXC1_GSE73234<-changenames(data=FOXC1_GSE73234, anno=cbind(rownames(FOXC1_GSE73234), adf[match( rownames(FOXC1_GSE73234), adf[,1]),2]))

#cond1 is the ctrl, cond2 is the trt (FOXC1 OE)
DE_FOXC1_GSE73234<-DE(cond1=FOXC1_GSE73234[,1:2], cond2=FOXC1_GSE73234[,3:4], FOXC1_GSE73234)

genes_up_FOXC1_GSE73234<-rownames(DE_FOXC1_GSE73234)[which(DE_FOXC1_GSE73234$P.Value<0.05 & DE_FOXC1_GSE73234$logFC>0)]
genes_down_FOXC1_GSE73234<-rownames(DE_FOXC1_GSE73234)[which(DE_FOXC1_GSE73234$P.Value<0.05 & DE_FOXC1_GSE73234$logFC<0)]

ETup_7<-fishertest_basal(genes_up_FOXC1_GSE73234, FOXC1_GSE73234)
ETdn_7<-fishertest_basal(genes_down_FOXC1_GSE73234, FOXC1_GSE73234)
ETall_7<-fishertest_basal(c(genes_down_FOXC1_GSE73234,genes_up_FOXC1_GSE73234), FOXC1_GSE73234)

############

#####FOXC1 GSE31912 removed as no replicates...
##########?


#################################

#######################

###########################
EZH2_GSE36939=read.csv("data/datasets/GSE36939_series_matrix.txt", sep="\t", row.names = 1)

Affy<-read.csv("data/datasets/AffyHG1st.txt")
adf <- data.frame(probe=as.character(Affy$AFFY.HuGene.1.0.st.v1.probe),
                  gene=as.character(Affy$HGNC.symbol))

EZH2_GSE36939<-changenames(data=EZH2_GSE36939, anno=cbind(rownames(EZH2_GSE36939), adf[match( rownames(EZH2_GSE36939), adf[,1]),2]))

#cond1 is the ctrl, cond2 is the trt (shEZH2), HCC70
DE_EZH2_GSE36939_HCC70<-DE(cond1=EZH2_GSE36939[,1:2], cond2=EZH2_GSE36939[,3:4], EZH2_GSE36939)

genes_up_EZH2_GSE36939_HCC70<-rownames(DE_EZH2_GSE36939_HCC70)[which(DE_EZH2_GSE36939_HCC70$P.Value<0.05 & DE_EZH2_GSE36939_HCC70$logFC>0)]
genes_down_EZH2_GSE36939_HCC70<-rownames(DE_EZH2_GSE36939_HCC70)[which(DE_EZH2_GSE36939_HCC70$P.Value<0.05 & DE_EZH2_GSE36939_HCC70$logFC<0)]

ETup_8<-fishertest_basal(genes_up_EZH2_GSE36939_HCC70, EZH2_GSE36939)
ETdn_8<-fishertest_basal(genes_down_EZH2_GSE36939_HCC70, EZH2_GSE36939)
ETall_8<-fishertest_basal(c(genes_down_EZH2_GSE36939_HCC70,genes_up_EZH2_GSE36939_HCC70), EZH2_GSE36939)

#cond1 is the ctrl, cond2 is the trt (shEZH2), MD468
DE_EZH2_GSE36939_MD468<-DE(cond1=EZH2_GSE36939[,5:6], cond2=EZH2_GSE36939[,7:8], EZH2_GSE36939)

genes_up_EZH2_GSE36939_MD468<-rownames(DE_EZH2_GSE36939_MD468)[which(DE_EZH2_GSE36939_MD468$P.Value<0.05 & DE_EZH2_GSE36939_MD468$logFC>0)]
genes_down_EZH2_GSE36939_MD468<-rownames(DE_EZH2_GSE36939_MD468)[which(DE_EZH2_GSE36939_MD468$P.Value<0.05 & DE_EZH2_GSE36939_MD468$logFC<0)]

ETup_9<-fishertest_basal(genes_up_EZH2_GSE36939_MD468, EZH2_GSE36939)
ETdn_9<-fishertest_basal(genes_down_EZH2_GSE36939_MD468, EZH2_GSE36939)
ETall_9<-fishertest_basal(c(genes_down_EZH2_GSE36939_MD468,genes_up_EZH2_GSE36939_MD468), EZH2_GSE36939)


#########HMGA1
huex <- read.csv("data/datasets/GPL5175-3188.txt", sep="\t", stringsAsFactors = F, header=T)

mylist<-strsplit(huex$gene_assignment, " ")
genes<-sapply(mylist,function(x) x[3])

###Annotation
HMGA1=read.csv("data/datasets/GSE45483_series_matrix.txt", sep="\t", row.names = 1, header=T)

HMGA1<-changenames(data=HMGA1, anno=cbind(rownames(HMGA1), genes[match( rownames(HMGA1), huex[,1])]))

#cond1 is the ctrl, cond2 is the trt (shHMGA1)
DE_HMGA1<-DE(cond1=HMGA1[,1:3], cond2=HMGA1[,4:6], HMGA1)

genes_up_HMGA1<-rownames(DE_HMGA1)[which(DE_HMGA1$P.Value<0.05 & DE_HMGA1$logFC>0)]
genes_down_HMGA1<-rownames(DE_HMGA1)[which(DE_HMGA1$P.Value<0.05 & DE_HMGA1$logFC<0)]

ETup_10<-fishertest_basal(genes_up_HMGA1, HMGA1)
ETdn_10<-fishertest_basal(genes_down_HMGA1, HMGA1)
ETall_10<-fishertest_basal(c(genes_down_HMGA1,genes_up_HMGA1), HMGA1)

#########HMGA1
HMGA1_GSE35525=read.csv("data/datasets/GSE35525_series_matrix.txt", sep="\t", row.names = 1, header=T)
library(hgu133plus2.db)

x <- hgu133plus2SYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
adf <- data.frame(probe=names(vals), gene=vals)


HMGA1_GSE35525<-changenames(data=HMGA1_GSE35525, anno=cbind(rownames(HMGA1_GSE35525), adf[match( rownames(HMGA1_GSE35525), adf[,1]),2]))

#cond1 is the ctrl, cond2 is the trt (shHMGA1)
DE_HMGA1_GSE35525<-DE(cond1=HMGA1_GSE35525[,1:4], cond2=HMGA1_GSE35525[,5:8], HMGA1_GSE35525)

genes_up_HMGA1_GSE35525<-rownames(DE_HMGA1_GSE35525)[which(DE_HMGA1_GSE35525$P.Value<0.05 & DE_HMGA1_GSE35525$logFC>0)]
genes_down_HMGA1_GSE35525<-rownames(DE_HMGA1_GSE35525)[which(DE_HMGA1_GSE35525$P.Value<0.05 & DE_HMGA1_GSE35525$logFC<0)]

ETup_11<-fishertest_basal(genes_up_HMGA1_GSE35525, HMGA1_GSE35525)
ETdn_11<-fishertest_basal(genes_down_HMGA1_GSE35525, HMGA1_GSE35525)
ETall_11<-fishertest_basal(c(genes_down_HMGA1_GSE35525,genes_up_HMGA1_GSE35525), HMGA1_GSE35525)


#########TCF7L1_GSE38893 removed as only one replicate of the control

#####SSRP1
SSRP1=read.csv("data/datasets/GSE92281_series_matrix.txt", sep="\t", row.names = 1, header=F, stringsAsFactors = F)
SSRP1mat<-as.data.frame(lapply(SSRP1,function(x) as.numeric((x))))
colnames(SSRP1mat)<-colnames(SSRP1)
rownames(SSRP1mat)<-rownames(SSRP1)
SSRP1<-SSRP1mat

GPL<-read.csv("data/datasets/GPL10558-50081.txt", sep="\t", header=T)
adf<-data.frame(probe=as.character(GPL$Probe_Id),
                gene=as.character(GPL$Symbol))

SSRP1<-changenames(data=SSRP1, anno=cbind(rownames(SSRP1), adf[match( rownames(SSRP1), adf[,1]),2]))

#cond1 is the ctrl, cond2 is the trt (sh SSRP1)
DE_SSRP1_T47D<-DE(cond1=SSRP1[,5:6], cond2=SSRP1[,7:8], SSRP1)

genes_up_SSRP1_T47D<-rownames(DE_SSRP1_T47D)[which(DE_SSRP1_T47D$P.Value<0.05 & DE_SSRP1_T47D$logFC>0)]
genes_down_SSRP1_T47D<-rownames(DE_SSRP1_T47D)[which(DE_SSRP1_T47D$P.Value<0.05 & DE_SSRP1_T47D$logFC<0)]

ETup_12<-fishertest_basal(genes_up_SSRP1_T47D, SSRP1)
ETdn_12<-fishertest_basal(genes_down_SSRP1_T47D, SSRP1)
ETall_12<-fishertest_basal(c(genes_down_SSRP1_T47D, genes_up_SSRP1_T47D), SSRP1)

DE_SSRP1_MCF7<-DE(cond1=SSRP1[,9:10], cond2=SSRP1[,11:12], SSRP1)

genes_up_SSRP1_MCF7<-rownames(DE_SSRP1_MCF7)[which(DE_SSRP1_MCF7$P.Value<0.05 & DE_SSRP1_MCF7$logFC>0)]
genes_down_SSRP1_MCF7<-rownames(DE_SSRP1_MCF7)[which(DE_SSRP1_MCF7$P.Value<0.05 & DE_SSRP1_MCF7$logFC<0)]

ETup_13<-fishertest_basal(genes_up_SSRP1_MCF7, SSRP1)
ETdn_13<-fishertest_basal(genes_down_SSRP1_MCF7, SSRP1)
ETall_13<-fishertest_basal(c(genes_down_SSRP1_MCF7,genes_up_SSRP1_MCF7), SSRP1)

############
#########ELF5
ELF5_GSE30405=read.csv("data/datasets/GSE30405_series_matrix.txt", sep="\t", row.names = 1, header=T)
Affy<-read.csv("data/datasets/AffyHG1st.txt")
adf<-data.frame(probe=as.character(Affy$AFFY.HuGene.1.0.st.v1.probe),
                gene=as.character(Affy$HGNC.symbol))

ELF5_GSE30405<-changenames(data=ELF5_GSE30405, anno=cbind(rownames(ELF5_GSE30405), adf[match( rownames(ELF5_GSE30405), adf[,1]),2]))

#cond1 is the ctrl, cond2 is the trt (shELF5)
DE_ELF5_GSE30405_HCC<-DE(cond1=ELF5_GSE30405[,1:3], cond2=ELF5_GSE30405[,4:6], SSRP1)

genes_up_ELF5_GSE30405_HCC<-rownames(DE_ELF5_GSE30405_HCC)[which(DE_ELF5_GSE30405_HCC$P.Value<0.05 & DE_ELF5_GSE30405_HCC$logFC>0)]
genes_down_ELF5_GSE30405_HCC<-rownames(DE_ELF5_GSE30405_HCC)[which(DE_ELF5_GSE30405_HCC$P.Value<0.05 & DE_ELF5_GSE30405_HCC$logFC<0)]

ETup_14<-fishertest_basal(genes_up_ELF5_GSE30405_HCC, ELF5_GSE30405)
ETdn_14<-fishertest_basal(genes_down_ELF5_GSE30405_HCC, ELF5_GSE30405)
ETall_14<-fishertest_basal(c(genes_down_ELF5_GSE30405_HCC,genes_up_ELF5_GSE30405_HCC), ELF5_GSE30405)


DE_ELF5_GSE30405_MCF7<-DE(cond1=ELF5_GSE30405[,7:8], cond2=ELF5_GSE30405[,9:10], SSRP1)

genes_up_ELF5_GSE30405_MCF7<-rownames(DE_ELF5_GSE30405_MCF7)[which(DE_ELF5_GSE30405_MCF7$P.Value<0.05 & DE_ELF5_GSE30405_MCF7$logFC>0)]
genes_down_ELF5_GSE30405_MCF7<-rownames(DE_ELF5_GSE30405_MCF7)[which(DE_ELF5_GSE30405_MCF7$P.Value<0.05 & DE_ELF5_GSE30405_MCF7$logFC<0)]

ETup_15<-fishertest_basal(genes_up_ELF5_GSE30405_MCF7, ELF5_GSE30405)
ETdn_15<-fishertest_basal(genes_down_ELF5_GSE30405_MCF7, ELF5_GSE30405)
ETall_15<-fishertest_basal(c(genes_down_ELF5_GSE30405_MCF7,genes_up_ELF5_GSE30405_MCF7) , ELF5_GSE30405)

DE_ELF5_GSE30405_T47D<-DE(cond1=ELF5_GSE30405[,11:12], cond2=ELF5_GSE30405[,13:14], SSRP1)

genes_up_ELF5_GSE30405_T47D<-rownames(DE_ELF5_GSE30405_T47D)[which(DE_ELF5_GSE30405_T47D$P.Value<0.05 & DE_ELF5_GSE30405_T47D$logFC>0)]
genes_down_ELF5_GSE30405_T47D<-rownames(DE_ELF5_GSE30405_T47D)[which(DE_ELF5_GSE30405_T47D$P.Value<0.05 & DE_ELF5_GSE30405_T47D$logFC<0)]

ETup_16<-fishertest_basal(genes_up_ELF5_GSE30405_T47D, ELF5_GSE30405)
ETdn_16<-fishertest_basal(genes_down_ELF5_GSE30405_T47D, ELF5_GSE30405)
ETall_16<-fishertest_basal(c(genes_down_ELF5_GSE30405_MCF7,genes_up_ELF5_GSE30405_T47D), ELF5_GSE30405)


#####################
#########ZNF165
ZNF165_GSE63984=read.csv("data/datasets/GSE63984_series_matrix.txt", sep="\t", row.names = 1, header=T)

huex <- read.csv("data/datasets/GPL11532-32230.txt", sep="\t", stringsAsFactors = F, header=T)

mylist<-strsplit(huex$gene_assignment, " ")
genes<-sapply(mylist,function(x) x[3])

###Annotation

ZNF165_GSE63984<-changenames(data=ZNF165_GSE63984, anno=cbind(rownames(ZNF165_GSE63984), genes[match( rownames(ZNF165_GSE63984), huex[,1])]))

#cond1 is the ctrl, cond2 is the trt (siZNF165), SUM159
DE_ZNF165_GSE63984<-DE(cond1=ZNF165_GSE63984[,c(1,5,9)], cond2=ZNF165_GSE63984[,c(3,7,11)], ZNF165_GSE63984)

genes_up_ZNF165_GSE63984<-rownames(DE_ZNF165_GSE63984)[which(DE_ZNF165_GSE63984$P.Value<0.05 & DE_ZNF165_GSE63984$logFC>0)]
genes_down_ZNF165_GSE63984<-rownames(DE_ZNF165_GSE63984)[which(DE_ZNF165_GSE63984$P.Value<0.05 & DE_ZNF165_GSE63984$logFC<0)]

ETup_17<-fishertest_basal(genes_up_ZNF165_GSE63984, ZNF165_GSE63984)
ETdn_17<-fishertest_basal(genes_down_ZNF165_GSE63984, ZNF165_GSE63984)
ETall_17<-fishertest_basal(c(genes_down_ZNF165_GSE63984,genes_up_ZNF165_GSE63984), ZNF165_GSE63984)

#cond1 is the ctrl, cond2 is the trt (siZNF165), WHIM12
DE_ZNF165_GSE63984_WHIM12<-DE(cond1=ZNF165_GSE63984[,c(1,5,9)+12], cond2=ZNF165_GSE63984[,c(3,7,11)+12], ZNF165_GSE63984)

genes_up_ZNF165_GSE63984_WHIM12<-rownames(DE_ZNF165_GSE63984_WHIM12)[which(DE_ZNF165_GSE63984_WHIM12$P.Value<0.05 & DE_ZNF165_GSE63984_WHIM12$logFC>0)]
genes_down_ZNF165_GSE63984_WHIM12<-rownames(DE_ZNF165_GSE63984_WHIM12)[which(DE_ZNF165_GSE63984_WHIM12$P.Value<0.05 & DE_ZNF165_GSE63984_WHIM12$logFC<0)]

ETup_18<-fishertest_basal(genes_up_ZNF165_GSE63984_WHIM12, ZNF165_GSE63984)
ETdn_18<-fishertest_basal(genes_down_ZNF165_GSE63984_WHIM12, ZNF165_GSE63984)
ETall_18<-fishertest_basal(c(genes_down_ZNF165_GSE63984_WHIM12,genes_up_ZNF165_GSE63984_WHIM12), ZNF165_GSE63984)


#########PARP1
PARP1_GSE34817=read.csv("data/datasets/GSE34817_series_matrix.txt", sep="\t", row.names = 1, header=T, stringsAsFactors = F)
library(illuminaHumanv4.db)

x <- illuminaHumanv4SYMBOL
# Get the probe identifiers that are mapped to a gene name
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
vals <- sapply(xx, as.vector)
adf <- data.frame(probe=names(vals), gene=vals)

PARP1_GSE34817<-changenames(data=PARP1_GSE34817, anno=cbind(rownames(PARP1_GSE34817), adf[match( rownames(PARP1_GSE34817), adf[,1]),2]))

#cond1 is the ctrl, cond2 is the trt (PARP1 inhibitor AG014699 10nM)
DE_PARP1_GSE34817_10nM<-DE(cond1=PARP1_GSE34817[,10:12], cond2=PARP1_GSE34817[,7:9], PARP1_GSE34817)

genes_up_PARP1_GSE34817_10nM<-rownames(DE_PARP1_GSE34817_10nM)[which(DE_PARP1_GSE34817_10nM$P.Value<0.05 & DE_PARP1_GSE34817_10nM$logFC>0)]
genes_down_PARP1_GSE34817_10nM<-rownames(DE_PARP1_GSE34817_10nM)[which(DE_PARP1_GSE34817_10nM$P.Value<0.05 & DE_PARP1_GSE34817_10nM$logFC<0)]

ETup_19<-fishertest_basal(genes_up_PARP1_GSE34817_10nM, PARP1_GSE34817)
ETdn_19<-fishertest_basal(genes_down_PARP1_GSE34817_10nM, PARP1_GSE34817)
ETall_19<-fishertest_basal(c(genes_down_PARP1_GSE34817_10nM, genes_up_PARP1_GSE34817_10nM), PARP1_GSE34817)

#cond1 is the ctrl, cond2 is the trt (PARP1 inhibitor AG014699 100nM)
DE_PARP1_GSE34817_100nM<-DE(cond1=PARP1_GSE34817[,10:12], cond2=PARP1_GSE34817[,1:3], PARP1_GSE34817)

genes_up_PARP1_GSE34817_100nM<-rownames(DE_PARP1_GSE34817_100nM)[which(DE_PARP1_GSE34817_100nM$P.Value<0.05 & DE_PARP1_GSE34817_100nM$logFC>0)]
genes_down_PARP1_GSE34817_100nM<-rownames(DE_PARP1_GSE34817_100nM)[which(DE_PARP1_GSE34817_100nM$P.Value<0.05 & DE_PARP1_GSE34817_100nM$logFC<0)]

ETup_20<-fishertest_basal(genes_up_PARP1_GSE34817_100nM, PARP1_GSE34817)
ETdn_20<-fishertest_basal(genes_down_PARP1_GSE34817_100nM, PARP1_GSE34817)
ETall_20<-fishertest_basal(c(genes_down_PARP1_GSE34817_100nM,genes_up_PARP1_GSE34817_100nM), PARP1_GSE34817)


####PLOT
#removed long term treatments (day10)
DE_lists<-c("FOXM1_GSE2222", "FOXM1_GSE55204", "FOXM1_GSE25741", "PTTG1_GSE48928", "EZH2_OE_GSE103242",
            "EZH2_GSE48979", "FOXC1_OE_GSE73234", "EZH2_GSE36939_HCC70", "EZH2_GSE36939_MD468", 
            "HMGA1_GSE45483", "HMGA1_GSE35525", "SSRP1_GSE92281_T47D", "SSRP1_GSE92281_MCF7", 
            
            "ELF5_GSE30405_HCC", "ELF5_GSE30405_MCF7" , "ELF5_GSE30405_T47D", "ZNF165_GSE63984", "ZNF165_GSE63984_WHIM12", "PARP1_GSE34817_10nM", "PARP1_GSE34817_100nM"
)


mat<-rbind(ETdn_1, ETdn_2, ETdn_3, ETdn_4, ETup_5, ETdn_6, ETup_7,
           ETdn_8, ETdn_9, ETdn_10, ETdn_11, ETdn_12, ETdn_13, ETdn_14, ETdn_15, ETdn_16, ETdn_17, ETdn_18, ETdn_19, ETdn_20 )

rownames(mat)<-DE_lists

mat<-mat[,-which(colnames(mat)=="b_Unconnected")]

mat_adj<-mat
for(col in 1:ncol(mat_adj)){
  mat_adj[,col]<-p.adjust(mat[,col], method="BH")
}
mat_adj[mat_adj<2.2*10^(-16)]<-2.2*10^(-16)

library(metap)
istwo <- rep(F, nrow(mat_adj))
toinvert <- rep(F, nrow(mat_adj))
pmerged<-c()
for(col in 1:ncol(mat_adj)){
  pmerged<-c(pmerged, sumlog(two2one(mat_adj[,col], two = istwo, invert = toinvert))$p)
}
names(pmerged)<-colnames(mat_adj)


library(pheatmap)

colnames(mat_adj)<-colnames(mat)
rownames(mat_adj)<-rownames(mat)
anno_p<-data.frame(mergedp= -log10(pmerged), npositive=colSums(mat_adj<0.05))
rownames(anno_p)<-colnames(mat_adj)

graphics.off()
png("results/2025/hubs_datasets_summary_down.png",res=300, 4000,3000)
pheatmap(-log10(t(mat_adj)), cellwidth=15, cellheight=15, keep.dendro=T, annotation_row = anno_p)
dev.off()

##########################################
############################################


mat<-rbind(ETup_1, ETup_2, ETup_3, ETup_4, ETdn_5, ETup_6, ETdn_7,
           ETup_8, ETdn_9, ETup_10, ETup_11, ETup_12, ETup_13, ETup_14, ETup_15, ETup_16,
           ETup_17, ETup_18, ETup_19, ETup_20 )

rownames(mat)<-DE_lists

mat<-mat[,-which(colnames(mat)=="b_Unconnected")]

mat_adj<-mat
for(col in 1:ncol(mat_adj)){
  mat_adj[,col]<-p.adjust(mat[,col], method="BH")
}
mat_adj[mat_adj<2.2*10^(-16)]<-2.2*10^(-16)

library(metap)
istwo <- rep(F, nrow(mat_adj))
toinvert <- rep(F, nrow(mat_adj))
pmerged<-c()
for(col in 1:ncol(mat_adj)){
  pmerged<-c(pmerged, sumlog(two2one(mat_adj[,col], two = istwo, invert = toinvert))$p)
}



library(pheatmap)

colnames(mat_adj)<-colnames(mat)
rownames(mat_adj)<-rownames(mat)
anno_p<-data.frame(mergedp= -log10(pmerged), npositive=colSums(mat_adj<0.05))
rownames(anno_p)<-colnames(mat_adj)

graphics.off()
png("results/2025/hubs_datasets_summary_up.png",res=300, 4000,3000)
pheatmap(-log10(t(mat_adj)), cellwidth=15, cellheight=15, keep.dendro=T, annotation_row = anno_p)
dev.off()



###############################

mat<-rbind(ETall_1, ETall_2, ETall_3, ETall_4, ETall_5, ETall_6, ETall_7,
           ETall_8, ETall_9, ETall_10, ETall_11, ETall_12, ETall_13, ETall_14, ETall_15, ETall_16,
           ETall_17, ETall_18, ETall_19, ETall_20 )



rownames(mat)<-DE_lists

mat<-mat[,-which(colnames(mat)=="b_Unconnected")]

mat_adj<-mat
for(col in 1:ncol(mat_adj)){
  mat_adj[,col]<-p.adjust(mat[,col], method="BH")
}
mat_adj[mat_adj<2.2*10^(-16)]<-2.2*10^(-16)

library(metap)
istwo <- rep(F, nrow(mat_adj))
toinvert <- rep(F, nrow(mat_adj))
pmerged<-c()
for(col in 1:ncol(mat_adj)){
  pmerged<-c(pmerged, sumlog(two2one(mat_adj[,col], two = istwo, invert = toinvert))$p)
}



library(pheatmap)

colnames(mat_adj)<-colnames(mat)
rownames(mat_adj)<-rownames(mat)
anno_p<-data.frame(mergedp= -log10(pmerged), npositive=colSums(mat_adj<0.05))
rownames(anno_p)<-colnames(mat_adj)

graphics.off()
png("results/2025/hubs_datasets_summary_all.png",res=300, 4000,3000)
pheatmap(-log10(t(mat_adj)), cellwidth=15, cellheight=15, keep.dendro=T, annotation_row = anno_p)
dev.off()



