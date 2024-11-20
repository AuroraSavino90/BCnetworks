library(DESeq2)

counts<-read.table("data/Poli_E2F3_KO-RNAseq-v1-run241011/RNAseq/dataset/v1-run241011/GEP.count", header = T, row.names = 1)
RPM<-t(t(counts)/colSums(counts))*1000000

counts<-counts[rowSums(counts>=10)>2,]
RPM<-RPM[rownames(counts),]
RPMlog<-log2(RPM+1)


#######################################################
########### DEGs
#########################################################

#clone can't be used as a covariate as we have only one clone for controls, leading to nested variables

############## pooling clones

metadata<-data.frame(KO.gene=rep(c("EV", "E2F3", "E2F3"), each=3),
                     clone=rep(c(100, 1, 3), each=3))

dds468 <- DESeqDataSetFromMatrix(countData = counts,
                                 colData = metadata,
                                 design= ~ KO.gene)
dds468 <- DESeq(dds468)
dds468 <- results(dds468, contrast = c("KO.gene","E2F3", "EV"))

##########################
### Enrichment test
###########################

######network data to load
load("../../R analyses/METABRIC networks/centrality_basal.RData")
load("../../R analyses/METABRIC networks/metabric.RData")
load("../../R analyses/METABRIC networks/meta.RData")

############## pooling clones
counts_name<-counts
ft<-fishertest_alldat(alldat=c("dds468"),
                      names_alldat=c("MDAMB468 E2F3"))

graphics.off()
png(paste("results/",date, "/Enrich_up_E2F3.png", sep=""), res=300, 1500, 2500)
pheatmap(-log10(ft[[1]]),cellwidth=15, cellheight=15)
dev.off()

png(paste("results/",date, "/Enrich_down_E2F3.png", sep=""), res=300, 1500, 2500)
pheatmap(-log10(ft[[2]]),cellwidth=15, cellheight=15)
dev.off()

png(paste("results/",date, "/Enrich_all_E2F3.png", sep=""), res=300, 1500, 2500)
pheatmap(-log10(ft[[3]]),cellwidth=15, cellheight=15)
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
pcaproj_hubs$KO.gene<-factor(pcaproj_hubs$KO.gene, levels=c("EV", "E2F3"))

####plot b_E2F_targets across KOs
png(paste("results/",date, "/bE2F_project_E2F.png", sep=""), res=300, 1500, 1500)
ggplot(pcaproj_hubs, aes(x=KO.gene, y=b_E2F_targets, fill=KO.gene))+geom_boxplot()+theme_classic()
dev.off()


#######compute the cohens'd

cd<-matrix(nrow=1, ncol=19)
  ctrl<-subset(pcaproj_hubs, KO.gene=="EV")
  trt<-subset(pcaproj_hubs, KO.gene=="E2F3")
  
  for(i in 1:19){
    cd[1,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
  }
  
  colnames(cd)<-colnames(pcaproj_hubs)[1:19]

  ########## plot changes in MEs (Cohen's d) vs modules' correlation with b_E2F_targets
  ## for each KO, the most affected modules are either the most highly or lowly correlated with b_E2F_targets
  
  df<-data.frame(cd=c(t(cd)), condition=rep(c("MDAMB468 E2F3"), each= 19),corr=rep(cc["b_E2F_targets",colnames(cd)],1),
                 module=rep(colnames(cd),1))
  
  png(paste("results/",date, "/CohenVScorr_bE2F_E2F3.png", sep=""),res=300, 4500, 2000)
  ggplot(df, aes(x=corr, y=cd, label=module))+geom_point(size=2)+geom_smooth(method = lm)+stat_cor(label.x=-0.5, label.y = 17)+geom_text_repel(max.overlaps = 5)+theme_bw()+theme(strip.text=element_text(size = 12, face = "bold"))
  dev.off()
  

  
  cd<-matrix(nrow=2, ncol=19)
  ctrl<-subset(pcaproj_hubs, KO.gene=="EV"&clone%in% c(100,1))
  trt<-subset(pcaproj_hubs, KO.gene=="E2F3"&clone%in% c(100,1))
  
  for(i in 1:19){
    cd[1,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
  }
  
  ctrl<-subset(pcaproj_hubs, KO.gene=="EV"&clone%in% c(100,3))
  trt<-subset(pcaproj_hubs, KO.gene=="E2F3"&clone%in% c(100,3))
  
  for(i in 1:19){
    cd[2,i]<-(mean(trt[,i])-mean(ctrl[,i]))/sqrt((var(trt[,i])+var(ctrl[,i]))/2)
  }
  
  colnames(cd)<-colnames(pcaproj_hubs)[1:19]
  
  df<-data.frame(cd=c(t(cd)), condition=rep(c("MDAMB468 E2F3 1", "MDAMB468 E2F3 3"), each= 19),corr=rep(cc["b_E2F_targets",colnames(cd)],1),
                 module=rep(colnames(cd),1))
  
  
  png(paste("results/",date, "/CohenVScorr_bE2F_E2F_clones.png", sep=""),res=300, 6000, 2000)
  ggplot(df, aes(x=corr, y=cd, label=module))+geom_point(size=2)+facet_grid(~condition)+geom_smooth(method = lm)+stat_cor(label.x=-0.5, label.y = 30)+geom_text_repel(max.overlaps = 5)+
    scale_y_continuous(limits = c(-30, 30))+
    theme_bw()+theme(strip.text=element_text(size = 12, face = "bold"))
  dev.off()
  
