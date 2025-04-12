###Check surv data!! Ci sono degli zeri...

library(ggplot2)

load("results/2025/centrality_basal.RData")
load("results/2025/centrality_global.RData")
load("results/2025/surv_data.RData")

surv_data_global_5years_poor<-surv_data_global_5years
surv_data_global_5years_poor$pval[surv_data_global_5years$direction=="Good"]<-1
surv_data_global_5years_poor$pval[surv_data_global_5years_poor$pval==1]<-NA
surv_data_global_5years_poor$pval[which(surv_data_global_5years_poor$pval==0)]<-2.2*10^(-16)

for(i in 1:length(unique(centrality_global$module))){
  if(length(na.omit(-log10(surv_data_global_5years_poor$pval[centrality_global$module==unique(centrality_global$module)[i]])))>0){
png(paste("results/2025/DSS5years_global_",unique(centrality_global$module)[i],"_poor.png", sep=""), res=300, 1100, 1100)

df<-data.frame(kWithin=centrality_global$kWithin[centrality_global$module==unique(centrality_global$module)[i]] ,
               pvalue=-log10(surv_data_global_5years_poor$pval[centrality_global$module==unique(centrality_global$module)[i]]))
df$pvalue[which(df$pvalue>15.65758)]<-15.65758
p<-ggplot(df, aes(x=kWithin, y=pvalue))+geom_point()+geom_smooth(method="lm")+theme_bw()+ggtitle(label = unique(centrality_global$module)[i])+ylab(label = "-log10(p-value)")
print(p)
dev.off()
}
}


tests_poor_pval<-list()
tests_poor_dir<-list()
for(i in 1:length(unique(centrality_global$module))){
  if(length(na.omit(-log10(surv_data_global_5years_poor$pval[centrality_global$module==unique(centrality_global$module)[i]])))>0){
  kW<-centrality_global$kWithin[centrality_global$module==unique(centrality_global$module)[i]]
  pval<-(-log10(surv_data_global_5years_poor$pval[centrality_global$module==unique(centrality_global$module)[i]]))
  pval[which(pval>15.65758)]<-15.65758
  
   tests_poor_pval[[i]]<-cor.test(kW , pval)[[3]]
  tests_poor_dir[[i]]<-cor.test(kW , pval)[[4]]
  
  }
 
}

names(tests_poor_pval)<-unique(centrality_global$module)
names(tests_poor_dir)<-unique(centrality_global$module)


###########
surv_data_global_5years_good<-surv_data_global_5years
surv_data_global_5years_good$pval[surv_data_global_5years$direction=="Poor"]<-1
surv_data_global_5years_good$pval[surv_data_global_5years_good$pval==1]<-NA
surv_data_global_5years_good$pval[which(surv_data_global_5years_good$pval==0)]<-2.2*10^(-16)

for(i in 1:length(unique(centrality_global$module))){
  if(length(na.omit(-log10(surv_data_global_5years_good$pval[centrality_global$module==unique(centrality_global$module)[i]])))>0){
    png(paste("results/2025/DSS5years_global_",unique(centrality_global$module)[i],"_good.png", sep=""), res=300, 1100, 1100)
    df<-data.frame(kWithin=centrality_global$kWithin[centrality_global$module==unique(centrality_global$module)[i]] ,
                   pvalue=-log10(surv_data_global_5years_good$pval[centrality_global$module==unique(centrality_global$module)[i]]))
    df$pvalue[which(df$pvalue>15.65758)]<-15.65758
    p<-ggplot(df, aes(x=kWithin, y=pvalue))+geom_point()+geom_smooth(method="lm")+theme_bw()+ggtitle(label = unique(centrality_global$module)[i])+ylab(label = "-log10(p-value)")
    print(p)
    dev.off()
  }
}


tests_good_pval<-list()
tests_good_dir<-list()
for(i in 21:21){
  if(length(na.omit(-log10(surv_data_global_5years_good$pval[centrality_global$module==unique(centrality_global$module)[i]])))>0){
    kW<-centrality_global$kWithin[centrality_global$module==unique(centrality_global$module)[i]]
    pval<-(-log10(surv_data_global_5years_good$pval[centrality_global$module==unique(centrality_global$module)[i]]))
    pval[which(pval>15.65758)]<-15.65758
    
    tests_good_pval[[i]]<-cor.test(kW , pval)[[3]]
    tests_good_dir[[i]]<-cor.test(kW , pval)[[4]]
    
  }
  
}

names(tests_good_pval)<-unique(centrality_global$module)
names(tests_good_dir)<-unique(centrality_global$module)

tests_poor_pval[[21]]<-1
tests_poor_pval[[4]]<-1
tests_poor_pval[[11]]<-1
tests_poor_pval[[19]]<-1

tests_poor_dir[[21]]<-0
tests_poor_dir[[4]]<-0
tests_poor_dir[[11]]<-0
tests_poor_dir[[19]]<-0

tests_good_pval[[6]]<-1
tests_good_pval[[7]]<-1
tests_good_pval[[10]]<-1
tests_good_pval[[16]]<-1
tests_good_pval[[18]]<-1
tests_good_pval[[20]]<-1

tests_good_dir[[6]]<-0
tests_good_dir[[7]]<-0
tests_good_dir[[10]]<-0
tests_good_dir[[16]]<-0
tests_good_dir[[18]]<-0
tests_good_dir[[20]]<-0

pval_mat<-matrix(nrow=21, ncol=2)
rownames(pval_mat)<-names(tests_good_pval)
pval_mat[,1]<-unlist(tests_poor_pval)
pval_mat[,2]<-unlist(tests_good_pval)

cor_mat<-matrix(nrow=21, ncol=2)
rownames(cor_mat)<-names(tests_good_dir)
cor_mat[,1]<-unlist(tests_poor_dir)
cor_mat[,2]<-unlist(tests_good_dir)

pval_mat[pval_mat>0.05]<-1
pval_mat<-pval_mat[-1,]
cor_mat<-cor_mat[-1,]

#textMatrix[moduleTraitPvalue_tot<2.2*10^(-16)]<-paste(signif(moduleTraitCor_tot[moduleTraitPvalue_tot<2.2*10^(-16)], 2), "\n(",
#                                                      "<2.2e-16)", sep = "");


pval_mat<-pval_mat[order((cor_mat[,1]-cor_mat[,2]), decreasing=T),]
cor_mat<-cor_mat[order((cor_mat[,1]-cor_mat[,2]), decreasing=T),]

textMatrix = paste(signif(cor_mat, 2), "\n(",
                   signif(pval_mat, 1), ")", sep = "");

png("Correlation_kWithinSurvival_metabric_global.png", res = 450, width=2500, height = 4000)
par(mar = c(6, 14, 3, 2.2));
labeledHeatmap(Matrix = cor_mat,
               xLabels = c("Poor prognosis", "Good prognosis"),
               yLabels = rownames(pval_mat),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = matrix(textMatrix, nrow=20),
               setStdMargins = FALSE,
               cex.text = 0.8,
               cex.lab.x=1,
               cex.lab.y=1,
               zlim = c(-1,1),
               main = paste("Correlation kWithin - survival"))
dev.off()


cor_mat[pval_mat>0.05]<-NA

pheatmap(cor_mat, cluster_rows = F, cluster_cols = F)

paletteLength <- 50
# use floor and ceiling to deal with even/odd length pallettelengths
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(unlist(cor_mat), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(unlist(cor_mat), na.rm=T)/paletteLength, max(unlist(cor_mat), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

png("Correlation_kWithinSurvival_metabric_global_pheat.png", res=300, 1200, 2000)
pheatmap(cor_mat[order(rowSums(cbind(cor_mat[,1],-cor_mat[,2]), na.rm=T), decreasing = T),],  cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F)
dev.off()

