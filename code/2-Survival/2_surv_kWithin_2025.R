library(ggplot2)

load(file="results/2025/centrality_global.RData")
load(file="results/2025/centrality_basal.RData")
load(file="results/2025/surv_data.RData")

##############################
##### Global modules
##############################

surv_data_global_5years_poor<-surv_data_global_5years
surv_data_global_5years_poor$pval[surv_data_global_5years$direction=="Good"]<-NA
surv_data_global_5years_poor$pval[which(surv_data_global_5years_poor$pval==0)]<-2.2*10^(-16)

###represent the relationship between kwithin and significance of the relationship with survival (poor)
# module by module

for(i in 1:length(unique(centrality_global$module))){
  if(length(na.omit(-log10(surv_data_global_5years_poor$pval[centrality_global$module==unique(centrality_global$module)[i]])))>=10){
png(paste("results/2025/DSS5years_global_",unique(centrality_global$module)[i],"_poor.png", sep=""), res=300, 1100, 1100)

df<-data.frame(kWithin=centrality_global$kWithin[centrality_global$module==unique(centrality_global$module)[i]] ,
               pvalue=-log10(surv_data_global_5years_poor$pval[centrality_global$module==unique(centrality_global$module)[i]]))
df$pvalue[which(df$pvalue> -log10(2.2*10^(-16)))]<- -log10(2.2*10^(-16))
p<-ggplot(df, aes(x=kWithin, y=pvalue))+geom_point()+geom_smooth(method="lm")+theme_bw()+ggtitle(label = unique(centrality_global$module)[i])+ylab(label = "-log10(p-value)")
print(p)
dev.off()
}
}

###test the relationship between kwithin and significance of the relationship with survival (poor)
# module by module

tests_poor_pval<-as.list(rep(NA, length(unique(centrality_global$module))))
tests_poor_dir<- as.list(rep(NA, length(unique(centrality_global$module))))

for(i in 1:length(unique(centrality_global$module))){
  if(length(na.omit(-log10(surv_data_global_5years_poor$pval[centrality_global$module==unique(centrality_global$module)[i]])))>=10){
  
    kW<-centrality_global$kWithin[centrality_global$module==unique(centrality_global$module)[i]]
  pval<-(-log10(surv_data_global_5years_poor$pval[centrality_global$module==unique(centrality_global$module)[i]]))
  pval[which(pval> -log10(2.2*10^(-16)))]<- -log10(2.2*10^(-16))
  
  tests_poor_pval[[i]]<-cor.test(kW , pval)[[3]]
  tests_poor_dir[[i]]<-cor.test(kW , pval)[[4]]
  
  }
 
}

names(tests_poor_pval)<-unique(centrality_global$module)
names(tests_poor_dir)<-unique(centrality_global$module)


###########
surv_data_global_5years_good<-surv_data_global_5years
surv_data_global_5years_good$pval[surv_data_global_5years$direction=="Poor"]<-NA
surv_data_global_5years_good$pval[which(surv_data_global_5years_good$pval==0)]<-2.2*10^(-16)

###represent the relationship between kwithin and significance of the relationship with survival (good)
# module by module

for(i in 1:length(unique(centrality_global$module))){
  if(length(na.omit(-log10(surv_data_global_5years_good$pval[centrality_global$module==unique(centrality_global$module)[i]])))>=10){
    png(paste("results/2025/DSS5years_global_",unique(centrality_global$module)[i],"_good.png", sep=""), res=300, 1100, 1100)
    df<-data.frame(kWithin=centrality_global$kWithin[centrality_global$module==unique(centrality_global$module)[i]] ,
                   pvalue=-log10(surv_data_global_5years_good$pval[centrality_global$module==unique(centrality_global$module)[i]]))
    df$pvalue[which(df$pvalue> -log10(2.2*10^(-16)))]<- -log10(2.2*10^(-16))
    p<-ggplot(df, aes(x=kWithin, y=pvalue))+geom_point()+geom_smooth(method="lm")+theme_bw()+ggtitle(label = unique(centrality_global$module)[i])+ylab(label = "-log10(p-value)")
    print(p)
    dev.off()
  }
}

###test the relationship between kwithin and significance of the relationship with survival (poor)
# module by module

tests_good_pval<-as.list(rep(NA, length(unique(centrality_global$module))))
tests_good_dir<-as.list(rep(NA, length(unique(centrality_global$module))))
for(i in 1:length(unique(centrality_global$module))){
  if(length(na.omit(-log10(surv_data_global_5years_good$pval[centrality_global$module==unique(centrality_global$module)[i]])))>=10){
    kW<-centrality_global$kWithin[centrality_global$module==unique(centrality_global$module)[i]]
    pval<-(-log10(surv_data_global_5years_good$pval[centrality_global$module==unique(centrality_global$module)[i]]))
    pval[which(pval> -log10(2.2*10^(-16)))]<- -log10(2.2*10^(-16))
    
    tests_good_pval[[i]]<-cor.test(kW , pval)[[3]]
    tests_good_dir[[i]]<-cor.test(kW , pval)[[4]]
    
  }
  
}

names(tests_good_pval)<-unique(centrality_global$module)
names(tests_good_dir)<-unique(centrality_global$module)



pval_mat<-matrix(nrow=length(unique(centrality_global$module)), ncol=2)
rownames(pval_mat)<-names(tests_good_pval)
pval_mat[,1]<-unlist(tests_poor_pval)
pval_mat[,2]<-unlist(tests_good_pval)

cor_mat<-matrix(nrow=length(unique(centrality_global$module)), ncol=2)
rownames(cor_mat)<-names(tests_good_dir)
cor_mat[,1]<-unlist(tests_poor_dir)
cor_mat[,2]<-unlist(tests_good_dir)

cor_mat[pval_mat>0.05]<-NA
pval_mat<-pval_mat[-which(rownames(pval_mat)=="Unconnected"),]
cor_mat<-cor_mat[-which(rownames(cor_mat)=="Unconnected"),]


library(pheatmap)
paletteLength <- 50
# use floor and ceiling to deal with even/odd length pallettelengths
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(unlist(cor_mat), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(unlist(cor_mat), na.rm=T)/paletteLength, max(unlist(cor_mat), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

png("results/2025/Correlation_kWithinSurvival_metabric_global_pheat.png", res=300, 1200, 2000)
pheatmap(cor_mat[order(rowSums(cbind(cor_mat[,1],-cor_mat[,2]), na.rm=T), decreasing = T),],  cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F)
dev.off()

###################################
######## Basal modules
######################################

surv_data_basal_5years_poor<-surv_data_Basal_5years
surv_data_basal_5years_poor$pval[surv_data_Basal_5years$direction=="Good"]<-NA
surv_data_basal_5years_poor$pval[which(surv_data_basal_5years_poor$pval==0)]<-2.2*10^(-16)

###represent the relationship between kwithin and significance of the relationship with survival (poor)
# module by module

for(i in 1:length(unique(centrality_global$module))){
  if(length(na.omit(-log10(surv_data_basal_5years_poor$pval[centrality_basal$module==unique(centrality_basal$module)[i]])))>=10){
    png(paste("results/2025/DSS5years_basal_",unique(centrality_basal$module)[i],"_poor.png", sep=""), res=300, 1100, 1100)
    
    df<-data.frame(kWithin=centrality_basal$kWithin[centrality_basal$module==unique(centrality_basal$module)[i]] ,
                   pvalue=-log10(surv_data_basal_5years_poor$pval[centrality_basal$module==unique(centrality_basal$module)[i]]))
    df$pvalue[which(df$pvalue> -log10(2.2*10^(-16)))]<- -log10(2.2*10^(-16))
    p<-ggplot(df, aes(x=kWithin, y=pvalue))+geom_point()+geom_smooth(method="lm")+theme_bw()+ggtitle(label = unique(centrality_basal$module)[i])+ylab(label = "-log10(p-value)")
    print(p)
    dev.off()
  }
}

###test the relationship between kwithin and significance of the relationship with survival (poor)
# module by module

tests_poor_pval<-as.list(rep(NA, length(unique(centrality_basal$module))))
tests_poor_dir<- as.list(rep(NA, length(unique(centrality_basal$module))))

for(i in 1:length(unique(centrality_basal$module))){
  if(length(na.omit(-log10(surv_data_basal_5years_poor$pval[centrality_basal$module==unique(centrality_basal$module)[i]])))>=10){
    
    kW<-centrality_basal$kWithin[centrality_basal$module==unique(centrality_basal$module)[i]]
    pval<-(-log10(surv_data_basal_5years_poor$pval[centrality_basal$module==unique(centrality_basal$module)[i]]))
    pval[which(pval> -log10(2.2*10^(-16)))]<- -log10(2.2*10^(-16))
    
    tests_poor_pval[[i]]<-cor.test(kW , pval)[[3]]
    tests_poor_dir[[i]]<-cor.test(kW , pval)[[4]]
    
  }
  
}

names(tests_poor_pval)<-unique(centrality_basal$module)
names(tests_poor_dir)<-unique(centrality_basal$module)


###########
surv_data_basal_5years_good<-surv_data_Basal_5years
surv_data_basal_5years_good$pval[surv_data_Basal_5years$direction=="Poor"]<-NA
surv_data_basal_5years_good$pval[which(surv_data_basal_5years_good$pval==0)]<-2.2*10^(-16)

###represent the relationship between kwithin and significance of the relationship with survival (good)
# module by module

for(i in 1:length(unique(centrality_basal$module))){
  if(length(na.omit(-log10(surv_data_basal_5years_good$pval[centrality_basal$module==unique(centrality_basal$module)[i]])))>=10){
    png(paste("results/2025/DSS5years_basal_",unique(centrality_basal$module)[i],"_good.png", sep=""), res=300, 1100, 1100)
    df<-data.frame(kWithin=centrality_basal$kWithin[centrality_basal$module==unique(centrality_basal$module)[i]] ,
                   pvalue=-log10(surv_data_basal_5years_good$pval[centrality_basal$module==unique(centrality_basal$module)[i]]))
    df$pvalue[which(df$pvalue> -log10(2.2*10^(-16)))]<- -log10(2.2*10^(-16))
    p<-ggplot(df, aes(x=kWithin, y=pvalue))+geom_point()+geom_smooth(method="lm")+theme_bw()+ggtitle(label = unique(centrality_basal$module)[i])+ylab(label = "-log10(p-value)")
    print(p)
    dev.off()
  }
}

###test the relationship between kwithin and significance of the relationship with survival (poor)
# module by module

tests_good_pval<-as.list(rep(NA, length(unique(centrality_basal$module))))
tests_good_dir<-as.list(rep(NA, length(unique(centrality_basal$module))))
for(i in 1:length(unique(centrality_basal$module))){
  if(length(na.omit(-log10(surv_data_basal_5years_good$pval[centrality_basal$module==unique(centrality_basal$module)[i]])))>=10){
    kW<-centrality_basal$kWithin[centrality_basal$module==unique(centrality_basal$module)[i]]
    pval<-(-log10(surv_data_basal_5years_good$pval[centrality_basal$module==unique(centrality_basal$module)[i]]))
    pval[which(pval> -log10(2.2*10^(-16)))]<- -log10(2.2*10^(-16))
    
    tests_good_pval[[i]]<-cor.test(kW , pval)[[3]]
    tests_good_dir[[i]]<-cor.test(kW , pval)[[4]]
    
  }
  
}

names(tests_good_pval)<-unique(centrality_basal$module)
names(tests_good_dir)<-unique(centrality_basal$module)



pval_mat<-matrix(nrow=length(unique(centrality_basal$module)), ncol=2)
rownames(pval_mat)<-names(tests_good_pval)
pval_mat[,1]<-unlist(tests_poor_pval)
pval_mat[,2]<-unlist(tests_good_pval)

cor_mat<-matrix(nrow=length(unique(centrality_basal$module)), ncol=2)
rownames(cor_mat)<-names(tests_good_dir)
cor_mat[,1]<-unlist(tests_poor_dir)
cor_mat[,2]<-unlist(tests_good_dir)

cor_mat[pval_mat>0.05]<-NA
pval_mat<-pval_mat[-which(rownames(pval_mat)=="b_Unconnected"),]
cor_mat<-cor_mat[-which(rownames(cor_mat)=="b_Unconnected"),]


library(pheatmap)
paletteLength <- 50
# use floor and ceiling to deal with even/odd length pallettelengths
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(unlist(cor_mat), na.rm=T), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(unlist(cor_mat), na.rm=T)/paletteLength, max(unlist(cor_mat), na.rm=T), length.out=floor(paletteLength/2)))
length(myBreaks) == length(paletteLength) + 1

png("results/2025/Correlation_kWithinSurvival_metabric_basal_pheat.png", res=300, 1500, 2000)
pheatmap(cor_mat[order(rowSums(cbind(cor_mat[,1],-cor_mat[,2]), na.rm=T), decreasing = T),],  cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F)
dev.off()

