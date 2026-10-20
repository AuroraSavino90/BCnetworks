library(ggplot2)
library(pheatmap)
dir.create("results/2025_fix", showWarnings=F)

load(file="results/2025/centrality_global.RData")
load(file="results/2025/centrality_basal.RData")
load(file="results/2025_fix/surv_data.RData")

###relationship between kWithin and -log10(p) of the association with survival (log-rank, 5-year DSS)
# separately for genes associated with poor and with good prognosis (direction from 1_surv_each_gene.R)
# Spearman correlation, no cap on -log10(p), BH correction across all tests (modules x 2)
# genes are matched by name between centrality and survival tables

save_plot<-function(p, file){
  png(paste(file, ".png", sep=""), res=300, 1100, 1100)
  print(p)
  dev.off()
  pdf(paste(file, ".pdf", sep=""), 5,5)
  print(p)
  dev.off()
}

kWithin_vs_survival<-function(centrality, surv_tab, unconnected, label){
  centrality<-as.data.frame(centrality)
  centrality$kWithin<-as.numeric(centrality$kWithin)
  if(any(surv_tab$pval==0, na.rm=T)) warning("p-values equal to 0: -log10(p) is infinite for these genes")
  genes<-intersect(rownames(centrality), rownames(surv_tab))
  message(label, ": ", length(genes), "/", nrow(centrality), " network genes with survival data")
  modules<-setdiff(unique(centrality$module), unconnected)
  res<-lapply(modules, function(m){
    out<-data.frame(module=m)
    for(dir in c("Poor", "Good")){
      g<-genes[centrality[genes, "module"]==m & surv_tab[genes, "direction"]==dir]
      df<-data.frame(kWithin=centrality[g, "kWithin"], pvalue=-log10(surv_tab[g, "pval"]))
      df<-df[complete.cases(df),]
      out[, paste("n_genes", tolower(dir), sep="_")]<-nrow(df)
      out[, paste("rho", tolower(dir), sep="_")]<-NA
      out[, paste("pval", tolower(dir), sep="_")]<-NA
      if(nrow(df)<10) next
      
      p<-ggplot(df, aes(x=kWithin, y=pvalue))+geom_point()+geom_smooth(method="lm")+theme_bw()+ggtitle(label = m)+ylab(label = paste(dir, "prognosis\n-log10(p-value)"))
      save_plot(p, paste("results/2025_fix/DSS5years_", label, "_", m, "_", tolower(dir), sep=""))
      
      ct<-cor.test(df$kWithin, df$pvalue, method="spearman", exact=F)
      out[, paste("rho", tolower(dir), sep="_")]<-unname(ct$estimate)
      out[, paste("pval", tolower(dir), sep="_")]<-ct$p.value
    }
    out
  })
  res<-do.call(rbind, res)
  fdr<-p.adjust(c(res$pval_poor, res$pval_good), method="BH")
  res$FDR_poor<-fdr[1:nrow(res)]
  res$FDR_good<-fdr[(nrow(res)+1):(2*nrow(res))]
  res<-res[order(rowSums(cbind(res$rho_poor, -res$rho_good), na.rm=T), decreasing=T),]
  write.csv(res, paste("results/2025_fix/Correlation_kWithinSurvival_metabric_", label, ".csv", sep=""), row.names=F)
  message(label, ": ", sum(res$FDR_poor<0.05, na.rm=T), "/", nrow(res), " modules with FDR<0.05 (poor); ",
          sum(res$FDR_good<0.05, na.rm=T), "/", nrow(res), " (good); ",
          sum(res$FDR_poor<0.05 | res$FDR_good<0.05, na.rm=T), "/", nrow(res), " (either)")
  res
}

plot_heatmap<-function(res, label, height){
  cor_mat<-cbind("Poor prognosis"=res$rho_poor, "Good prognosis"=res$rho_good)
  rownames(cor_mat)<-res$module
  fdr_mat<-cbind(res$FDR_poor, res$FDR_good)
  cor_mat[is.na(fdr_mat) | fdr_mat>=0.05]<-NA
  if(all(is.na(cor_mat))) return(invisible(NULL))
  # module order as in the original figure (computed on the significant correlations only)
  cor_mat<-cor_mat[order(rowSums(cbind(cor_mat[,1], -cor_mat[,2]), na.rm=T), decreasing=T),]
  paletteLength <- 50
  myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
  # symmetric breaks: equal colour intensity for equal |rho|
  M <- max(abs(cor_mat), na.rm=T)
  myBreaks <- seq(-M, M, length.out=paletteLength + 1)
  
  png(paste("results/2025_fix/Correlation_kWithinSurvival_metabric_", label, "_pheat.png", sep=""), res=300, 1200, height)
  pheatmap(cor_mat, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F)
  dev.off()
  pdf(paste("results/2025_fix/Correlation_kWithinSurvival_metabric_", label, "_pheat.pdf", sep=""), 10,10)
  pheatmap(cor_mat, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, cluster_rows = F, cluster_cols = F)
  dev.off()
}

##############################
##### Global modules
##############################
res_global<-kWithin_vs_survival(centrality_global, surv_data_global_5years, "Unconnected", "global")
plot_heatmap(res_global, "global", 2000)

###################################
######## Basal modules
######################################
res_basal<-kWithin_vs_survival(centrality_basal, surv_data_Basal_5years, "b_Unconnected", "basal")
plot_heatmap(res_basal, "basal", 2000)