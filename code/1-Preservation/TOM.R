#Network
TOM_basal=TOMsimilarityFromExpr( t(metabric)[meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal", ], power = 6, networkType = "signed" )


TOM_LumA=TOMsimilarityFromExpr( t(metabric)[meta$NOT_IN_OSLOVAL_Pam50Subtype=="LumA", ], power = 6, networkType = "signed"  )


TOM_LumB=TOMsimilarityFromExpr( t(metabric)[meta$NOT_IN_OSLOVAL_Pam50Subtype=="LumB",], power = 6, networkType = "signed"  )


TOM_Her2=TOMsimilarityFromExpr( t(metabric)[meta$NOT_IN_OSLOVAL_Pam50Subtype=="Her2", ], power = 6, networkType = "signed"  )


TOM_Normal=TOMsimilarityFromExpr( t(metabric)[meta$NOT_IN_OSLOVAL_Pam50Subtype=="Normal",], power = 6, networkType = "signed" )




save(TOM_basal, file = "TOM_basal_tot.RData")
save(TOM_LumA, file = "TOM_LumA_tot.RData")
save(TOM_LumB, file = "TOM_LumB_tot.RData")
save(TOM_Her2, file = "TOM_Her2_tot.RData")
save(TOM_Normal, file = "TOM_Normal_tot.RData")


for(i in 3:length(unique(moduleColors_metabric))){
  
  a<-TOM_basal[moduleColors_metabric==unique(moduleColors_metabric)[i],moduleColors_metabric==unique(moduleColors_metabric)[i]]
  b<-TOM_LumA[moduleColors_metabric==unique(moduleColors_metabric)[i],moduleColors_metabric==unique(moduleColors_metabric)[i]]
  c<-TOM_LumB[moduleColors_metabric==unique(moduleColors_metabric)[i],moduleColors_metabric==unique(moduleColors_metabric)[i]]
  d<-TOM_Her2[moduleColors_metabric==unique(moduleColors_metabric)[i],moduleColors_metabric==unique(moduleColors_metabric)[i]]
  e<-TOM_Normal[moduleColors_metabric==unique(moduleColors_metabric)[i],moduleColors_metabric==unique(moduleColors_metabric)[i]]
  
  
  pdf(paste(unique(moduleColors_metabric)[i],"TOMinsubtypes.pdf"))
  boxplot(log10(as.vector(a)), log10(as.vector(b)), log10(as.vector(c)), log10(as.vector(d)), log10(as.vector(e)),
          names=c("Basal", "LumA", "LumB", "Her2", "Normal"), las=2, main=unique(moduleColors_metabric)[i])
  dev.off()
  #boxplot(log10(as.vector(a)), log10(as.vector(d)), log10(as.vector(e)))
  #t.test(log10(as.vector(a)), log10(as.vector(e)))
  
  #c(mean(log10(as.vector(a))), mean(log10(as.vector(b))), mean(log10(as.vector(c))), mean(log10(as.vector(d))), mean(log10(as.vector(e))))
}

