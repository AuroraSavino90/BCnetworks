load("data/Alldegrees_metabric_global.RData")
centrality_global <- Alldegrees
load("data/Alldegrees_basal.RData")
centrality_basal <- Alldegrees_basal

library(WGCNA)
moduleColors = labels2colors(net_metabric$colors)
moduleColors_basal = labels2colors(net_metabric_basal$colors)

#rename modules
alt_names<-c("Unconnected",
             "Estrogen_response",
             "NC1",
             "Allograft_rejection",
             "E2F_targets",
             "KRAS_dn1",
             "TNFA_signaling_via_NFKB",
             "MYC_targets1",
             "KRAS_dn2",
             "Oxidative_phosphorylation1",
             "Oxidative_phosphorylation2",
             "EMT",
             "MYC_targets2",
             "NC2",
             "NC3",
             "Pancreas_beta_cells",
             "NC4",
             "NC5",
             "Adipogenesis",
             "Interpheron_alpha",
             "NC6")

alt_names_b<-c("b_Unconnected",
               "b_NC1",
               "b_E2F_targets",
               "b_Allograft_rejection",
               "b_EMT1",
               "b_KRAS_dn",
               "b_Pancreas_beta_cells",
               "b_Myc_targets1",
               "b_Oxidative_phosphorylation1",
               "b_Myc_targets2",
               "b_Inflammatory_response",
               "b_Oxidative_phosphorylation2",
               "b_Estrogen_response",
               "b_Estrogen_response_early",
               "b_Myogenesis",
               "b_NC2",
               "b_Interferon_alpha", 
               "b_EMT2",
               "b_NC3",
               "b_Interferon_gamma",
               "b_NC4")


moduleColors<-alt_names[match(moduleColors, labels2colors(0:20))]

centrality_global <- Alldegrees
centrality_global$module<-moduleColors

rank<-rep(0, length(moduleColors))
for(i in 1:length(unique(moduleColors))){
  rank[moduleColors==unique(moduleColors)[i]]<-rank(-centrality_global$kWithin[moduleColors==unique(moduleColors)[i]])
  
}

centrality_global$rank<-rank

save(centrality_global, file="centrality_global.RData")



moduleColors_basal<-alt_names_b[match(moduleColors_basal, labels2colors(0:20))]

centrality_basal <- Alldegrees_basal
centrality_basal$module<-moduleColors_basal

rank<-rep(0, length(moduleColors_basal))
for(i in 1:length(unique(moduleColors_basal))){
  rank[moduleColors_basal==unique(moduleColors_basal)[i]]<-rank(-centrality_basal$kWithin[moduleColors_basal==unique(moduleColors_basal)[i]])
  
}

centrality_basal$rank<-rank

save(centrality_basal, file="centrality_basal.RData")
