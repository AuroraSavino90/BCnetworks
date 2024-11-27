#centrality comprende solo i geni in comune tra i due dataset
ssMod<-function(data_orig=metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], centrality=centrality_basal, data, modules){
diff_score<-list()

#si selezionano i dati metabric solo riguardo i geni in comune con il secondo dataset
inboth<-intersect(rownames(centrality), rownames(data))
data_orig<-data_orig[inboth,]
data<-data[inboth,]
ind<-0
for(m in which(unique(centrality$module) %in% modules)){
  ind<-ind+1
  #geni del modulo che sono in entrambi i dataset
    genes<-rownames(centrality)[(centrality$module==unique(centrality$module)[m])]
    genes<-intersect(genes, inboth)
    
    ##per calcolare la reference si usa lo stesso dataset: metabric
    #i ranghi vanno calcolati rispetto a tutti i geni, non solo quelli del modulo
    rank_orig<-apply(data_orig,2,rank)
    
    #si selezionano solo i geni del modulo
    rank_orig<-rank_orig[genes,]
    
    ##dataset con cui ci si sta confrontando
    #i ranghi vanno calcolati rispetto a tutti i geni, non solo quelli del modulo
    rank_data<-apply(data,2,rank)
    
    #si selezionano solo i geni del modulo
    rank_data<-rank_data[genes,]
    
    ###calcolato il riferimento basal
    score_long_ref<-matrix(nrow=(length(genes)^2-length(genes))/2, ncol=ncol(data_orig), byrow=F)
    for(i in 1:ncol(data_orig)){
      ##calcolo indice per il campione i
      matrix_rank1<-matrix(rank_orig[,i],nrow=length(genes), ncol=length(genes), byrow=F)
      matrix_rank2<-matrix(rank_orig[,i],nrow=length(genes), ncol=length(genes), byrow=T)
      
      
      score_long_ref[,i]<-unlist(matrix_rank1[upper.tri(matrix_rank1)]-matrix_rank2[upper.tri(matrix_rank2)])
    }
    
    score_long_ref_mean<-rowMeans(score_long_ref)
    
    #to avoid Inf values, add 0.1 to score_long_ref_mean and score_long
    
    #per ogni campione, si calcolano le distanze tra i ranghi dei geni
    
    #si tengono solo le coppie di geni che hanno una certa variabilità nelle differenze dei ranghi
    ####NO!!!, anzi se sono costanti vuol dire che nella reference è importante che la differenza sia quella...
    #dev<-apply(score_long_ref,1,sd)
   # ind<-which(dev<quantile(dev, 0.05) & dev>0)
   # score_long_ref<-score_long_ref[ind,]
    
    score_long<-matrix(nrow=nrow(score_long_ref), ncol=ncol(data), byrow=F)
    for(i in 1:ncol(data)){
      ##calcolo indice per il campione i
      matrix_rank1<-matrix(rank_data[,i],nrow=length(genes), ncol=length(genes), byrow=F)
      matrix_rank2<-matrix(rank_data[,i],nrow=length(genes), ncol=length(genes), byrow=T)
      
      score_long[,i]<-unlist(matrix_rank1[upper.tri(matrix_rank1)]-matrix_rank2[upper.tri(matrix_rank2)])
    }
    
    ##confronto dei network con la reference
    diff_score[[ind]]<-(score_long-score_long_ref_mean)/(score_long_ref_mean)
    
    rm(score_long)
  rm(score_long_ref)
  gc()
}

names(diff_score)<-modules

return(diff_score)
}
