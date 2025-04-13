library(STRINGdb)
string_db <- STRINGdb$new( version="12", species=9606,
                           score_threshold=200, input_directory="")

load("data/RData/metabric.RData")
load("results/2025/centrality_global.RData")
load("results/2025/centrality_basal.RData")

set.seed(0683650)
scoremod<-c()
scorerand<-matrix(nrow=100, ncol=length(unique(centrality_global$module)))
for(i in 2:length(unique(centrality_global$module))){
df<-data.frame(gene=rownames(centrality_global)[centrality_global$module==unique(centrality_global$module)[i]])
Modmapped <- string_db$map( df, "gene", removeUnmappedRows = TRUE )
int<-string_db$get_interactions(Modmapped$STRING_id)
scoremod<-c(scoremod, sum(int[,3]))
for(j in 1:100){
dfrand<-data.frame(gene=rownames(centrality_global)[sample(1:nrow(centrality_global), length(rownames(centrality_global)[centrality_global$module==unique(centrality_global$module)[i]]), replace = F)])
Randmapped <- string_db$map( dfrand, "gene", removeUnmappedRows = TRUE )
intrand<-string_db$get_interactions(Randmapped$STRING_id)
scorerand[j,i]<-sum(intrand[,3])
}
}
save(scorerand, file="results/2025/scorerand.RData")

library(ggplot2)
df<-data.frame(scorerand=colMeans(scorerand[,-1]), scoremod, module=unique(centrality_global$module)[-1])
png("results/2025/STRING_PPI.png", res=300, 1300, 1300)
ggplot(df, aes(scorerand, scoremod))+geom_point(size=4)+ geom_abline(intercept = 0, slope = 1, color="red", 
                                                               linetype="dashed", size=1.5)+theme_bw()+xlab("STRING score random genes")+ylab("STRING score module genes")
dev.off()

pdf("results/2025/STRING_PPI_labels.pdf", 10, 10)
ggplot(df, aes(scorerand, scoremod, label=module))+geom_point(size=4)+geom_text()+ geom_abline(intercept = 0, slope = 1, color="red", 
                                                                 linetype="dashed", size=1.5)+theme_bw()+xlab("STRING score random genes")+ylab("STRING score module genes")
dev.off()

