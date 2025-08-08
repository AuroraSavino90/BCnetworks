library(WGCNA)

#load pre-computed moduletrait relationships

load("results/2025/moduleTraitPvalue_METABRIC.RData")
mtPval_metabric<-moduleTraitPvalue_tot

load("results/2025/moduleTraitPvalue_MAINZ.RData")
mtPval_Mainz<-moduleTraitPvalue_tot

load("results/2025/moduleTraitPvalue_UNT.RData")
mtPval_Unt<-moduleTraitPvalue_tot

load("results/2025/moduleTraitPvalue_UPP.RData")
mtPval_Upp<-moduleTraitPvalue_tot

load("results/2025/moduleTraitPvalue_TRANSBIG.RData")
mtPval_Transbig<-moduleTraitPvalue_tot

load("results/2025/moduleTraitPvalue_VDX.RData")
mtPval_Vdx<-moduleTraitPvalue_tot

load("results/2025/moduleTraitPvalue_NKI.RData")
mtPval_NKI<-moduleTraitPvalue_tot


load("results/2025/moduleTraitCor_METABRIC.RData")
mtCor_metabric<-moduleTraitCor_tot

load("results/2025/moduleTraitCor_MAINZ.RData")
mtCor_Mainz<-moduleTraitCor_tot

load("results/2025/moduleTraitCor_UNT.RData")
mtCor_Unt<-moduleTraitCor_tot

load("results/2025/moduleTraitCor_UPP.RData")
mtCor_Upp<-moduleTraitCor_tot

load("results/2025/moduleTraitCor_TRANSBIG.RData")
mtCor_Transbig<-moduleTraitCor_tot

load("results/2025/moduleTraitCor_VDX.RData")
mtCor_Vdx<-moduleTraitCor_tot

load("results/2025/moduleTraitCor_NKI.RData")
mtCor_NKI<-moduleTraitCor_tot

moduleTraitCor_tot<-cbind(mtCor_NKI,mtCor_Transbig,mtCor_Unt,mtCor_Upp,mtCor_Mainz,mtCor_Vdx)
moduleTraitPvalue_tot<-cbind(mtPval_NKI,mtPval_Transbig,mtPval_Unt,mtPval_Upp,mtPval_Mainz,mtPval_Vdx)
moduleTraitCor_tot<-cbind(mtCor_metabric[rownames(moduleTraitCor_tot),],moduleTraitCor_tot)
moduleTraitPvalue_tot<-cbind(mtPval_metabric[rownames(moduleTraitCor_tot),],moduleTraitPvalue_tot)

#remove Unconnected module
moduleTraitCor_tot<-moduleTraitCor_tot[-1,]
moduleTraitPvalue_tot<-moduleTraitPvalue_tot[-1,]

##order rows based on module's correlation with grade
moduleTraitPvalue_tot<-moduleTraitPvalue_tot[order(rowSums(moduleTraitCor_tot[,grep("grade", colnames(moduleTraitCor_tot))]), decreasing=T),]
moduleTraitCor_tot<-moduleTraitCor_tot[order(rowSums(moduleTraitCor_tot[,grep("grade", colnames(moduleTraitCor_tot))]), decreasing=T),]
# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames = rownames(moduleTraitCor_tot)

moduleTraitCor_tot_grade<-moduleTraitCor_tot[,grep("grade", colnames(moduleTraitCor_tot))]
moduleTraitPvalue_tot_grade<-moduleTraitPvalue_tot[,grep("grade", colnames(moduleTraitPvalue_tot))]


library(metap)
#p merged for high grade
pmerged_hg<-c()
for(row in 1:nrow(moduleTraitPvalue_tot_grade)){
  istwo <- rep(T, ncol(moduleTraitCor_tot_grade))
  toinvert <- ifelse(moduleTraitCor_tot_grade[row,]<0,T,F)
  pmerged_hg<-c(pmerged_hg, sumlog(two2one(moduleTraitPvalue_tot_grade[row,], two = istwo, invert = toinvert))$p)
}
names(pmerged_hg)<-rownames(moduleTraitPvalue_tot_grade)

#p merged for low grade
pmerged_lg<-c()
for(row in 1:nrow(moduleTraitPvalue_tot_grade)){
  istwo <- rep(T, ncol(moduleTraitCor_tot_grade))
  toinvert <- ifelse(moduleTraitCor_tot_grade[row,]>0,T,F)
  pmerged_lg<-c(pmerged_lg, sumlog(two2one(moduleTraitPvalue_tot_grade[row,], two = istwo, invert = toinvert))$p)
}
names(pmerged_lg)<-rownames(moduleTraitPvalue_tot_grade)

gradeavg<-rowMeans(moduleTraitCor_tot_grade)


pdf("results/2025/ModuleGradeRelationship_alldatasets.pdf", width=7, height = 7)
textMatrix = paste(signif(moduleTraitCor_tot_grade, 2), "\n(",
                   signif(moduleTraitPvalue_tot_grade, 1), ")", sep = "");

textMatrix[moduleTraitPvalue_tot_grade<2.2*10^(-16)]<-paste(signif(moduleTraitCor_tot_grade[moduleTraitPvalue_tot_grade<2.2*10^(-16)], 2), "\n(",
                                                            
                                                            
                                                            "<2.2e-16)", sep = "");
dim(textMatrix) = dim(moduleTraitCor_tot_grade)
par(mar = c(6, 14, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor_tot_grade,
               xLabels =c("METABRIC","NKI", "TRANSBIG", "UNT", "UPP", "MAINZ", "VDX"),
               
               yLabels = MEColorNames,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.x=0.6,
               cex.lab.y=0.6,
               zlim = c(-1,1),
               main = paste("Module - grade"))
dev.off();

moduleTraitCor_tot_age<-moduleTraitCor_tot[,grep("age", colnames(moduleTraitCor_tot))]
moduleTraitPvalue_tot_age<-moduleTraitPvalue_tot[,grep("age", colnames(moduleTraitCor_tot))]

#p merged for high age
pmerged_ha<-c()
for(row in 1:nrow(moduleTraitPvalue_tot_age)){
  istwo <- rep(T, ncol(moduleTraitCor_tot_age))
  toinvert <- ifelse(moduleTraitCor_tot_age[row,]<0,T,F)
  pmerged_ha<-c(pmerged_ha, sumlog(two2one(moduleTraitPvalue_tot_age[row,], two = istwo, invert = toinvert))$p)
}
names(pmerged_ha)<-rownames(moduleTraitPvalue_tot_age)

#p merged for low age
pmerged_la<-c()
for(row in 1:nrow(moduleTraitPvalue_tot_age)){
  istwo <- rep(T, ncol(moduleTraitCor_tot_age))
  toinvert <- ifelse(moduleTraitCor_tot_age[row,]>0,T,F)
  pmerged_la<-c(pmerged_la, sumlog(two2one(moduleTraitPvalue_tot_age[row,], two = istwo, invert = toinvert))$p)
}
names(pmerged_la)<-rownames(moduleTraitPvalue_tot_age)

ageavg<-rowMeans(moduleTraitCor_tot_age)

pdf("results/2025/ModuleAgeRelationship_alldatasets.pdf", 7,7)

textMatrix = paste(signif(moduleTraitCor_tot_age, 2), "\n(",
                   signif(moduleTraitPvalue_tot_age, 1), ")", sep = "");

textMatrix[moduleTraitPvalue_tot_age<2.2*10^(-16)]<-paste(signif(moduleTraitCor_tot_age[moduleTraitPvalue_tot_age<2.2*10^(-16)], 2), "\n(",
                                                            
                                                            
                                                            "<2.2e-16)", sep = "");
dim(textMatrix) = dim(moduleTraitCor_tot_age)
par(mar = c(6, 14, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor_tot_age,
               xLabels =c("METABRIC","NKI", "TRANSBIG", "UNT", "UPP", "MAINZ", "VDX"),
               
               yLabels = MEColorNames,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.x=0.6,
               cex.lab.y=0.6,
               zlim = c(-1,1),
               main = paste("Module - age"))
dev.off();

#####survival
load(file="results/2025/pvalue_merge_surv_global.RData")



#p merged for low survival
pmerged_ls<-c()
for(row in 1:nrow(pvalue_merge)){
  istwo <- rep(T, ncol(pvalue_merge))
  toinvert <- ifelse(pvalue_merge[row,]>0,T,F)
  ps<-pvalue_merge[row,]
  ps[ps<0]<- (-ps[ps<0])
  ps<-10^(-ps)
  pmerged_ls<-c(pmerged_ls, sumlog(two2one(ps, two = istwo, invert = toinvert))$p)
}
names(pmerged_ls)<-rownames(pvalue_merge)
pmerged_ls<- pmerged_ls[-which(names(pmerged_ls)=="Unconnected")]

pmerged_hs<-c()
for(row in 1:nrow(pvalue_merge)){
  istwo <- rep(T, ncol(pvalue_merge))
  toinvert <- ifelse(pvalue_merge[row,]<0,T,F)
  ps<-pvalue_merge[row,]
  ps[ps<0]<- (-ps[ps<0])
  ps<-10^(-ps)
  pmerged_hs<-c(pmerged_hs, sumlog(two2one(ps, two = istwo, invert = toinvert))$p)
}
names(pmerged_hs)<-rownames(pvalue_merge)
pmerged_hs<- pmerged_hs[-which(names(pmerged_hs)=="Unconnected")]

library(ggrepel)
df<-data.frame(pgrade= -log10(pmerged_hg), page= -log10(pmerged_la), psurv=-log10(pmerged_ls), grade=gradeavg, age=ageavg, module=rownames(moduleTraitPvalue_tot_age))

png("results/2025/global_p_highgrade.png", res=300, 1500,1500)
ggplot(df, aes(x=grade, y=pgrade, label=module))+geom_point()+geom_label_repel()+theme_classic()
dev.off()

pdf("results/2025/global_p_highgrade.pdf", 5,5)
ggplot(df, aes(x=grade, y=pgrade, label=module))+geom_point()+geom_label_repel()+theme_classic()
dev.off()


png("results/2025/global_p_lowage.png", res=300, 1500,1500)
ggplot(df, aes(x=age, y=page, label=module))+geom_point()+geom_label_repel()+theme_classic()
dev.off()

pdf("results/2025/global_p_lowage.pdf", 5,5)
ggplot(df, aes(x=age, y=page, label=module))+geom_point()+geom_label_repel()+theme_classic()
dev.off()


df$module<-factor(df$module, levels=df$module[(order(df$page+df$pgrade+df$psurv, decreasing=T))])

png("results/2025/global_p_aggressiveness.png", res=300, 2000,2000)
ggplot(df, aes(x=module, y=page+pgrade+psurv))+geom_bar(stat="identity")+ theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("results/2025/global_p_aggressiveness.pdf", 5,5)
ggplot(df, aes(x=module, y=page+pgrade+psurv))+geom_bar(stat="identity")+ theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


df<-data.frame(pgrade= -log10(pmerged_lg), page= -log10(pmerged_ha), psurv=-log10(pmerged_hs), grade=gradeavg, age=ageavg, module=rownames(moduleTraitPvalue_tot_age))

png("results/2025/global_p_lowgrade.png", res=300, 1500,1500)
ggplot(df, aes(x=grade, y=pgrade, label=module))+geom_point()+geom_label_repel()+theme_classic()
dev.off()

pdf("results/2025/global_p_lowgrade.pdf", 5, 5)
ggplot(df, aes(x=grade, y=pgrade, label=module))+geom_point()+geom_label_repel()+theme_classic()
dev.off()

png("results/2025/global_p_highage.png", res=300, 1500,1500)
ggplot(df, aes(x=age, y=page, label=module))+geom_point()+geom_label_repel()+theme_classic()
dev.off()

pdf("results/2025/global_p_highage.pdf", 5,5)
ggplot(df, aes(x=age, y=page, label=module))+geom_point()+geom_label_repel()+theme_classic()
dev.off()

df$module<-factor(df$module, levels=df$module[(order(df$page+df$pgrade+df$psurv, decreasing=T))])

png("results/2025/global_p_LOWaggressiveness.png", res=300, 2000,2000)
ggplot(df, aes(x=module, y=page+pgrade+psurv))+geom_bar(stat="identity")+ theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("results/2025/global_p_LOWaggressiveness.pdf", 5,5)
ggplot(df, aes(x=module, y=page+pgrade+psurv))+geom_bar(stat="identity")+ theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


#################################################
############# Basal
###################################################

load("results/2025/moduleTraitPvalue_METABRIC_basal.RData")
mtPval_metabric<-moduleTraitPvalue_tot

load("results/2025/moduleTraitPvalue_MAINZ_basal.RData")
mtPval_Mainz<-moduleTraitPvalue_tot

load("results/2025/moduleTraitPvalue_UNT_basal.RData")
mtPval_Unt<-moduleTraitPvalue_tot

load("results/2025/moduleTraitPvalue_UPP_basal.RData")
mtPval_Upp<-moduleTraitPvalue_tot

load("results/2025/moduleTraitPvalue_TRANSBIG_basal.RData")
mtPval_Transbig<-moduleTraitPvalue_tot

load("results/2025/moduleTraitPvalue_VDX_basal.RData")
mtPval_Vdx<-moduleTraitPvalue_tot

load("results/2025/moduleTraitPvalue_NKI_basal.RData")
mtPval_NKI<-moduleTraitPvalue_tot


load("results/2025/moduleTraitCor_METABRIC_basal.RData")
mtCor_metabric<-moduleTraitCor_tot

load("results/2025/moduleTraitCor_MAINZ_basal.RData")
mtCor_Mainz<-moduleTraitCor_tot

load("results/2025/moduleTraitCor_UNT_basal.RData")
mtCor_Unt<-moduleTraitCor_tot

load("results/2025/moduleTraitCor_UPP_basal.RData")
mtCor_Upp<-moduleTraitCor_tot

load("results/2025/moduleTraitCor_TRANSBIG_basal.RData")
mtCor_Transbig<-moduleTraitCor_tot

load("results/2025/moduleTraitCor_VDX_basal.RData")
mtCor_Vdx<-moduleTraitCor_tot

load("results/2025/moduleTraitCor_NKI_basal.RData")
mtCor_NKI<-moduleTraitCor_tot



moduleTraitCor_tot<-cbind(mtCor_NKI,mtCor_Transbig,mtCor_Unt,mtCor_Upp,mtCor_Mainz,mtCor_Vdx)
moduleTraitPvalue_tot<-cbind(mtPval_NKI,mtPval_Transbig,mtPval_Unt,mtPval_Upp,mtPval_Mainz,mtPval_Vdx)
moduleTraitCor_tot<-cbind(mtCor_metabric[rownames(moduleTraitCor_tot),],moduleTraitCor_tot)
moduleTraitPvalue_tot<-cbind(mtPval_metabric[rownames(moduleTraitCor_tot),],moduleTraitPvalue_tot)

#remove b_unconnected
moduleTraitCor_tot<-moduleTraitCor_tot[-1,]
moduleTraitPvalue_tot<-moduleTraitPvalue_tot[-1,]


moduleTraitPvalue_tot<-moduleTraitPvalue_tot[order(rowSums(moduleTraitCor_tot[,grep("grade", colnames(moduleTraitCor_tot))]), decreasing=T),]
moduleTraitCor_tot<-moduleTraitCor_tot[order(rowSums(moduleTraitCor_tot[,grep("grade", colnames(moduleTraitCor_tot))]), decreasing=T),]
# Convert numerical lables to colors for labeling of modules in the plot
MEColorNames = rownames(moduleTraitCor_tot)

moduleTraitCor_tot_grade<-moduleTraitCor_tot[,grep("grade", colnames(moduleTraitCor_tot))]
moduleTraitPvalue_tot_grade<-moduleTraitPvalue_tot[,grep("grade", colnames(moduleTraitPvalue_tot))]

moduleTraitCor_tot_age<-moduleTraitCor_tot[,grep("age", colnames(moduleTraitCor_tot))]
moduleTraitPvalue_tot_age<-moduleTraitPvalue_tot[,grep("age", colnames(moduleTraitPvalue_tot))]



pdf("results/2025/ModuleGradeRelationship_alldatasets_basal.pdf",7,7)
textMatrix = paste(signif(moduleTraitCor_tot_grade, 2), "\n(",
                   signif(moduleTraitPvalue_tot_grade, 1), ")", sep = "");

textMatrix[moduleTraitPvalue_tot_grade<2.2*10^(-16)]<-paste(signif(moduleTraitCor_tot_grade[moduleTraitPvalue_tot_grade<2.2*10^(-16)], 2), "\n(",
                                                            
                                                            
                                                            "<2.2e-16)", sep = "");
dim(textMatrix) = dim(moduleTraitCor_tot_grade)
par(mar = c(6, 14, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor_tot_grade,
               xLabels =c("METABRIC","NKI", "TRANSBIG", "UNT", "UPP", "MAINZ", "VDX"),
               
               yLabels = MEColorNames,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.x=0.6,
               cex.lab.y=0.6,
               zlim = c(-1,1),
               main = paste("Module - grade"))
dev.off();




pdf("results/2025/ModuleAgeRelationship_alldatasets_basal.pdf", 7,7)

textMatrix = paste(signif(moduleTraitCor_tot_age, 2), "\n(",
                   signif(moduleTraitPvalue_tot_age, 1), ")", sep = "");

textMatrix[moduleTraitPvalue_tot_age<2.2*10^(-16)]<-paste(signif(moduleTraitCor_tot_age[moduleTraitPvalue_tot_age<2.2*10^(-16)], 2), "\n(",
                                                          
                                                          
                                                          "<2.2e-16)", sep = "");
dim(textMatrix) = dim(moduleTraitCor_tot_age)
par(mar = c(6, 14, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor_tot_age,
               xLabels =c("METABRIC","NKI", "TRANSBIG", "UNT", "UPP", "MAINZ", "VDX"),
               
               yLabels = MEColorNames,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.x=0.6,
               cex.lab.y=0.6,
               zlim = c(-1,1),
               main = paste("Module - age"))
dev.off();


library(metap)
#p merged for high grade
pmerged_hg<-c()
for(row in 1:nrow(moduleTraitPvalue_tot_grade)){
  istwo <- rep(T, ncol(moduleTraitCor_tot_grade))
  toinvert <- ifelse(moduleTraitCor_tot_grade[row,]<0,T,F)
  pmerged_hg<-c(pmerged_hg, sumlog(two2one(moduleTraitPvalue_tot_grade[row,], two = istwo, invert = toinvert))$p)
}
names(pmerged_hg)<-rownames(moduleTraitPvalue_tot_grade)

#p merged for low grade
pmerged_lg<-c()
for(row in 1:nrow(moduleTraitPvalue_tot_grade)){
  istwo <- rep(T, ncol(moduleTraitCor_tot_grade))
  toinvert <- ifelse(moduleTraitCor_tot_grade[row,]>0,T,F)
  pmerged_lg<-c(pmerged_lg, sumlog(two2one(moduleTraitPvalue_tot_grade[row,], two = istwo, invert = toinvert))$p)
}
names(pmerged_lg)<-rownames(moduleTraitPvalue_tot_grade)

gradeavg<-rowMeans(moduleTraitCor_tot_grade)


#p merged for high age
pmerged_ha<-c()
for(row in 1:nrow(moduleTraitPvalue_tot_age)){
  istwo <- rep(T, ncol(moduleTraitCor_tot_age))
  toinvert <- ifelse(moduleTraitCor_tot_age[row,]<0,T,F)
  pmerged_ha<-c(pmerged_ha, sumlog(two2one(moduleTraitPvalue_tot_age[row,], two = istwo, invert = toinvert))$p)
}
names(pmerged_ha)<-rownames(moduleTraitPvalue_tot_age)

#p merged for low age
pmerged_la<-c()
for(row in 1:nrow(moduleTraitPvalue_tot_age)){
  istwo <- rep(T, ncol(moduleTraitCor_tot_age))
  toinvert <- ifelse(moduleTraitCor_tot_age[row,]>0,T,F)
  pmerged_la<-c(pmerged_la, sumlog(two2one(moduleTraitPvalue_tot_age[row,], two = istwo, invert = toinvert))$p)
}
names(pmerged_la)<-rownames(moduleTraitPvalue_tot_age)

ageavg<-rowMeans(moduleTraitCor_tot_age)

#####survival
load(file="results/2025/pvalue_merge_surv_basal.RData")

#p merged for low survival
pmerged_ls<-c()
for(row in 1:nrow(pvalue_merge_b)){
  istwo <- rep(T, ncol(pvalue_merge_b))
  toinvert <- ifelse(pvalue_merge[row,]>0,T,F)
  ps<-pvalue_merge_b[row,]
  ps[ps<0]<- (-ps[ps<0])
  ps<-10^(-ps)
  pmerged_ls<-c(pmerged_ls, sumlog(two2one(ps, two = istwo, invert = toinvert))$p)
}
names(pmerged_ls)<-rownames(pvalue_merge_b)
pmerged_ls<- pmerged_ls[-which(names(pmerged_ls)=="b_Unconnected")]

pmerged_hs<-c()
for(row in 1:nrow(pvalue_merge_b)){
  istwo <- rep(T, ncol(pvalue_merge_b))
  toinvert <- ifelse(pvalue_merge_b[row,]<0,T,F)
  ps<-pvalue_merge_b[row,]
  ps[ps<0]<- (-ps[ps<0])
  ps<-10^(-ps)
  pmerged_hs<-c(pmerged_hs, sumlog(two2one(ps, two = istwo, invert = toinvert))$p)
}
names(pmerged_hs)<-rownames(pvalue_merge_b)
pmerged_hs<- pmerged_hs[-which(names(pmerged_hs)=="b_Unconnected")]

library(ggrepel)
df<-data.frame(pgrade= -log10(pmerged_hg), page= -log10(pmerged_la), psurv=-log10(pmerged_ls), grade=gradeavg, age=ageavg, module=rownames(moduleTraitPvalue_tot_age))

png("results/2025/global_p_highgrade_b.png", res=300, 1500,1500)
ggplot(df, aes(x=grade, y=pgrade, label=module))+geom_point()+geom_label_repel()+theme_classic()
dev.off()

pdf("results/2025/global_p_highgrade_b.pdf", 5,5)
ggplot(df, aes(x=grade, y=pgrade, label=module))+geom_point()+geom_label_repel()+theme_classic()
dev.off()

png("results/2025/global_p_lowage_b.png", res=300, 1500,1500)
ggplot(df, aes(x=age, y=page, label=module))+geom_point()+geom_label_repel()+theme_classic()
dev.off()

pdf("results/2025/global_p_lowage_b.pdf", 5,5)
ggplot(df, aes(x=age, y=page, label=module))+geom_point()+geom_label_repel()+theme_classic()
dev.off()

df$module<-factor(df$module, levels=df$module[(order(df$page+df$pgrade+df$psurv, decreasing=T))])

png("results/2025/global_p_aggressiveness_b.png", res=300, 2000,2000)
ggplot(df, aes(x=module, y=page+pgrade+psurv))+geom_bar(stat="identity")+ theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("results/2025/global_p_aggressiveness_b.pdf", 5,5)
ggplot(df, aes(x=module, y=page+pgrade+psurv))+geom_bar(stat="identity")+ theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

df<-data.frame(pgrade= -log10(pmerged_lg), page= -log10(pmerged_ha), psurv=-log10(pmerged_hs), grade=gradeavg, age=ageavg, module=rownames(moduleTraitPvalue_tot_age))

png("results/2025/global_p_lowgrade_b.png", res=300, 1500,1500)
ggplot(df, aes(x=grade, y=pgrade, label=module))+geom_point()+geom_label_repel()+theme_classic()
dev.off()

pdf("results/2025/global_p_lowgrade_b.pdf", 5,5)
ggplot(df, aes(x=grade, y=pgrade, label=module))+geom_point()+geom_label_repel()+theme_classic()
dev.off()

png("results/2025/global_p_highage_b.png", res=300, 1500,1500)
ggplot(df, aes(x=age, y=page, label=module))+geom_point()+geom_label_repel()+theme_classic()
dev.off()

pdf("results/2025/global_p_highage_b.pdf", 5,5)
ggplot(df, aes(x=age, y=page, label=module))+geom_point()+geom_label_repel()+theme_classic()
dev.off()

df$module<-factor(df$module, levels=df$module[(order(df$page+df$pgrade+df$psurv, decreasing=T))])

png("results/2025/global_p_LOWaggressiveness_b.png", res=300, 2000,2000)
ggplot(df, aes(x=module, y=page+pgrade+psurv))+geom_bar(stat="identity")+ theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("results/2025/global_p_LOWaggressiveness_b.pdf", 5,5)
ggplot(df, aes(x=module, y=page+pgrade+psurv))+geom_bar(stat="identity")+ theme_classic()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


