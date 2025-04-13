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



png("results/2025/ModuleGradeRelationship_alldatasets.png", res = 300, width=2000, height = 2000)
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

png("results/2025/ModuleAgeRelationship_alldatasets.png", res = 300, width=2000, height = 2000)

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



png("results/2025/ModuleGradeRelationship_alldatasets_basal.png", res = 300, width=2000, height = 2000)
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




png("results/2025/ModuleAgeRelationship_alldatasets_basal.png", res = 300, width=2000, height = 2000)

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


