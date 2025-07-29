library(WGCNA)

load("data/RData/meta.RData")
load("data/RData/metabric.RData")
load("results/2025/centrality_global.RData")
load("data/Complete_METABRIC_Clinical_Survival_Data__DSS.rbin")
MEs= moduleEigengenes(t(metabric), centrality_global$module)$eigengenes
colnames(MEs)<-gsub("^ME","", colnames(MEs))


p_sub<-c()
for(i in 1:ncol(MEs)){p_sub<-c(p_sub, (summary(lm(as.numeric(meta$grade)~MEs[,i]+meta$NOT_IN_OSLOVAL_Pam50Subtype))[[4]][2,4]))}
coef_sub<-c()
for(i in 1:ncol(MEs)){coef_sub<-c(coef_sub, (summary(lm(as.numeric(meta$grade)~MEs[,i]+meta$NOT_IN_OSLOVAL_Pam50Subtype))[[4]][2,1]))}

p<-c()
for(i in 1:ncol(MEs)){p<-c(p, (summary(lm(as.numeric(meta$grade)~MEs[,i]))[[4]][2,4]))}

plot(-log10(p), -log10(p_sub))

p_sub_age<-c()
for(i in 1:ncol(MEs)){p_sub_age<-c(p_sub_age, (summary(lm(as.numeric(meta$age_at_diagnosis)~MEs[,i]+meta$NOT_IN_OSLOVAL_Pam50Subtype))[[4]][2,4]))}
coef_sub_age<-c()
for(i in 1:ncol(MEs)){coef_sub_age<-c(coef_sub_age, (summary(lm(as.numeric(meta$age)~MEs[,i]+meta$NOT_IN_OSLOVAL_Pam50Subtype))[[4]][2,1]))}

surv<-Complete_METABRIC_Clinical_Survival_Data__DSS
surv[surv[,1]>1825,"status"]<-0
surv[surv[,1]>1825,"time"]<-1825

library(survminer)
library(survival)
p_sub_surv<-c()
for(i in 1:ncol(MEs)){
  p_sub_surv<-c(p_sub_surv, summary(coxph(Surv(surv[,1], surv[,2])~MEs[,i]+meta$NOT_IN_OSLOVAL_Pam50Subtype))[[7]][1,5])
}


coef_sub_surv<-c()
for(i in 1:ncol(MEs)){
  coef_sub_surv<-c(coef_sub_surv, summary(coxph(Surv(surv[,1], surv[,2])~MEs[,i]+meta$NOT_IN_OSLOVAL_Pam50Subtype))[[7]][1,1])
}

collapse<-cbind(p.adjust(p_sub, method="fdr"), p.adjust(p_sub_age, method="fdr"), p.adjust(p_sub_surv, method="fdr"))
rownames(collapse)<-colnames(MEs)
colnames(collapse)<-c("Grade", "Age at onset", "Survival")

coef_collapse<-cbind(coef_sub, coef_sub_age, coef_sub_surv)
rownames(coef_collapse)<-colnames(MEs)
colnames(coef_collapse)<-c("Grade", "Age at onset", "Survival")

coef_collapse<-coef_collapse[-which(rownames(coef_collapse)=="Unconnected"),]
collapse<-collapse[-which(rownames(collapse)=="Unconnected"),]

#collapse[collapse<2.2*10^(-16)]<-2.2*10^(-16)
#pheatmap(-log10(collapse))

png("ModuleGradeRelationship_sub.png", res = 300, width=1500, height = 2000)

moduleTraitCor_tot_grade<-coef_collapse
moduleTraitPvalue_tot_grade<-collapse

moduleTraitPvalue_tot_grade<-moduleTraitPvalue_tot_grade[order(moduleTraitCor_tot_grade[,1], decreasing = T),]
moduleTraitCor_tot_grade<-moduleTraitCor_tot_grade[order(moduleTraitCor_tot_grade[,1], decreasing = T),]

textMatrix = paste(signif(moduleTraitCor_tot_grade, 2), "\n(",
                   signif(moduleTraitPvalue_tot_grade, 1), ")", sep = "");

textMatrix[moduleTraitPvalue_tot_grade<2.2*10^(-16)]<-paste(signif(moduleTraitCor_tot_grade[moduleTraitPvalue_tot_grade<2.2*10^(-16)], 2), "\n(",
                                                            
                                                            
                                                            "<2.2e-16)", sep = "");

dim(textMatrix) = dim(moduleTraitCor_tot_grade)
par(mar = c(6, 14, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor_tot_grade,
               xLabels =c("Grade","Age at onset", "Survival"),
               
               yLabels = rownames(moduleTraitPvalue_tot_grade),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.x=1,
               cex.lab.y=0.8,
               zlim = c(min(moduleTraitCor_tot_grade[,1]),max(moduleTraitCor_tot_grade[,1])),
               main = paste("Module - grade"))
dev.off();

