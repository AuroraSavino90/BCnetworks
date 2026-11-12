set.seed(123)
####setting paths
pathdata <- "data"
pathscript <- "pipelines"
resultPath<-'results/20250221/'

###loading input data
gene_annot <- read.csv("data/CRISPR/gene_identifiers_20241212.csv")
### gene_identifiers_20191101 downloaded from https://cog.sanger.ac.uk/cmp/download/gene_identifiers_20241212.csv on 20250221

CMP_annot <- read.csv("data/CRISPR/model_list_20241120.csv")
### model_list_20240110.csv downloaded from https://cog.sanger.ac.uk/cmp/download/model_list_20241120.csv on 20250129

#select only breast carcinoma cell lines
BC_lines<-CMP_annot$model_name[CMP_annot$cancer_type=="Breast Carcinoma"]


scaled_depFC<-read.csv('data/CRISPR//CRISPRGeneEffect.csv', row.names = 1, check.names = FALSE)
colnames(scaled_depFC)<-sub(" \\(.*","",colnames(scaled_depFC)) # "GENE (EntrezID)" -> "GENE"; keeps hyphenated symbols (e.g. HLA-A)
###scaled essentiality matrices downloaded from https://depmap.org/portal/data_page/?tab=allData on 20250129 (24Q4)

toremove<-which(is.na(CMP_annot$model_id[match(rownames(scaled_depFC),CMP_annot$BROAD_ID)]))
scaled_depFC<-scaled_depFC[-toremove,]
rownames(scaled_depFC)<-CMP_annot$model_id[match(rownames(scaled_depFC),CMP_annot$BROAD_ID)]

scaled_depFC<-t(scaled_depFC)

CMP_annot$model_name[CMP_annot$cancer_type=="Breast Carcinoma"][which(CMP_annot$model_id[CMP_annot$cancer_type=="Breast Carcinoma"] %in% colnames(scaled_depFC))]

BL<-c("HCC1143", "Hs-578-T", "BT-549", "MDA-MB-231", "MFM-223",
      "MDA-MB-157", "HCC38", "HCC70", "MDA-MB-468", "HCC1937", "HCC1187", "HCC1806",
      "CAL-120", "COLO-824", "DU-4475", "SUM-185PE", "SUM-159PT", "SUM-149PT",
      "SUM-229PE", "HMC-1-8", "VP229", "SUM-1315MO2", "MDA-MB-436", "HCC1395", 
      "SUM-102PT", "CAL-51")
BL_id<-CMP_annot$model_id[CMP_annot$model_name %in% BL]

###pan-cancer reference (mean gene effect in non-breast cancer lines) and all breast cancer lines, for the specificity analyses
pan_effect<-rowMeans(scaled_depFC[, setdiff(colnames(scaled_depFC), CMP_annot$model_id[CMP_annot$cancer_type %in% c("Breast Carcinoma", "Non-Cancerous")])], na.rm=TRUE)
BC_id<-intersect(CMP_annot$model_id[CMP_annot$cancer_type=="Breast Carcinoma"], colnames(scaled_depFC))
stopifnot(all(BL_id %in% BC_id))
depFC_BC<-scaled_depFC[,BC_id]
depFC_BC<-depFC_BC[rowSums(is.na(depFC_BC))==0 & !is.na(pan_effect[rownames(depFC_BC)]),]

###get CRISPR data for BL lines
scaled_depFC<-scaled_depFC[,BL_id]

#remove genes with missing values
scaled_depFC<-scaled_depFC[-which(rowSums(is.na(scaled_depFC))>0),]



###corelation between essentiality and centrality for each module
load(file="results/2025/centrality_basal.RData")

##keep only genes in common between centrality_basal and essentiality matrices
centrality_basal<-centrality_basal[which(rownames(centrality_basal) %in% rownames(scaled_depFC)),]

pval<-matrix(nrow=length(unique(centrality_basal$module)), ncol=ncol(scaled_depFC))
cc<-matrix(nrow=length(unique(centrality_basal$module)), ncol=ncol(scaled_depFC))
for(j in 1:length(unique(centrality_basal$module))){
  for(i in 1:ncol(scaled_depFC)){
    pval[j,i]<-cor.test(centrality_basal$kWithin[centrality_basal$module==unique(centrality_basal$module)[j]], scaled_depFC[rownames(centrality_basal)[centrality_basal$module==unique(centrality_basal$module)[j]],i])[[3]]
    cc[j,i]<-cor.test(centrality_basal$kWithin[centrality_basal$module==unique(centrality_basal$module)[j]], scaled_depFC[rownames(centrality_basal)[centrality_basal$module==unique(centrality_basal$module)[j]],i])[[4]]
  }
}
rownames(cc)<-unique(centrality_basal$module)
colnames(cc)<-colnames(scaled_depFC)
rownames(pval)<-unique(centrality_basal$module)
colnames(pval)<-colnames(scaled_depFC)

#remove unconnected
pval<-pval[rownames(pval)!="b_Unconnected",]
cc<-cc[rownames(cc)!="b_Unconnected",]

#compute average correlation
cc_mean<-rowMeans(cc)

#adjust pvalue
adj_pval<-matrix(p.adjust(pval, method = "BH"), nrow=nrow(pval))
colnames(adj_pval)<-colnames(pval)
rownames(adj_pval)<-rownames(pval)
#cut pvalue at 2.2*10^-16
adj_pval[adj_pval<2.2*10^(-16)]<-2.2*10^(-16)

#one test per module: per-line tests share genes and kWithin and gene effects are highly correlated across lines,
#so their p-values are not independent and are not combined (Fisher would be anticonservative)
#(1) kWithin vs mean gene effect across BL lines; (2) the same, partial correlation adjusted for the pan-cancer gene effect;
#(3) basal vs other breast cancer lines: difference in partial correlation, p from 1000 permutations of line labels
pcor<-function(k, y, z) cor(resid(lm(k~z)), resid(lm(y~z)))
mod_test<-do.call(rbind, lapply(rownames(cc), function(m){
  g<-rownames(centrality_basal)[centrality_basal$module==m]
  k<-centrality_basal[g, "kWithin"]
  t1<-cor.test(k, rowMeans(scaled_depFC[g, ]))
  t2<-cor.test(resid(lm(k~pan_effect[g])), resid(lm(rowMeans(scaled_depFC[g, ])~pan_effect[g])))
  gb<-intersect(g, rownames(depFC_BC)); kb<-centrality_basal[gb, "kWithin"]; zb<-pan_effect[gb]
  st<-function(lb) pcor(kb, rowMeans(depFC_BC[gb, lb]), zb)-pcor(kb, rowMeans(depFC_BC[gb, setdiff(BC_id, lb)]), zb)
  obs<-st(BL_id)
  null<-replicate(1000, st(sample(BC_id, length(BL_id))))
  data.frame(module=m, n_genes=length(g), r=t1$estimate, p=t1$p.value,
             partial_r_pan_adjusted=t2$estimate, p_partial=t2$p.value,
             partial_r_basal=pcor(kb, rowMeans(depFC_BC[gb, BL_id]), zb), partial_r_otherBC=pcor(kb, rowMeans(depFC_BC[gb, setdiff(BC_id, BL_id)]), zb),
             delta_partial_r=obs, p_basal_more_negative=(1+sum(null<=obs))/1001, p_two_sided=(1+sum(abs(null)>=abs(obs)))/1001)
}))
mod_test$fdr<-p.adjust(mod_test$p, method="BH")
mod_test$fdr_partial<-p.adjust(mod_test$p_partial, method="BH")
mod_test$fdr_basal_more_negative<-p.adjust(mod_test$p_basal_more_negative, method="BH")
write.csv(mod_test, "results/2025_fix/CRISPR_centrality_module_tests.csv", row.names=FALSE)

pmerged<-setNames(mod_test$p, mod_test$module)

cc[adj_pval>0.05]<-NA
cc<-t(cc)
rownames(cc)<-CMP_annot$model_name[match(rownames(cc), CMP_annot$model_id)]


library(pheatmap)

paletteLength <- 50
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
lim<-max(abs(unlist(cc)), na.rm=T)
myBreaks <- seq(-lim, lim, length.out=paletteLength+1) # symmetric scale: significant positive correlations are shown too

anno_p<-data.frame(meanR= cc_mean, minus_log10_p= -log10(pmerged[names(cc_mean)]))

# Define a continuous color gradient
continuous_colors <- colorRampPalette(c("white", "mediumorchid4"))(100)
rng<-range(cc_mean); nneg<-round(100*max(0, -rng[1])/diff(rng)) # diverging palette with white at meanR = 0
meanR_colors<-c(colorRampPalette(c("#4575B4", "white"))(nneg), colorRampPalette(c("white", "#D73027"))(100-nneg))

# Map the continuous annotation to colors
annotation_colors <- list(meanR = meanR_colors, minus_log10_p = continuous_colors)


graphics.off()
png("results/2025_fix/CRISPR_corr.png", res=300, 3000, 3000)
pheatmap((cc[,names(sort(colSums(cc, na.rm=T)))]), cluster_cols = F, cluster_rows = F, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, annotation_col = anno_p, annotation_colors = annotation_colors)
dev.off()


graphics.off()
pdf("results/2025_fix/CRISPR_corr.pdf", 10,10)
pheatmap((cc[,names(sort(colSums(cc, na.rm=T)))]), cluster_cols = F, cluster_rows = F, cellwidth=15, cellheight=15, breaks=myBreaks, color = myColor, annotation_col = anno_p, annotation_colors = annotation_colors)
dev.off()


###b_E2F_TARGETS is the pre-specified module (basal-specific connectivity, section 4): basal vs other BC test repeated with 10000 permutations for a stable p-value, reported without correction across modules
gE<-intersect(rownames(centrality_basal)[centrality_basal$module=="b_E2F_TARGETS"], rownames(depFC_BC)); kE<-centrality_basal[gE, "kWithin"]; zE<-pan_effect[gE]
stE<-function(lb) pcor(kE, rowMeans(depFC_BC[gE, lb]), zE)-pcor(kE, rowMeans(depFC_BC[gE, setdiff(BC_id, lb)]), zE)
obsE<-stE(BL_id)
nullE<-replicate(10000, stE(sample(BC_id, length(BL_id))))
e2f<-data.frame(module="b_E2F_TARGETS", n_genes=length(gE), delta_partial_r=obsE, n_perm=10000, p_basal_more_negative=(1+sum(nullE<=obsE))/10001, p_two_sided=(1+sum(abs(nullE)>=abs(obsE)))/10001)
write.csv(e2f, "results/2025_fix/CRISPR_bE2F_basal_vs_otherBC_10000perm.csv", row.names=FALSE)

###Supplementary figure: (A) kWithin-essentiality correlation before/after pan-cancer adjustment; (B) pan-cancer adjusted correlation in basal vs other BC lines
library(ggplot2)
library(ggrepel)
library(ggpubr)
ord<-mod_test$module[order(mod_test$r)]
dfA<-rbind(data.frame(module=mod_test$module, r=mod_test$r, fdr=mod_test$fdr, type="unadjusted"),
           data.frame(module=mod_test$module, r=mod_test$partial_r_pan_adjusted, fdr=mod_test$fdr_partial, type="pan-cancer adjusted"))
dfA$module<-factor(dfA$module, levels=ord)
dfA$type<-factor(dfA$type, levels=c("unadjusted", "pan-cancer adjusted"))
pA<-ggplot(dfA, aes(x=r, y=module))+geom_vline(xintercept=0, colour="grey60")+
  geom_line(aes(group=module), colour="grey70")+geom_point(aes(colour=type), size=2.5)+
  geom_text(data=dfA[dfA$fdr<0.05,], aes(label="*", colour=type), nudge_y=0.25, size=5, show.legend=FALSE)+
  scale_colour_manual(values=c("unadjusted"="#D73027", "pan-cancer adjusted"="#4575B4"))+
  labs(x="Pearson r (kWithin vs mean gene effect, basal-like lines)", y=NULL, colour=NULL, caption="* FDR < 0.05")+
  theme_classic()+theme(legend.position="top")
dfB<-mod_test
dfB$group<-ifelse(dfB$p_basal_more_negative<0.05, "more negative in basal-like (p<0.05)", ifelse(dfB$p_two_sided<0.05 & dfB$delta_partial_r>0, "more negative in other BC (p<0.05)", "n.s."))
dfB$label<-ifelse(dfB$group!="n.s.", dfB$module, "")
dfB$label[dfB$module=="b_E2F_TARGETS"]<-sprintf("b_E2F_TARGETS\np = %.3f (10000 perm.)", e2f$p_basal_more_negative)
lims<-range(c(dfB$partial_r_basal, dfB$partial_r_otherBC))
pB<-ggplot(dfB, aes(x=partial_r_otherBC, y=partial_r_basal))+geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey60")+
  geom_hline(yintercept=0, colour="grey85")+geom_vline(xintercept=0, colour="grey85")+
  geom_point(aes(colour=group), size=2.5)+geom_text_repel(aes(label=label), size=3, min.segment.length=0, box.padding=0.5, max.overlaps=Inf)+
  scale_colour_manual(values=c("more negative in basal-like (p<0.05)"="#D73027", "more negative in other BC (p<0.05)"="#4575B4", "n.s."="grey60"))+
  coord_equal(xlim=lims, ylim=lims)+
  labs(x="Partial r, other BC lines (n=23)", y="Partial r, basal-like lines (n=26)", colour=NULL)+
  theme_classic()+theme(legend.position="top", legend.direction="vertical")
pS<-ggarrange(pA, pB, labels=c("A", "B"), ncol=2, widths=c(1, 1))

png("results/2025_fix/CRISPR_pan_adjusted_specificity.png", res=300, 4200, 2200)
print(pS)
dev.off()

pdf("results/2025_fix/CRISPR_pan_adjusted_specificity.pdf", 14, 7.3)
print(pS)
dev.off()
