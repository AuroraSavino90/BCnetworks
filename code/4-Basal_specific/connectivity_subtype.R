library(WGCNA)
library(ggpubr)
load(file="data/RData/metabric.RData")
load(file="data/RData/meta.RData")
load(file="results/2025/centrality_basal.RData")

###moduleColors_basal in other subtypes
moduleColors_basal<-centrality_basal$module

## sample-label permutation test for kWithin (Basal vs each other subtype).
## Group sizes are kept in the permutations, so the null includes the effect of sample size on kWithin.
perm_kwithin<-function(X, subtype, mods, n_perm=1000){
  res<-NULL
  for(o in c("LumA", "LumB", "Her2")){
    idx<-which(subtype %in% c("Basal", o))
    for(m in setdiff(unique(mods), "b_Unconnected")){
      g<-names(mods)[which(mods==m)]
      mk<-function(s){A<-((1+cor(t(X[g, idx[s]]), use="p"))/2)^6; diag(A)<-0; mean(rowSums(A, na.rm=TRUE))} # mean kWithin, signed adjacency, power 6; constant genes (NA cor) count 0, as in intramodularConnectivity
      stat<-function(l) mk(l=="Basal")-mk(l==o)
      obs<-stat(subtype[idx])
      null<-replicate(n_perm, stat(sample(subtype[idx])))
      res<-rbind(res, data.frame(module=m, comparison=paste0("Basal_vs_", o), delta_kwithin=obs, null_mean_delta=mean(null), p_perm=(1+sum(null>=obs))/(1+n_perm)))
    }
  }
  res$fdr<-ave(res$p_perm, res$comparison, FUN=function(p) p.adjust(p, method="BH"))
  res
}
set.seed(1)

ADJ1=adjacency(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="LumA"]), power = 6, type= "signed")
Alldegrees_basalmod_lumA=intramodularConnectivity(ADJ1, moduleColors_basal)


ADJ1=adjacency(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="LumB"]), power = 6, type= "signed")
Alldegrees_basalmod_lumB=intramodularConnectivity(ADJ1, moduleColors_basal)


ADJ1=adjacency(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Her2"]), power = 6, type= "signed")
Alldegrees_basalmod_her2=intramodularConnectivity(ADJ1, moduleColors_basal)

ADJ1=adjacency(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]), power = 6, type= "signed")
Alldegrees_basalmod_basal=intramodularConnectivity(ADJ1, moduleColors_basal)

library(ggpubr)
for(i in 1:length(unique(moduleColors_basal))){
  
  df<-data.frame(kwithin=c(Alldegrees_basalmod_basal[which(moduleColors_basal==unique(moduleColors_basal)[i]),2], 
                           Alldegrees_basalmod_lumA[which(moduleColors_basal==unique(moduleColors_basal)[i]),2], 
                           Alldegrees_basalmod_lumB[which(moduleColors_basal==unique(moduleColors_basal)[i]),2], 
                           Alldegrees_basalmod_her2[which(moduleColors_basal==unique(moduleColors_basal)[i]),2]),
                 subtype=rep(c("Basal", "LumA", "LumB", "Her2"), each=length(which(moduleColors_basal==unique(moduleColors_basal)[i]))))
  
  my_comparisons <- list( c("Basal", "LumA"), c("Basal", "LumB"), c("Basal", "Her2") )
  p <- ggboxplot(df, x = "subtype", y = "kwithin")
  #  Add p-value
  p <- p + stat_compare_means(comparisons = my_comparisons, ref.group = "Basal", 
                              method = "wilcox.test", method.args = list(alternative = "greater", paired = TRUE)) # Add pairwise comparisons p-value
  
  
  png(paste("results/2025_fix/",unique(moduleColors_basal)[i],"moduleColorsbasal_inothersubtypes_kwithin.png",sep=""), res=300, 1000, 1000)
 print(p)
 dev.off()
 
 pdf(paste("results/2025_fix/",unique(moduleColors_basal)[i],"moduleColorsbasal_inothersubtypes_kwithin.pdf",sep=""), 5,5)
 print(p)
 dev.off()
}
perm_metabric<-perm_kwithin(metabric, meta$NOT_IN_OSLOVAL_Pam50Subtype, setNames(moduleColors_basal, rownames(metabric)))
write.csv(perm_metabric, "results/2025_fix/kwithin_permutation_METABRIC.csv", row.names=FALSE)


####module eigengenes across subtypes

MEs_basal= moduleEigengenes(t(metabric), centrality_basal$module)$eigengenes
colnames(MEs_basal)<-gsub("^ME","", colnames(MEs_basal))


for(i in 1:length(unique(moduleColors_basal))){
  
  df<-data.frame(ME=c(MEs_basal[meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal",unique(moduleColors_basal)[i]], 
                      MEs_basal[meta$NOT_IN_OSLOVAL_Pam50Subtype=="LumA",unique(moduleColors_basal)[i]], 
                      MEs_basal[meta$NOT_IN_OSLOVAL_Pam50Subtype=="LumB",unique(moduleColors_basal)[i]], 
                      MEs_basal[meta$NOT_IN_OSLOVAL_Pam50Subtype=="Her2",unique(moduleColors_basal)[i]]),
                 subtype=rep(c("Basal", "LumA", "LumB", "Her2"), times=table(meta$NOT_IN_OSLOVAL_Pam50Subtype)[c("Basal", "LumA", "LumB", "Her2")]))
  
  my_comparisons <- list( c("Basal", "LumA"), c("Basal", "LumB"), c("Basal", "Her2") )
  p <- ggboxplot(df, x = "subtype", y = "ME")
  #  Add p-value
  p <- p + stat_compare_means(comparisons = my_comparisons, ref.group = "Basal", 
                              method = "wilcox.test", method.args = list(alternative = "two.sided")) # Add pairwise comparisons p-value
  
  
  png(paste("results/2025_fix/",unique(moduleColors_basal)[i],"moduleColorsbasal_inothersubtypes_ME.png",sep=""), res=300, 1000, 1000)
  print(p)
  dev.off()
  pdf(paste("results/2025_fix/",unique(moduleColors_basal)[i],"moduleColorsbasal_inothersubtypes_ME.pdf",sep=""), 5,5)
  print(p)
  dev.off()
}

#######################
### in TCGA
########################
load("../../../Data/Data_downloaded/TCGA/TCGA-BRCA_Primary Tumor_dedupl.RData")
load("../../../Data/Data_downloaded/TCGA/TCGA-BRCA_clinical.RData")

library(TCGAbiolinks)
BRCA_subtype <- TCGAquery_subtype(tumor = "BRCA")

anno<-data.frame(barcode=colnames(TCGA_dedupl), 
                 subtype=BRCA_subtype$BRCA_Subtype_PAM50[match(substr(colnames(TCGA_dedupl),1,12), BRCA_subtype$patient)])

incommon<-intersect(rownames(metabric), rownames(TCGA_dedupl))
TCGA_dedupl_ic<-TCGA_dedupl[incommon,]
names(moduleColors_basal)<-rownames(metabric)
moduleColors_basal_ic<-moduleColors_basal[incommon]


ADJ1=adjacency(t(TCGA_dedupl_ic[,which(anno$subtype=="LumA")]), power = 6, type= "signed")
Alldegrees_basalmod_lumA=intramodularConnectivity(ADJ1, moduleColors_basal_ic)


ADJ1=adjacency(t(TCGA_dedupl_ic[,which(anno$subtype=="LumB")]), power = 6, type= "signed")
Alldegrees_basalmod_lumB=intramodularConnectivity(ADJ1, moduleColors_basal_ic)


ADJ1=adjacency(t(TCGA_dedupl_ic[,which(anno$subtype=="Her2")]), power = 6, type= "signed")
Alldegrees_basalmod_her2=intramodularConnectivity(ADJ1, moduleColors_basal_ic)

ADJ1=adjacency(t(TCGA_dedupl_ic[,which(anno$subtype=="Basal")]), power = 6, type= "signed")
Alldegrees_basalmod_basal=intramodularConnectivity(ADJ1, moduleColors_basal_ic)


for(i in 1:length(unique(moduleColors_basal_ic))){
  
  df<-data.frame(kwithin=c(Alldegrees_basalmod_basal[which(moduleColors_basal_ic==unique(moduleColors_basal_ic)[i]),2], 
                           Alldegrees_basalmod_lumA[which(moduleColors_basal_ic==unique(moduleColors_basal_ic)[i]),2], 
                           Alldegrees_basalmod_lumB[which(moduleColors_basal_ic==unique(moduleColors_basal_ic)[i]),2], 
                           Alldegrees_basalmod_her2[which(moduleColors_basal_ic==unique(moduleColors_basal_ic)[i]),2]),
                 subtype=rep(c("Basal", "LumA", "LumB", "Her2"), each=length(which(moduleColors_basal_ic==unique(moduleColors_basal_ic)[i]))))
  
  my_comparisons <- list( c("Basal", "LumA"), c("Basal", "LumB"), c("Basal", "Her2") )
  p <- ggboxplot(df, x = "subtype", y = "kwithin")
  #  Add p-value
  p <- p + stat_compare_means(comparisons = my_comparisons, ref.group = "Basal", 
                              method = "wilcox.test", method.args = list(alternative = "greater", paired = TRUE)) # Add pairwise comparisons p-value
  
  
  png(paste("results/2025_fix/",unique(moduleColors_basal_ic)[i],"moduleColorsbasal_inothersubtypes_kwithin_TCGA.png",sep=""), res=300, 1000, 1000)
  print(p)
  dev.off()
  pdf(paste("results/2025_fix/",unique(moduleColors_basal_ic)[i],"moduleColorsbasal_inothersubtypes_kwithin_TCGA.pdf",sep=""), 5,5)
  print(p)
  dev.off()
}
perm_tcga<-perm_kwithin(TCGA_dedupl_ic, anno$subtype, moduleColors_basal_ic)
write.csv(perm_tcga, "results/2025_fix/kwithin_permutation_TCGA.csv", row.names=FALSE)

####module eigengenes across subtypes

metabric_ic<-metabric[incommon,]
centrality_basal_ic<-centrality_basal[incommon,]

#####################################################################
###function to project MEs computed on a dataset (original_data) on new transcriptomic data (newdata)
#####################################################################

pcaproject <- function(
    newdata,
    original_data,
    modules,
    ME = NULL,
    min_common_genes = 2L,
    min_overlap_fraction = 0.50
) {
  original_data <- as.matrix(original_data)
  newdata <- as.matrix(newdata)
  
  if (is.null(rownames(original_data)) || is.null(rownames(newdata))) {
    stop("Both original_data and newdata must have gene symbols as row names.")
  }
  
  if (length(modules) != nrow(original_data)) {
    stop("'modules' must contain one module label for every row of original_data.")
  }
  
  # Align module labels explicitly to original_data gene order.
  if (!is.null(names(modules))) {
    missing_module_genes <- setdiff(rownames(original_data), names(modules))
    if (length(missing_module_genes) > 0) {
      stop(
        "Missing module assignments for ",
        length(missing_module_genes),
        " genes in original_data."
      )
    }
    modules <- modules[rownames(original_data)]
  }
  
  modules <- as.character(modules)
  module_names <- levels(factor(modules[!is.na(modules)]))
  
  pcaproj <- matrix(
    NA_real_,
    nrow = ncol(newdata),
    ncol = length(module_names),
    dimnames = list(colnames(newdata), module_names)
  )
  
  diagnostics <- vector("list", length(module_names))
  names(diagnostics) <- module_names
  
  for (module_name in module_names) {
    module_genes <- rownames(original_data)[modules == module_name]
    common_genes <- intersect(module_genes, rownames(newdata))
    
    n_module_genes <- length(module_genes)
    overlap_fraction <- length(common_genes) / n_module_genes
    
    if (length(common_genes) < min_common_genes) {
      warning(
        "Module ", module_name,
        " has fewer than ", min_common_genes,
        " common genes; projection set to NA."
      )
      diagnostics[[module_name]] <- data.frame(
        module = module_name,
        n_module_genes = n_module_genes,
        n_common_genes = length(common_genes),
        overlap_fraction = overlap_fraction,
        n_genes_used = 0L,
        stringsAsFactors = FALSE
      )
      next
    }
    
    if (overlap_fraction < min_overlap_fraction) {
      warning(
        "Module ", module_name,
        " has only ", round(100 * overlap_fraction, 1),
        "% gene overlap between original_data and newdata."
      )
    }
    
    # Refit the reference PC on the genes that can actually be projected.
    # This avoids treating absent genes as zero contributions.
    training <- t(original_data[common_genes, , drop = FALSE])
    
    gene_sd <- apply(training, 2, sd, na.rm = TRUE)
    valid_genes <- colnames(training)[
      is.finite(gene_sd) & gene_sd > 0 &
        colSums(!is.finite(training)) == 0
    ]
    
    if (length(valid_genes) < min_common_genes) {
      warning(
        "Module ", module_name,
        " has fewer than ", min_common_genes,
        " finite, variable common genes; projection set to NA."
      )
      diagnostics[[module_name]] <- data.frame(
        module = module_name,
        n_module_genes = n_module_genes,
        n_common_genes = length(common_genes),
        overlap_fraction = overlap_fraction,
        n_genes_used = length(valid_genes),
        stringsAsFactors = FALSE
      )
      next
    }
    
    training <- training[, valid_genes, drop = FALSE]
    
    # WGCNA::moduleEigengenes uses scale = TRUE by default. Therefore the
    # reference PCA and the new-data projection must use the same gene-wise
    # centering and scaling estimated in METABRIC.
    pca <- prcomp(
      training,
      center = TRUE,
      scale. = TRUE,
      rank. = 1
    )
    
    # prcomp$x[,1] = U[,1] * d1, whereas moduleEigengenes returns U[,1].
    singular_value_1 <- pca$sdev[1] * sqrt(nrow(training) - 1)
    
    if (!is.finite(singular_value_1) || singular_value_1 <= 0) {
      warning("Invalid first singular value for module ", module_name, ".")
      next
    }
    
    # Orient the projected PC exactly as the corresponding WGCNA eigengene.
    if (!is.null(ME) && module_name %in% colnames(ME)) {
      common_samples <- intersect(rownames(pca$x), rownames(ME))
      
      if (length(common_samples) >= 3) {
        reference_score <- pca$x[common_samples, 1] / singular_value_1
        orientation_cor <- suppressWarnings(
          cor(
            reference_score,
            ME[common_samples, module_name],
            use = "pairwise.complete.obs"
          )
        )
        
        if (is.finite(orientation_cor) && orientation_cor < 0) {
          pca$rotation[, 1] <- -pca$rotation[, 1]
        }
      }
    } else {
      # Fallback matching WGCNA's default align = "along average".
      reference_score <- pca$x[, 1] / singular_value_1
      average_expression <- rowMeans(scale(training), na.rm = TRUE)
      orientation_cor <- suppressWarnings(
        cor(reference_score, average_expression, use = "pairwise.complete.obs")
      )
      
      if (is.finite(orientation_cor) && orientation_cor < 0) {
        pca$rotation[, 1] <- -pca$rotation[, 1]
      }
    }
    
    new_samples_by_genes <- t(newdata[valid_genes, , drop = FALSE])
    
    if (any(!is.finite(new_samples_by_genes))) {
      stop(
        "Non-finite values found in newdata for module ", module_name,
        ". Clean or impute TPMlog before projection."
      )
    }
    
    new_scaled <- sweep(
      new_samples_by_genes,
      2,
      pca$center,
      FUN = "-"
    )
    new_scaled <- sweep(
      new_scaled,
      2,
      pca$scale,
      FUN = "/"
    )
    
    pcaproj[, module_name] <- as.numeric(
      new_scaled %*% pca$rotation[, 1, drop = FALSE]
    ) / singular_value_1
    
    diagnostics[[module_name]] <- data.frame(
      module = module_name,
      n_module_genes = n_module_genes,
      n_common_genes = length(common_genes),
      overlap_fraction = overlap_fraction,
      n_genes_used = length(valid_genes),
      stringsAsFactors = FALSE
    )
  }
  
  attr(pcaproj, "projection_diagnostics") <- do.call(rbind, diagnostics)
  pcaproj
}

validate_pcaproject <- function(
    original_data,
    modules,
    ME,
    minimum_abs_correlation = 0.999
) {
  projected_reference <- pcaproject(
    newdata = original_data,
    original_data = original_data,
    modules = modules,
    ME = ME,
    min_overlap_fraction = 1
  )
  
  common_modules <- intersect(
    colnames(projected_reference),
    colnames(ME)
  )
  
  validation <- data.frame(
    module = common_modules,
    correlation_with_WGCNA_ME = vapply(
      common_modules,
      function(module_name) {
        suppressWarnings(
          cor(
            projected_reference[, module_name],
            ME[, module_name],
            use = "pairwise.complete.obs"
          )
        )
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
  
  validation$abs_correlation <- abs(
    validation$correlation_with_WGCNA_ME
  )
  
  failed <- validation[
    !is.finite(validation$abs_correlation) |
      validation$abs_correlation < minimum_abs_correlation,
    ,
    drop = FALSE
  ]
  
  if (nrow(failed) > 0) {
    print(failed)
    stop(
      "The PCA projection does not reproduce the reference WGCNA eigengenes ",
      "for one or more modules."
    )
  }
  
  validation
}

rownames(MEs_basal)<-colnames(metabric)
pcaproj<-pcaproject(newdata=TCGA_dedupl_ic, original_data=metabric, modules=centrality_basal$module, ME=MEs_basal)[,unique(centrality_basal_ic$module)]

colnames(pcaproj)<-unique(centrality_basal_ic$module)


for(i in 1:length(unique(moduleColors_basal_ic))){
  
  df<-data.frame(ME=c(pcaproj[which(anno$subtype=="Basal"),which(unique(centrality_basal_ic$module)==unique(moduleColors_basal_ic)[i])], 
                           pcaproj[which(anno$subtype=="LumA"),which(unique(centrality_basal_ic$module)==unique(moduleColors_basal_ic)[i])], 
                           pcaproj[which(anno$subtype=="LumB"),which(unique(centrality_basal_ic$module)==unique(moduleColors_basal_ic)[i])], 
                           pcaproj[which(anno$subtype=="Her2"),which(unique(centrality_basal_ic$module)==unique(moduleColors_basal_ic)[i])]),
                 subtype=rep(c("Basal", "LumA", "LumB", "Her2"), times=table(anno$subtype)[c("Basal", "LumA", "LumB", "Her2")]))
  
  my_comparisons <- list( c("Basal", "LumA"), c("Basal", "LumB"), c("Basal", "Her2") )
  p <- ggboxplot(df, x = "subtype", y = "ME")
  #  Add p-value
  p <- p + stat_compare_means(comparisons = my_comparisons, ref.group = "Basal", 
                              method = "wilcox.test", method.args = list(alternative = "two.sided")) # Add pairwise comparisons p-value
  
  
  png(paste("results/2025_fix/",unique(moduleColors_basal_ic)[i],"moduleColorsbasal_inothersubtypes_ME_TCGA.png",sep=""), res=300, 1000, 1000)
  print(p)
  dev.off()
  pdf(paste("results/2025_fix/",unique(moduleColors_basal_ic)[i],"moduleColorsbasal_inothersubtypes_ME_TCGA.pdf",sep=""), 5,5)
  print(p)
  dev.off()
}

#####################################################################
### basal-specific modules: METABRIC and TCGA combined (Fisher and Stouffer), Basal vs each other subtype
### one-sided (Basal higher) p-values, BH across modules x comparisons; pass = FDR<0.05 in all three comparisons
#####################################################################
library(metap)
comb_p<-function(d, p1, p2){
  p<-pmin(pmax(as.matrix(d[,c(p1, p2)]), 1e-300), 1-1e-10)
  d$p_fisher<-apply(p, 1, function(x) sumlog(x)$p)
  d$p_stouffer<-apply(p, 1, function(x) sumz(x)$p)
  d$fdr_fisher<-ave(d$p_fisher, d$comparison, FUN=function(p) p.adjust(p, method="BH")) # BH across modules, separately per comparison (as in section 3)
  d$fdr_stouffer<-ave(d$p_stouffer, d$comparison, FUN=function(p) p.adjust(p, method="BH"))
  d
}
# kWithin: permutation p-values are already one-sided
kw<-comb_p(merge(perm_metabric, perm_tcga, by=c("module", "comparison"), suffixes=c("_METABRIC", "_TCGA")), "p_perm_METABRIC", "p_perm_TCGA")
write.csv(kw, "results/2025_fix/kwithin_combined_METABRIC_TCGA.csv", row.names=FALSE)

# ME: two-sided Wilcoxon on samples, converted to one-sided (Basal higher) with two2one
me_test<-function(ME, subtype){
  res<-NULL
  for(o in c("LumA", "LumB", "Her2")){
    for(m in setdiff(colnames(ME), "b_Unconnected")){
      x<-ME[which(subtype=="Basal"), m]; y<-ME[which(subtype==o), m]
      w<-wilcox.test(x, y, alternative="two.sided")
      res<-rbind(res, data.frame(module=m, comparison=paste0("Basal_vs_", o), delta_median=median(x)-median(y),
                                 p_two_sided=w$p.value, p_basal_higher=two2one(w$p.value, invert=w$statistic<length(x)*length(y)/2)))
    }
  }
  res
}
stopifnot(identical(rownames(MEs_basal), rownames(meta)), identical(rownames(pcaproj), anno$barcode))
me<-comb_p(merge(me_test(MEs_basal, meta$NOT_IN_OSLOVAL_Pam50Subtype), me_test(pcaproj, anno$subtype), by=c("module", "comparison"), suffixes=c("_METABRIC", "_TCGA")),
           "p_basal_higher_METABRIC", "p_basal_higher_TCGA")
write.csv(me, "results/2025_fix/ME_basal_vs_subtypes_combined_METABRIC_TCGA.csv", row.names=FALSE)

sel<-merge(aggregate(cbind(kwithin_higher_stouffer=fdr_stouffer<0.05, kwithin_higher_fisher=fdr_fisher<0.05,
                           kwithin_higher_each_dataset=delta_kwithin_METABRIC>0 & fdr_METABRIC<0.05 & delta_kwithin_TCGA>0 & fdr_TCGA<0.05)~module, kw, all),
           aggregate(cbind(ME_higher_stouffer=fdr_stouffer<0.05, ME_higher_fisher=fdr_fisher<0.05)~module, me, all), by="module")
write.csv(sel, "results/2025_fix/basal_specific_modules.csv", row.names=FALSE)
