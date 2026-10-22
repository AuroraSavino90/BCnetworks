library(breastCancerNKI)
library(Biobase)
library(WGCNA)
library(genefu)
library(TCGAbiolinks)
library(pheatmap)
library(survcomp)
data(scmod2.robust)
data(pam50.robust)
# pcaproject() / validate_pcaproject(): copied verbatim from 2_DEGs_sh_final_with_enrichment_tables.R
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

# collapse probes to one row per gene symbol (probe with highest mean expression)
collapse_probes <- function(expr, symbols) {
  expr <- as.matrix(expr)
  symbols <- trimws(as.character(symbols))
  if (length(symbols) != nrow(expr)) {
    stop("'symbols' must have one entry per row of 'expr'.")
  }
  
  keep <- !is.na(symbols) & symbols != "" & symbols != "---" &
    !grepl("///", symbols, fixed = TRUE)
  
  probe_mean <- rowMeans(expr, na.rm = TRUE)
  probe_mean[!is.finite(probe_mean)] <- -Inf
  
  idx <- which(keep)
  ord <- idx[order(symbols[idx], -probe_mean[idx], seq_along(idx))]
  chosen <- ord[!duplicated(symbols[ord])]
  
  out <- expr[chosen, , drop = FALSE]
  rownames(out) <- symbols[chosen]
  
  attr(out, "collapse_audit") <- data.frame(
    n_probes = nrow(expr),
    n_probes_without_unique_symbol = sum(!keep),
    n_genes = nrow(out),
    n_genes_with_multiple_probes = sum(table(symbols[idx]) > 1),
    rule = "max_mean_probe_per_symbol",
    stringsAsFactors = FALSE
  )
  out
}

#############load network objects and data
dir.create("results/2025_fix", showWarnings=FALSE)

load("data/RData/metabric.RData")
load("data/RData/meta.RData")
load("data/RData/net_metabric_Basal.RData")
load("results/2025/centrality_basal.RData")
load("results/2025/centrality_global.RData")
load("data/Complete_METABRIC_Clinical_Survival_Data__DSS.rbin")

library("breastCancerMAINZ")
library("breastCancerTRANSBIG")
library("breastCancerUPP")
library("breastCancerUNT")
library("breastCancerNKI")
library("breastCancerVDX")


dn <- c("transbig", "unt", "upp", "mainz", "nki", "vdx")
dn.platform <- c("affy", "affy", "affy", "affy", "agilent", "affy")
feat<-c("grade", "age")

#############compute module eigengenes

MEs= moduleEigengenes(t(metabric), centrality_global$module)$eigengenes
colnames(MEs)<-gsub("^ME","", colnames(MEs))
MEs_basal= moduleEigengenes(t(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]), centrality_basal$module)$eigengenes
colnames(MEs_basal)<-gsub("^ME","", colnames(MEs_basal))
rownames(MEs)<-colnames(metabric)
rownames(MEs_basal)<-colnames(metabric)[meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]
validate_pcaproject(metabric, centrality_global$module, MEs)
validate_pcaproject(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], centrality_basal$module, MEs_basal)

#######################
### FUNCTION: module vs clinical feature
######################

moduletrait<-function(features, metadata, original_dataset=metabric, new_dataset,genes_annot, centrality, MEs, new_dataset_name="TRANSBIG", basal=F){
  
  probes<-new_dataset
  new_dataset<-collapse_probes(new_dataset, rownames(new_dataset))
  if(basal){
    PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=t(probes),
                                    annot=genes_annot)
    PAM50<-PAM50Preds$subtype
    new_dataset<-new_dataset[,PAM50=="Basal"]
    metadata<-metadata[PAM50=="Basal",]
  }
  
  new_dataset<-new_dataset[which(rowSums(is.na(new_dataset))==0),]
  
  # MEs approximation for new_dataset
  # perform principal components analysis
  # project new data onto the PCA space
  
  pcaproj<-pcaproject(newdata=new_dataset, original_data=original_dataset, modules=centrality$module, ME=MEs)[,unique(centrality$module)]
  
  colnames(pcaproj)<-unique(centrality$module)
  
  
  ##Correlation between the projected MEs and clinical features
  
  moduleTraitCor = list();
  moduleTraitPvalue = list();
  # Calculate the correlations
  for (set in features)
  {
    moduleTraitCor[[set]] = cor(pcaproj, as.numeric(as.character(metadata[,set])), use = "p");
    moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], sum(!is.na(as.numeric(as.character(metadata[,set])))));
  }
  
  moduleTraitCor_tot<-moduleTraitCor[[1]]
  moduleTraitPvalue_tot<-moduleTraitPvalue[[1]]
  
  for (set in 2:length(features)){
    moduleTraitCor_tot<-cbind(moduleTraitCor_tot,moduleTraitCor[[set]])
    moduleTraitPvalue_tot<-cbind(moduleTraitPvalue_tot, moduleTraitPvalue[[set]])
  }
  colnames(moduleTraitCor_tot)<-features
  colnames(moduleTraitPvalue_tot)<-features
  
  if(basal){
    save(moduleTraitCor_tot, file=paste("results/2025_fix/moduleTraitCor_",new_dataset_name,"_basal.RData", sep=""))
    save(moduleTraitPvalue_tot, file=paste("results/2025_fix/moduleTraitPvalue_", new_dataset_name, "_basal.RData", sep=""))
    
  } else {
    save(moduleTraitCor_tot, file=paste("results/2025_fix/moduleTraitCor_",new_dataset_name,".RData", sep=""))
    save(moduleTraitPvalue_tot, file=paste("results/2025_fix/moduleTraitPvalue_", new_dataset_name, ".RData", sep=""))
  }
}

####################################
############ FUNCTION: MODULE vs SURVIVAL
#####################################

modulesurv<-function(features, surv, original_dataset=metabric, new_dataset,genes_annot, centrality, MEs, basal=F){
  
  probes<-new_dataset
  new_dataset<-collapse_probes(new_dataset, rownames(new_dataset))
  if(basal){
    PAM50Preds<-molecular.subtyping(sbt.model = "pam50",data=t(probes),
                                    annot=genes_annot)
    PAM50<-PAM50Preds$subtype
    new_dataset<-new_dataset[,PAM50=="Basal"]
    metadata<-metadata[PAM50=="Basal",]
    surv<-surv[PAM50=="Basal",]
  }
  
  new_dataset<-new_dataset[which(rowSums(is.na(new_dataset))==0),]
  
  # MEs approximation for new_dataset
  # perform principal components analysis
  # project new data onto the PCA space
  
  pcaproj<-pcaproject(newdata=new_dataset, original_data=original_dataset, modules=centrality$module, ME=MEs)[,unique(centrality$module)]
  
  colnames(pcaproj)<-unique(centrality$module)
  
  ################# relationship between projected ME and survival
  
  
  surv[which(surv[,1]>1825),"status"]<-0
  surv[which(surv[,1]>1825),"time"]<-1825
  
  pvalue_good<-c()
  pvalue_poor<-c()
  
  for(module in unique(centrality$module)){
    a=data.frame(ID=colnames(new_dataset),OS=surv$time, IND=surv$status)
    b=data.frame(eigen_gene=pcaproj[,module])
    c=cbind(a,b)
    eigen_class<-cut(pcaproj[,module], breaks=quantile(pcaproj[,module],probs=c(0,0.5,1), na.rm=T))
    
    
    pp<-survdiff(Surv(OS,IND)~eigen_class, data=c)
    pval<-1-pchisq(pp$chisq, length(pp$n)-1)
    
    if(pp$obs[1]>pp$exp[1]){
      pvalue_good<-c(pvalue_good, pval)
      pvalue_poor<-c(pvalue_poor, NA)
    } else {
      pvalue_good<-c(pvalue_good, NA)
      pvalue_poor<-c(pvalue_poor, pval)
    }
  }
  
  names(pvalue_poor)<-unique(centrality$module)
  names(pvalue_good)<-unique(centrality$module)
  
  return(list(pvalue_poor, pvalue_good))
  
}

#################################
####### All datasets
#################################

pvalue_poor_tot<-matrix(nrow=length(unique(centrality_global$module)), ncol=8)
pvalue_good_tot<-matrix(nrow=length(unique(centrality_global$module)), ncol=8)
pvalue_poor_tot_b<-matrix(nrow=length(unique(centrality_basal$module)), ncol=8)
pvalue_good_tot_b<-matrix(nrow=length(unique(centrality_basal$module)), ncol=8)

#######Correlation with grade

i<-1
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
data <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(data)<-dannot[rownames(data),3]

metadata<-pData(dd)
surv<-data.frame(status=metadata$e.rfs, time=metadata$t.rfs)

moduletrait(features=feat, metadata=metadata, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
            centrality=centrality_global, MEs=MEs, new_dataset_name="TRANSBIG", basal=F)
moduletrait(features=feat, metadata=metadata, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
            centrality=centrality_basal, MEs=MEs_basal, new_dataset_name="TRANSBIG", basal=T)


p<-modulesurv(features=feat, surv=surv, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
              centrality=centrality_global, MEs=MEs, basal=F)

pvalue_poor_tot[,2]<-p[[1]]
pvalue_good_tot[,2]<-p[[2]]

p<-modulesurv(features=feat, surv=surv, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
              centrality=centrality_basal, MEs=MEs_basal, basal=T)

pvalue_poor_tot_b[,2]<-p[[1]]
pvalue_good_tot_b[,2]<-p[[2]]


i<-2
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
data <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(data)<-dannot[rownames(data),3]

metadata<-pData(dd)
surv<-data.frame(status=metadata$e.rfs, time=metadata$t.rfs)


moduletrait(features=feat, metadata=metadata, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
            centrality=centrality_global, MEs=MEs, new_dataset_name="UNT", basal=F)
moduletrait(features=feat, metadata=metadata, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
            centrality=centrality_basal, MEs=MEs_basal, new_dataset_name="UNT", basal=T)
p<-modulesurv(features=feat, surv=surv, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
              centrality=centrality_global, MEs=MEs, basal=F)

pvalue_poor_tot[,3]<-p[[1]]
pvalue_good_tot[,3]<-p[[2]]

p<-modulesurv(features=feat, surv=surv, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
              centrality=centrality_basal, MEs=MEs_basal, basal=T)
pvalue_poor_tot_b[,3]<-p[[1]]
pvalue_good_tot_b[,3]<-p[[2]]

i<-3
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
data <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(data)<-dannot[rownames(data),3]

metadata<-pData(dd)
surv<-data.frame(status=metadata$e.rfs, time=metadata$t.rfs)


moduletrait(features=feat, metadata=metadata, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
            centrality=centrality_global, MEs=MEs, new_dataset_name="UPP", basal=F)
moduletrait(features=feat, metadata=metadata, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
            centrality=centrality_basal, MEs=MEs_basal, new_dataset_name="UPP", basal=T)
p<-modulesurv(features=feat, surv=surv, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
              centrality=centrality_global, MEs=MEs, basal=F)

pvalue_poor_tot[,4]<-p[[1]]
pvalue_good_tot[,4]<-p[[2]]

p<-modulesurv(features=feat, surv=surv, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
              centrality=centrality_basal, MEs=MEs_basal, basal=T)
pvalue_poor_tot_b[,4]<-p[[1]]
pvalue_good_tot_b[,4]<-p[[2]]


i<-4
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
data <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(data)<-dannot[rownames(data),3]

metadata<-pData(dd)
surv<-data.frame(status=metadata$e.dmfs, time=metadata$t.dmfs)

moduletrait(features=feat, metadata=metadata, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
            centrality=centrality_global, MEs=MEs, new_dataset_name="MAINZ", basal=F)
moduletrait(features=feat, metadata=metadata, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
            centrality=centrality_basal, MEs=MEs_basal, new_dataset_name="MAINZ", basal=T)



p<-modulesurv(features=feat, surv=surv, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
              centrality=centrality_global, MEs=MEs, basal=F)

pvalue_poor_tot[,5]<-p[[1]]
pvalue_good_tot[,5]<-p[[2]]

p<-modulesurv(features=feat, surv=surv, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
              centrality=centrality_basal, MEs=MEs_basal, basal=T)
pvalue_poor_tot_b[,5]<-p[[1]]
pvalue_good_tot_b[,5]<-p[[2]]


i<-5
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
data <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(data)<-dannot[rownames(data),6]

metadata<-pData(dd)
surv<-data.frame(status=metadata$e.rfs, time=metadata$t.rfs)

moduletrait(features=feat, metadata=metadata, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
            centrality=centrality_global, MEs=MEs, new_dataset_name="NKI", basal=F)
moduletrait(features=feat, metadata=metadata, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
            centrality=centrality_basal, MEs=MEs_basal, new_dataset_name="NKI", basal=T)

p<-modulesurv(features=feat, surv=surv, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
              centrality=centrality_global, MEs=MEs, basal=F)

pvalue_poor_tot[,1]<-p[[1]]
pvalue_good_tot[,1]<-p[[2]]

p<-modulesurv(features=feat, surv=surv, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
              centrality=centrality_basal, MEs=MEs_basal, basal=T)
pvalue_poor_tot_b[,1]<-p[[1]]
pvalue_good_tot_b[,1]<-p[[2]]



i<-6
dd <- get(data(list=dn[i]))
#Extract expression set, pData, fData for each dataset
data <- (exprs(dd))
dannot <- featureData(dd)@data
rownames(data)<-dannot[rownames(data),3]

metadata<-pData(dd)
surv<-data.frame(status=metadata$e.dmfs, time=metadata$t.dmfs)

moduletrait(features=feat, metadata=metadata, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
            centrality=centrality_global, MEs=MEs, new_dataset_name="VDX", basal=F)
moduletrait(features=feat, metadata=metadata, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
            centrality=centrality_basal, MEs=MEs_basal, new_dataset_name="VDX", basal=T)


p<-modulesurv(features=feat, surv=surv, original_dataset=metabric, genes_annot = dannot, new_dataset=data, 
              centrality=centrality_global, MEs=MEs, basal=F)

pvalue_poor_tot[,6]<-p[[1]]
pvalue_good_tot[,6]<-p[[2]]

p<-modulesurv(features=feat, surv=surv, original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = dannot, new_dataset=data, 
              centrality=centrality_basal, MEs=MEs_basal, basal=T)
pvalue_poor_tot_b[,6]<-p[[1]]
pvalue_good_tot_b[,6]<-p[[2]]

#TCGA

load("../../../Data/Data_downloaded/TCGA/TCGA-BRCA_Primary Tumor_dedupl.RData")
load("../../../Data/Data_downloaded/TCGA/TCGA-BRCA_clinical.RData")

library(TCGAbiolinks)
BRCA_subtype <- TCGAquery_subtype(tumor = "BRCA")

anno<-data.frame(barcode=colnames(TCGA_dedupl), 
                 subtype=BRCA_subtype$BRCA_Subtype_PAM50[match(substr(colnames(TCGA_dedupl),1,12), BRCA_subtype$patient)])
alldata_basal<-TCGA_dedupl[,which(anno$subtype=="Basal")]

# one row per patient (GDCquery_clinic); the follow-up table has several rows per patient
tcga_endpoint<-"OS"   # "OS" (GDCquery_clinic) or "PFI" (TCGA-CDR, Liu et al. 2018; sensitivity)
clin<-GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
m<-match(substr(colnames(TCGA_dedupl),1,12), clin$submitter_id)
surv<-data.frame(status=as.numeric(clin$vital_status[m]=="Dead"))
surv$time<-as.numeric(ifelse(surv$status==1, clin$days_to_death[m], clin$days_to_last_follow_up[m]))
if(tcga_endpoint=="PFI"){
  cdr<-openxlsx::read.xlsx("data/TCGA-CDR-SupplementalTableS1.xlsx"); m<-match(substr(colnames(TCGA_dedupl),1,12), cdr$bcr_patient_barcode)
  surv<-data.frame(status=as.numeric(cdr$PFI[m]), time=as.numeric(cdr$PFI.time[m]))
}

p<-modulesurv(features=feat, surv=surv[,c("status", "time")], original_dataset=metabric, genes_annot = dannot, new_dataset=TCGA_dedupl, 
              centrality=centrality_global, MEs=MEs, basal=F)


pvalue_poor_tot[,7]<-p[[1]]
pvalue_good_tot[,7]<-p[[2]]

p<-modulesurv(features=feat, surv=surv[which(anno$subtype=="Basal"),c("status", "time")], original_dataset=metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], genes_annot = NA, new_dataset=TCGA_dedupl[,which(anno$subtype=="Basal")], 
              centrality=centrality_basal, MEs=MEs_basal, basal=F)


pvalue_poor_tot_b[,7]<-p[[1]]
pvalue_good_tot_b[,7]<-p[[2]]


#################################
#########METABRIC global
#################################

###moduletrait
feat<-c("grade", "age_at_diagnosis")

moduleTraitCor = list();
moduleTraitPvalue = list();
# Calculate the correlations
for (set in feat){
  moduleTraitCor[[set]] = cor(MEs, as.numeric(as.character(meta[,set])), use = "p");
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], sum(!is.na(as.numeric(as.character(meta[,set])))));
}

moduleTraitCor_tot<-moduleTraitCor[[1]]
moduleTraitPvalue_tot<-moduleTraitPvalue[[1]]

for (set in 2:length(feat)){
  moduleTraitCor_tot<-cbind(moduleTraitCor_tot,moduleTraitCor[[set]])
  moduleTraitPvalue_tot<-cbind(moduleTraitPvalue_tot, moduleTraitPvalue[[set]])
}
colnames(moduleTraitCor_tot)<-feat
colnames(moduleTraitPvalue_tot)<-feat

save(moduleTraitCor_tot, file="results/2025_fix/moduleTraitCor_METABRIC.RData")
save(moduleTraitPvalue_tot, file="results/2025_fix/moduleTraitPvalue_METABRIC.RData")

##########survival

pvalue_good<-c()
pvalue_poor<-c()

surv<-Complete_METABRIC_Clinical_Survival_Data__DSS
surv[surv[,1]>1825,"status"]<-0
surv[surv[,1]>1825,"time"]<-1825


for(module in unique(centrality_global$module)){
  a=data.frame(ID=colnames(metabric),OS=surv[,1], IND=surv[,2])
  b=data.frame(ID=rownames(MEs), eigen_gene=MEs[,module])
  c=cbind(a,b)
  
  eigen_class<-cut(MEs[,module], breaks=quantile(MEs[,module],probs=c(0,1,0.5), na.rm=T))
  
  pp<-survdiff(Surv(OS,IND)~eigen_class, data=c)
  pval<-1-pchisq(pp$chisq, length(pp$n)-1)
  
  if(pp$obs[1]>pp$exp[1]){
    pvalue_good<-c(pvalue_good, pval)
    pvalue_poor<-c(pvalue_poor, NA)
  } else {
    pvalue_good<-c(pvalue_good, NA)
    pvalue_poor<-c(pvalue_poor, pval)
  }
  
}

names(pvalue_poor)<-unique(centrality_global$module)
names(pvalue_good)<-unique(centrality_global$module)

pvalue_poor_tot[,8]<-pvalue_poor
pvalue_good_tot[,8]<-pvalue_good

#################################
#########METABRIC basal
#################################

###moduletrait
moduleTraitCor = list();
moduleTraitPvalue = list();
# Calculate the correlations
for (set in feat){
  moduleTraitCor[[set]] = cor(MEs_basal, as.numeric(as.character(meta[meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal",set])), use = "p");
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], sum(!is.na(as.numeric(as.character(meta[meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal",set])))));
}

moduleTraitCor_tot<-moduleTraitCor[[1]]
moduleTraitPvalue_tot<-moduleTraitPvalue[[1]]

for (set in 2:length(feat)){
  moduleTraitCor_tot<-cbind(moduleTraitCor_tot,moduleTraitCor[[set]])
  moduleTraitPvalue_tot<-cbind(moduleTraitPvalue_tot, moduleTraitPvalue[[set]])
}
colnames(moduleTraitCor_tot)<-feat
colnames(moduleTraitPvalue_tot)<-feat

save(moduleTraitCor_tot, file="results/2025_fix/moduleTraitCor_METABRIC_basal.RData")
save(moduleTraitPvalue_tot, file="results/2025_fix/moduleTraitPvalue_METABRIC_basal.RData")


############survival

pvalue_good<-c()
pvalue_poor<-c()

surv<-Complete_METABRIC_Clinical_Survival_Data__DSS
surv[surv[,1]>1825,"status"]<-0
surv[surv[,1]>1825,"time"]<-1825


for(module in unique(centrality_basal$module)){
  a=data.frame(ID=colnames(metabric[,meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"]),OS=surv[meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal",1], IND=surv[meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal",2])
  b=data.frame(ID=rownames(MEs_basal), eigen_gene=MEs_basal[,module])
  c=cbind(a,b)
  
  eigen_class<-cut(MEs_basal[,module], breaks=quantile(MEs_basal[,module],probs=c(0,1,0.5), na.rm=T))
  
  pp<-survdiff(Surv(OS,IND)~eigen_class, data=c)
  pval<-1-pchisq(pp$chisq, length(pp$n)-1)
  
  if(pp$obs[1]>pp$exp[1]){
    pvalue_good<-c(pvalue_good, pval)
    pvalue_poor<-c(pvalue_poor, NA)
  } else {
    pvalue_good<-c(pvalue_good, NA)
    pvalue_poor<-c(pvalue_poor, pval)
  }
  
}

names(pvalue_poor)<-unique(centrality_basal$module)
names(pvalue_good)<-unique(centrality_basal$module)

pvalue_poor_tot_b[,8]<-pvalue_poor
pvalue_good_tot_b[,8]<-pvalue_good



###########################Overall representation of survival data
##########################Diversi tipi di survival nei diversi dataset -> DSS, DMFS, RFS

pvalue_good_tot[pvalue_good_tot<2.2*10^(-16)]<-2.2*10^(-16)
pvalue_poor_tot[pvalue_poor_tot<2.2*10^(-16)]<-2.2*10^(-16)

pvalue_merge<-(-log10(pvalue_good_tot))
pvalue_merge[is.na(pvalue_poor_tot)==F]<-(log10(pvalue_poor_tot)[is.na(pvalue_poor_tot)==F])

rownames(pvalue_merge)<-unique(centrality_global$module)
colnames(pvalue_merge)<-c("NKI", "TRANSBIG", "UNT", "UPP", "MAINZ", "VDX", "TCGA","METABRIC")

save(pvalue_merge, file="results/2025_fix/pvalue_merge_surv_global.RData")

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths


png("results/2025_fix/Suvr_7datasets_5years.png", res=500, 4000, 3000)
pheatmap(pvalue_merge[-which(rownames(pvalue_merge)=="Unconnected"),], cellwidth = 15, cellheight = 15, color = myColor)
dev.off()
pdf("results/2025_fix/Suvr_7datasets_5years.pdf", 20, 15)
pheatmap(pvalue_merge[-which(rownames(pvalue_merge)=="Unconnected"),], cellwidth = 15, cellheight = 15, color = myColor)
dev.off()

pvalue_good_tot_b[pvalue_good_tot_b<2.2*10^(-16)]<-2.2*10^(-16)
pvalue_poor_tot_b[pvalue_poor_tot_b<2.2*10^(-16)]<-2.2*10^(-16)

pvalue_merge_b<-(-log10(pvalue_good_tot_b))
pvalue_merge_b[is.na(pvalue_poor_tot_b)==F]<-(log10(pvalue_poor_tot_b)[is.na(pvalue_poor_tot_b)==F])

rownames(pvalue_merge_b)<-unique(centrality_basal$module)
colnames(pvalue_merge_b)<-c("NKI", "TRANSBIG", "UNT", "UPP", "MAINZ", "VDX", "TCGA","METABRIC")

save(pvalue_merge_b, file="results/2025_fix/pvalue_merge_surv_basal.RData")

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths


png("results/2025_fix/Suvr_7datasets_5years_b.png", res=500, 4000, 3000)
pheatmap(pvalue_merge_b[-which(rownames(pvalue_merge_b)=="b_Unconnected"),], cellwidth = 15, cellheight = 15, color = myColor)
dev.off()
pdf("results/2025_fix/Suvr_7datasets_5years_b.pdf",  20, 15)
pheatmap(pvalue_merge_b[-which(rownames(pvalue_merge_b)=="b_Unconnected"),], cellwidth = 15, cellheight = 15, color = myColor)
dev.off()