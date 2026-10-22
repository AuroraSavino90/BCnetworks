### Extension of 4_subtype_regression.R to the external cohorts (D4).
### Same models as in METABRIC: grade ~ ME + PAM50, age ~ ME + PAM50,
### Cox(5-year) ~ ME + PAM50, with projected eigengenes (R/pcaproject.R),
### PAM50 calls obtained as in 1_clinicalFeatures_7datasets.R.
### Combination across external cohorts: Fisher and Stouffer (R/combine_p.R).

library(Biobase)
library(WGCNA)
library(genefu)
data(pam50.robust)
library(survival)
library(TCGAbiolinks)
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

# Fisher (metap::sumlog) and Stouffer (metap::sumz) combinations, all cohorts and external only
combine_directional <- function(p, effect) {
  ok <- is.finite(p) & is.finite(effect)
  if (sum(ok) < 2) return(c(k = sum(ok), p_fisher = NA, p_stouffer = NA))
  p1 <- two2one(p[ok], two = rep(TRUE, sum(ok)), invert = effect[ok] < 0)
  c(k = sum(ok), p_fisher = sumlog(p1)$p, p_stouffer = sumz(p1)$p)
}

# p, effect: modules x datasets; effect > 0 = tested direction.
combine_matrix <- function(p, effect, test, discovery_col = 1) {
  out <- list()
  for (cohort_set in c("all_cohorts", "external_only")) {
    cols <- if (cohort_set == "all_cohorts") seq_len(ncol(p)) else setdiff(seq_len(ncol(p)), discovery_col)
    res <- t(sapply(rownames(p), function(m) combine_directional(p[m, cols], effect[m, cols])))
    df <- data.frame(module = rownames(p), test = test, cohort_set = cohort_set, res, row.names = NULL)
    df$fdr_fisher <- p.adjust(df$p_fisher, method = "BH")
    df$fdr_stouffer <- p.adjust(df$p_stouffer, method = "BH")
    out[[cohort_set]] <- df
  }
  do.call(rbind, out)
}

# Inputs as built in 2_pheatmap_clinical.R:
#   grade_p/grade_r, age_p/age_r: modules x datasets, METABRIC in column 1;
#   surv_signed: pvalue_merge(_b), signed -log10 p (> 0 = good prognosis),
#   METABRIC in column "METABRIC". p-values truncated at 2.2e-16 upstream.
combine_p_tables <- function(grade_p, grade_r, age_p, age_r, surv_signed) {
  modules <- rownames(grade_p)
  surv_p <- 10^(-abs(surv_signed[modules, , drop = FALSE]))
  surv_poor <- -sign(surv_signed[modules, , drop = FALSE])  # > 0 = poor prognosis
  surv_disc <- which(colnames(surv_signed) == "METABRIC")
  rownames(age_p) <- rownames(age_r) <- modules
  rbind(
    combine_matrix(grade_p, grade_r, "high_grade"),
    combine_matrix(grade_p, -grade_r, "low_grade"),
    combine_matrix(age_p, -age_r, "early_onset"),
    combine_matrix(age_p, age_r, "late_onset"),
    combine_matrix(surv_p, surv_poor, "poor_prognosis", surv_disc),
    combine_matrix(surv_p, -surv_poor, "good_prognosis", surv_disc)
  )
}

dir.create("results/2025_fix", showWarnings=FALSE)
load("data/RData/meta.RData")
load("data/RData/metabric.RData")
load("results/2025/centrality_global.RData")
MEs= moduleEigengenes(t(metabric), centrality_global$module)$eigengenes
colnames(MEs)<-gsub("^ME","", colnames(MEs))
rownames(MEs)<-colnames(metabric)

datasets<-list(TRANSBIG=c("breastCancerTRANSBIG","transbig",3,"rfs"), UNT=c("breastCancerUNT","unt",3,"rfs"),
               UPP=c("breastCancerUPP","upp",3,"rfs"), MAINZ=c("breastCancerMAINZ","mainz",3,"dmfs"),
               NKI=c("breastCancerNKI","nki",6,"rfs"), VDX=c("breastCancerVDX","vdx",3,"dmfs"))

fit_models<-function(ME, subtype, grade, age, surv, dataset){
  modules<-setdiff(colnames(ME), "Unconnected")
  surv[which(surv[,"time"]>1825),"status"]<-0
  surv[which(surv[,"time"]>1825),"time"]<-1825
  res<-list()
  for(m in modules){
    x<-ME[,m]
    if(sum(!is.na(grade))>=30){
      s<-summary(lm(grade~x+subtype))[[4]]
      res[[length(res)+1]]<-data.frame(dataset=dataset, feature="Grade", module=m, coef=s[2,1], p=s[2,4])
    }
    if(sum(!is.na(age))>=30){
      s<-summary(lm(age~x+subtype))[[4]]
      res[[length(res)+1]]<-data.frame(dataset=dataset, feature="Age at onset", module=m, coef=s[2,1], p=s[2,4])
    }
    s<-summary(coxph(Surv(surv[,"time"], surv[,"status"])~x+subtype))[[7]]
    res[[length(res)+1]]<-data.frame(dataset=dataset, feature="Survival", module=m, coef=s[1,1], p=s[1,5])
  }
  do.call(rbind, res)
}

results<-list()
for(ds in names(datasets)){
  cfg<-datasets[[ds]]
  library(cfg[1], character.only=TRUE)
  dd <- get(data(list=cfg[2]))
  data <- exprs(dd)
  dannot <- featureData(dd)@data
  rownames(data)<-dannot[rownames(data),as.numeric(cfg[3])]
  metadata<-pData(dd)
  PAM50<-molecular.subtyping(sbt.model = "pam50", data=t(data), annot=dannot)$subtype
  expr<-collapse_probes(data, rownames(data))
  expr<-expr[which(rowSums(is.na(expr))==0),]
  ME<-pcaproject(newdata=expr, original_data=metabric, modules=centrality_global$module, ME=MEs)
  surv<-cbind(status=metadata[[paste0("e.",cfg[4])]], time=metadata[[paste0("t.",cfg[4])]])
  results[[ds]]<-fit_models(ME, as.character(PAM50), as.numeric(as.character(metadata$grade)),
                            as.numeric(as.character(metadata$age)), surv, ds)
}

# TCGA: survival only (as in 1_clinicalFeatures_7datasets.R)
load("../../../Data/Data_downloaded/TCGA/TCGA-BRCA_Primary Tumor_dedupl.RData")
BRCA_subtype <- TCGAquery_subtype(tumor = "BRCA")
subtype<-BRCA_subtype$BRCA_Subtype_PAM50[match(substr(colnames(TCGA_dedupl),1,12), BRCA_subtype$patient)]
clin<-GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")
m<-match(substr(colnames(TCGA_dedupl),1,12), clin$submitter_id)
status<-as.numeric(clin$vital_status[m]=="Dead")
surv<-cbind(status=status, time=as.numeric(ifelse(status==1, clin$days_to_death[m], clin$days_to_last_follow_up[m])))
ME<-pcaproject(newdata=TCGA_dedupl, original_data=metabric, modules=centrality_global$module, ME=MEs)
results$TCGA<-fit_models(ME, as.character(subtype), rep(NA,ncol(TCGA_dedupl)), rep(NA,ncol(TCGA_dedupl)), surv, "TCGA")

res<-do.call(rbind, results)
res$fdr<-ave(res$p, res$dataset, res$feature, FUN=function(p) p.adjust(p, method="fdr"))
write.csv(res, "results/2025_fix/subtype_regression_external_per_dataset.csv", row.names=FALSE)

# Combination across external cohorts (positive and negative association), Fisher and Stouffer
comb<-list()
for(f in unique(res$feature)){
  sub<-res[res$feature==f,]
  p<-tapply(sub$p, list(sub$module, sub$dataset), identity)
  coef<-tapply(sub$coef, list(sub$module, sub$dataset), identity)
  comb[[f]]<-rbind(combine_matrix(p, coef, paste(f, "positive"), discovery_col=integer(0)),
                   combine_matrix(p, -coef, paste(f, "negative"), discovery_col=integer(0)))
}
comb<-do.call(rbind, comb)
comb<-comb[comb$cohort_set=="all_cohorts",]   # all cohorts here are external
write.csv(comb, "results/2025_fix/subtype_regression_external_combined.csv", row.names=FALSE)