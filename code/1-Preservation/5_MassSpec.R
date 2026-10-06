#Testing within module correlations in proteomic data
## Permutation test: median within-module correlation (unique pairs)
## vs 1000 random protein sets of the same size, drawn from proteins mapped to network genes
load(file="results/2025/centrality_global.RData")
load(file="results/2025/centrality_basal.RData")

outdir <- "results/2025_fix/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

B         <- 1000   # number of random sets
min_n     <- 3      # minimum number of mapped proteins to test a module
min_valid <- 0.5    # minimum fraction of samples with non-missing values to keep a protein

network_genes <- rownames(centrality_global)   # null pool: all network genes, Unconnected included

################################
## Data preparation
################################
## one row per gene: if several rows map to the same gene, the one with most non-missing values is kept
collapse_by_gene <- function(expr, gene){
  keep <- !is.na(gene) & gene != "" & gene %in% network_genes
  expr <- expr[keep, , drop = FALSE]; gene <- gene[keep]
  o <- order(-rowSums(!is.na(expr)))
  expr <- expr[o, , drop = FALSE]; gene <- gene[o]
  first <- !duplicated(gene)
  expr <- expr[first, , drop = FALSE]
  rownames(expr) <- gene[first]
  expr
}

## MaxQuant proteinGroups: remove reverse/contaminant hits, keep per-sample ratios only, log2
prep_maxquant <- function(path){
  d <- read.csv(path, sep = "\t", stringsAsFactors = FALSE, check.names = TRUE)
  bad <- rep(FALSE, nrow(d))
  for (cc in intersect(c("Reverse", "Potential.contaminant", "Contaminant", "Only.identified.by.site"), colnames(d)))
    bad <- bad | d[[cc]] == "+"
  ratio_cols <- grep("^Ratio\\.H\\.L\\.normalized\\.", colnames(d))   # excludes the summary column "Ratio.H.L.normalized"
  ratio <- as.matrix(d[, ratio_cols])
  ratio <- log2(ratio)
  ratio[!is.finite(ratio)] <- NA
  gene_col <- if ("Gene.names" %in% colnames(d)) "Gene.names" else colnames(d)[7]
  message(basename(path), ": gene column = ", gene_col, "; ", length(ratio_cols), " samples")
  gene <- sub(";.*$", "", as.character(d[[gene_col]]))   # protein groups with multiple genes: first gene
  ok <- !bad & rowMeans(!is.na(ratio)) >= min_valid
  collapse_by_gene(ratio[ok, , drop = FALSE], gene[ok])
}

## RPPA TCGA L4 (already log scale): columns 1-4 are metadata, antibodies in columns
prep_rppa <- function(path){
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  expr <- t(as.matrix(d[, 5:ncol(d)]))
  ok <- rowMeans(!is.na(expr)) >= min_valid
  collapse_by_gene(expr[ok, , drop = FALSE], rownames(expr)[ok])
}

################################
## Permutation test
################################
upper_median <- function(C, idx) { m <- C[idx, idx]; median(m[upper.tri(m)], na.rm = TRUE) }

perm_test_modules <- function(expr, modules, seed){
  set.seed(seed)
  C <- cor(t(expr), use = "pairwise.complete.obs")   # computed once on the whole pool
  mods <- setdiff(unique(modules), "Unconnected")
  res <- data.frame(module = mods, n = NA, median_r = NA, null_mean = NA, null_sd = NA,
                    delta = NA, z = NA, p = NA, stringsAsFactors = FALSE)
  for (k in seq_along(mods)) {
    idx <- which(rownames(expr) %in% names(modules)[modules == mods[k]])
    res$n[k] <- length(idx)
    if (length(idx) < min_n) next
    obs  <- upper_median(C, idx)
    null <- replicate(B, upper_median(C, sample(nrow(C), length(idx))))
    res$median_r[k]  <- obs
    res$null_mean[k] <- mean(null, na.rm = TRUE)
    res$null_sd[k]   <- sd(null, na.rm = TRUE)
    res$delta[k]     <- obs - res$null_mean[k]
    res$z[k]         <- res$delta[k] / res$null_sd[k]
    res$p[k]         <- (1 + sum(null >= obs, na.rm = TRUE)) / (1 + sum(!is.na(null)))   # empirical, one-sided
  }
  res$FDR <- p.adjust(res$p, method = "BH")
  res
}

modules_global <- setNames(centrality_global$module, rownames(centrality_global))

datasets <- list(
  PXD002619 = prep_maxquant("data/Mass Spec BC/proteinGroups_PXD002619.txt"),
  PXD000815 = prep_maxquant("data/Mass Spec BC/proteinGroups_PXD000815.txt"),
  RPPA_TCGA = prep_rppa("data/Mass Spec BC/RPPA TCGA/TCGA-BRCA-L4.csv"))
seeds <- c(PXD002619 = 54987656, PXD000815 = 59465496, RPPA_TCGA = 2648696)
print(sapply(datasets, dim))   # proteins x samples after filtering

results <- lapply(names(datasets), function(nm) {
  r <- perm_test_modules(datasets[[nm]], modules_global, seeds[nm])
  r$dataset <- nm
  r
})
results <- do.call(rbind, results)
write.csv(results, paste0(outdir, "protein_preservation_global_modules.csv"), row.names = FALSE)
print(results)

################################
## Heatmap: observed median minus null mean; * FDR < 0.05
################################
library(pheatmap)
delta <- tapply(results$delta, list(results$module, results$dataset), identity)
fdr   <- tapply(results$FDR,   list(results$module, results$dataset), identity)
stars <- ifelse(!is.na(fdr) & fdr < 0.05, "*", "")

paletteLength <- 50
myColor <- colorRampPalette(c("#4575B4", "white", "#D73027"))(paletteLength)
lim <- max(abs(delta), na.rm = TRUE)
myBreaks <- seq(-lim, lim, length.out = paletteLength + 1)

png(paste0(outdir, "Protein_preservation_global_modules.png"), res = 300, 1700, 2000)
pheatmap(delta, cluster_cols = F, cluster_rows = F,   # NAs for untestable modules prevent clustering
         cellwidth = 15, cellheight = 15, breaks = myBreaks, color = myColor,
         display_numbers = stars, na_col = "grey90")
dev.off()

################################
## Heatmap: -log10(empirical p); * FDR < 0.05
################################
## minimum p = 1/(B+1), hence maximum -log10(p) = log10(B+1) (about 3 with B = 1000)
p_mat <- tapply(results$p, list(results$module, results$dataset), identity)
logp  <- -log10(p_mat)
myColorP  <- colorRampPalette(c("white", "#D73027"))(paletteLength)
myBreaksP <- seq(0, log10(B + 1), length.out = paletteLength + 1)

png(paste0(outdir, "Protein_preservation_global_modules_log10p.png"), res = 300, 1700, 2000)
pheatmap(logp, cluster_cols = F, cluster_rows = F, cellwidth = 15, cellheight = 15,
         breaks = myBreaksP, color = myColorP, display_numbers = stars, na_col = "grey90")
dev.off()