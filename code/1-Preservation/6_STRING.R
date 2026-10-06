## Protein-protein interaction connectivity of global modules in STRING (v12)
## - interactions scored using the "experiments" and "database" channels only
##   (the coexpression channel is excluded to avoid circularity with expression-based modules)
## - modules and random sets restricted to genes mapped to STRING, so that random sets have the same number of mapped proteins
## - empirical p-values against two nulls: (i) uniform random sets, (ii) degree-matched random sets (controls for study bias)
library(STRINGdb)
library(data.table)
library(Matrix)
library(ggplot2)

load("results/2025/centrality_global.RData")

outdir <- "results/2025_fix/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

B         <- 1000   # number of random sets per null
threshold <- 400    # minimum recombined experiments+database score (STRING "medium confidence")
n_bins    <- 10     # degree bins for the degree-matched null
min_n     <- 3      # minimum number of mapped genes to test a module

################################
## 1. Map all network genes to STRING once
################################
string_db <- STRINGdb$new(version = "12", species = 9606, score_threshold = 0, input_directory = "")
network_genes <- rownames(centrality_global)
mapped <- string_db$map(data.frame(gene = network_genes), "gene", removeUnmappedRows = TRUE)
mapped <- mapped[!duplicated(mapped$STRING_id) & !duplicated(mapped$gene), ]   # one-to-one gene/protein mapping
message("Network genes mapped to STRING: ", nrow(mapped), " / ", length(network_genes))
pool <- mapped$STRING_id                                    # null pool: all mapped network genes, Unconnected included
gene_modules <- centrality_global[mapped$gene, "module"]

################################
## 2. Channel-specific scores (experiments + database), recombined as in STRING
################################
links_file <- "data/9606.protein.links.detailed.v12.0.txt.gz"
if (!file.exists(links_file))
  download.file("https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz",
                links_file, mode = "wb")
links <- fread(links_file, select = c("protein1", "protein2", "experimental", "database"))
links <- links[protein1 %in% pool & protein2 %in% pool & (experimental > 0 | database > 0)]

## STRING combination rule: remove the prior from each channel, combine as independent evidence, add the prior back
prior <- 0.041
recombine <- function(...) {
  p_not <- Reduce(`*`, lapply(list(...), function(s) {
    s <- pmax((s / 1000 - prior) / (1 - prior), 0)
    1 - s
  }))
  ((1 - p_not) * (1 - prior) + prior) * 1000
}
links[, score := recombine(experimental, database)]
links <- links[score >= threshold]

## symmetric sparse matrix of scores among pool proteins (each pair is listed twice in the STRING file)
A <- sparseMatrix(i = match(links$protein1, pool), j = match(links$protein2, pool), x = links$score,
                  dims = c(length(pool), length(pool)))
degree <- rowSums(A > 0)
deg_bin <- cut(rank(degree, ties.method = "random"), n_bins, labels = FALSE)   # equal-size degree bins

################################
## 3. Permutation tests
################################
set_score <- function(idx) sum(A[idx, idx]) / 2      # total score of interacting pairs within the set

sample_degree_matched <- function(idx) {
  unlist(lapply(split(idx, deg_bin[idx]), function(v) {
    b <- deg_bin[v[1]]
    sample(which(deg_bin == b), length(v))
  }))
}

set.seed(683650)
mods <- setdiff(unique(gene_modules), "Unconnected")
res <- data.frame(module = mods, n_mapped = NA, observed = NA,
                  null_uniform_mean = NA, fold_uniform = NA, p_uniform = NA,
                  null_degree_mean = NA, fold_degree = NA, p_degree = NA, stringsAsFactors = FALSE)
for (k in seq_along(mods)) {
  idx <- which(gene_modules == mods[k])
  res$n_mapped[k] <- length(idx)
  if (length(idx) < min_n) next
  obs <- set_score(idx)
  null_u <- replicate(B, set_score(sample(length(pool), length(idx))))
  null_d <- replicate(B, set_score(sample_degree_matched(idx)))
  res$observed[k]          <- obs
  res$null_uniform_mean[k] <- mean(null_u)
  res$fold_uniform[k]      <- obs / mean(null_u)
  res$p_uniform[k]         <- (1 + sum(null_u >= obs)) / (1 + B)   # empirical, one-sided
  res$null_degree_mean[k]  <- mean(null_d)
  res$fold_degree[k]       <- obs / mean(null_d)
  res$p_degree[k]          <- (1 + sum(null_d >= obs)) / (1 + B)
}
res$FDR_uniform <- p.adjust(res$p_uniform, method = "BH")
res$FDR_degree  <- p.adjust(res$p_degree,  method = "BH")
write.csv(res, paste0(outdir, "STRING_PPI_global_modules.csv"), row.names = FALSE)
print(res)

################################
## 4. Plots: observed score vs mean of each null (log scale); filled points FDR < 0.05
################################
plot_null <- function(null_col, fdr_col, xlab, file) {
  df <- data.frame(null = res[[null_col]], obs = res$observed, module = res$module,
                   sig = ifelse(res[[fdr_col]] < 0.05, "FDR < 0.05", "n.s."))
  df <- df[!is.na(df$obs), ]
  p <- ggplot(df, aes(null, obs, label = module, shape = sig)) + geom_point(size = 4) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1.5) +
    scale_x_log10() + scale_y_log10() + scale_shape_manual(values = c("FDR < 0.05" = 16, "n.s." = 1)) +
    theme_bw() + xlab(xlab) + ylab("STRING score module genes") + labs(shape = NULL)
  pdf(file, 10, 10); print(p + geom_text(vjust = -1)); dev.off()
}
plot_null("null_uniform_mean", "FDR_uniform", "STRING score random genes (mean)",
          paste0(outdir, "STRING_PPI_uniform_null.pdf"))
plot_null("null_degree_mean", "FDR_degree", "STRING score degree-matched random genes (mean)",
          paste0(outdir, "STRING_PPI_degree_matched_null.pdf"))
