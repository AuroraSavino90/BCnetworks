################################################################################
# CENTRALITY VS EXTERNAL PROGNOSIS — VERSION 2.1
#
# METABRIC supplies the global-network kWithin values.
# TRANSBIG, UNT, UPP, MAINZ, NKI and VDX supply independent survival effects.
#
# Outputs:
#   01_gene_cox_meta.csv
#   02_module_tests.csv  (poor/good x significance/effect size, as in Fig. 3)
#   03_random_module_tests.csv
#   04_global_null_summary.csv
#   figures in the same output directory
################################################################################

suppressPackageStartupMessages({
  library(Biobase)
  library(survival)
  library(dplyr)
  library(MASS)
  library(metafor)
  library(ggplot2)
  library(pheatmap)
})

set.seed(123)

CENTRALITY_FILE <- file.path("results", "2025", "centrality_global.RData")
OUT <- file.path("results", "2025_fix", "centrality_vs_external_prognosis")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

N_PERM <- 2000L
N_RANDOM <- 1000L
MIN_GENES <- 20L
MIN_COHORTS <- 3L
MIN_EVENTS <- 10L
MATCH_BLOCK <- 10L
N_BINS <- 5L
RECOMPUTE_COX <- FALSE

cohorts <- data.frame(
  dataset = c("TRANSBIG", "UNT", "UPP", "MAINZ", "NKI", "VDX"),
  package = c(
    "breastCancerTRANSBIG", "breastCancerUNT", "breastCancerUPP",
    "breastCancerMAINZ", "breastCancerNKI", "breastCancerVDX"
  ),
  object = c("transbig", "unt", "upp", "mainz", "nki", "vdx"),
  symbol_col = c(3L, 3L, 3L, 3L, 6L, 3L),
  time_col = c("t.rfs", "t.rfs", "t.rfs", "t.dmfs", "t.rfs", "t.dmfs"),
  event_col = c("e.rfs", "e.rfs", "e.rfs", "e.dmfs", "e.rfs", "e.dmfs"),
  endpoint = c("RFS", "RFS", "RFS", "DMFS", "RFS", "DMFS"),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\s*\\(.*\\)$", "", x)
  toupper(sub("\\.[0-9]+$", "", x))
}

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

to_event <- function(x) {
  if (is.numeric(x) || is.logical(x)) return(as.numeric(x))
  y <- toupper(trimws(as.character(x)))
  out <- suppressWarnings(as.numeric(y))
  out[y %in% c("TRUE", "YES", "EVENT", "DEAD", "RELAPSE")] <- 1
  out[y %in% c("FALSE", "NO", "CENSORED", "ALIVE", "NONE")] <- 0
  out
}

collapse_probes <- function(expr, symbols) {
  symbols <- clean_gene(symbols)
  v <- apply(expr, 1, var, na.rm = TRUE)
  keep <- !is.na(symbols) & symbols != "" & is.finite(v)
  expr <- expr[keep, , drop = FALSE]
  symbols <- symbols[keep]
  v <- v[keep]
  o <- order(symbols, -v)
  expr <- expr[o, , drop = FALSE]
  symbols <- symbols[o]
  keep <- !duplicated(symbols)
  expr <- expr[keep, , drop = FALSE]
  rownames(expr) <- symbols[keep]
  expr
}

load_dataset <- function(spec) {
  if (!requireNamespace(spec$package, quietly = TRUE)) {
    stop("Missing package: ", spec$package)
  }
  
  e <- new.env(parent = globalenv())
  data(list = spec$object, package = spec$package, envir = e)
  eset <- get(spec$object, envir = e)
  
  expr <- Biobase::exprs(eset)
  fd <- Biobase::fData(eset)
  pd <- Biobase::pData(eset)
  
  if (spec$symbol_col > ncol(fd)) {
    stop(spec$dataset, ": invalid symbol_col.")
  }
  
  expr <- collapse_probes(
    expr,
    fd[rownames(expr), spec$symbol_col]
  )
  
  pd <- pd[colnames(expr), , drop = FALSE]
  
  list(
    expr = expr,
    time = suppressWarnings(as.numeric(as.character(pd[[spec$time_col]]))),
    event = to_event(pd[[spec$event_col]])
  )
}

cox_gene <- function(x, time, event) {
  keep <- is.finite(x) & is.finite(time) & is.finite(event)
  if (sum(keep) < 20L || sum(event[keep] == 1) < MIN_EVENTS) return(NULL)
  
  x <- x[keep]
  if (!is.finite(sd(x)) || sd(x) == 0) return(NULL)
  
  fit <- tryCatch(
    survival::coxph(
      survival::Surv(time[keep], event[keep]) ~ scale(x),
      ties = "efron"
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  tab <- summary(fit)$coefficients
  
  data.frame(
    beta = tab[1, "coef"],
    se = tab[1, "se(coef)"],
    z = tab[1, "z"],
    p = tab[1, "Pr(>|z|)"],
    n = sum(keep),
    events = sum(event[keep] == 1),
    mean_expr = mean(x),
    sd_expr = sd(x),
    missing = 1 - sum(keep) / length(event)
  )
}

run_dataset <- function(spec) {
  cache <- file.path(OUT, paste0("cox_", spec$dataset, ".rds"))
  if (file.exists(cache) && !RECOMPUTE_COX) return(readRDS(cache))
  
  message("Cox models: ", spec$dataset)
  d <- load_dataset(spec)
  ans <- vector("list", nrow(d$expr))
  
  for (i in seq_len(nrow(d$expr))) {
    if (i %% 1000L == 0L) message(spec$dataset, ": ", i, "/", nrow(d$expr))
    a <- cox_gene(d$expr[i, ], d$time, d$event)
    if (!is.null(a)) {
      a$gene <- rownames(d$expr)[i]
      ans[[i]] <- a
    }
  }
  
  ans <- dplyr::bind_rows(ans) |>
    dplyr::mutate(
      dataset = spec$dataset,
      endpoint = spec$endpoint,
      mean_rank = rank(mean_expr) / dplyr::n(),
      sd_rank = rank(sd_expr) / dplyr::n()
    )
  
  saveRDS(ans, cache)
  ans
}

meta_gene <- function(d) {
  d <- d[is.finite(d$beta) & is.finite(d$se) & d$se > 0, , drop = FALSE]
  if (nrow(d) < MIN_COHORTS) return(NULL)
  
  fit <- tryCatch(
    metafor::rma.uni(
      yi = beta,
      sei = se,
      data = d,
      method = "REML",
      test = "knha"
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  
  beta <- as.numeric(fit$b[1])
  se <- as.numeric(fit$se[1])
  
  data.frame(
    meta_beta = beta,
    meta_se = se,
    meta_z = beta / se,
    meta_p = as.numeric(fit$pval[1]),
    CI_low = as.numeric(fit$ci.lb[1]),
    CI_high = as.numeric(fit$ci.ub[1]),
    I2 = fit$I2,
    n_cohorts = nrow(d),
    total_n = sum(d$n),
    total_events = sum(d$events),
    mean_rank = mean(d$mean_rank),
    sd_rank = mean(d$sd_rank),
    missing = mean(d$missing)
  )
}

robust_residuals <- function(d, outcome) {
  covars <- c("mean_rank", "sd_rank", "missing", "n_cohorts")
  
  covars <- covars[vapply(
    covars,
    function(v) {
      current_sd <- sd(d[[v]], na.rm = TRUE)
      is.finite(current_sd) && current_sd > 0
    },
    logical(1)
  )]
  
  required <- c("z_kWithin", outcome, covars)
  d <- d[complete.cases(d[, required, drop = FALSE]), , drop = FALSE]
  if (nrow(d) < MIN_GENES) return(NULL)
  
  if (length(covars) == 0L) {
    rx <- d$z_kWithin
    ry <- d[[outcome]]
  } else {
    rhs <- paste(covars, collapse = " + ")
    
    rx <- residuals(MASS::rlm(
      as.formula(paste("z_kWithin ~", rhs)),
      data = d,
      maxit = 100
    ))
    
    ry <- residuals(MASS::rlm(
      as.formula(paste(outcome, "~", rhs)),
      data = d,
      maxit = 100
    ))
  }
  
  # rlm residuals can retain matrix dimensions. Force plain numeric vectors;
  # otherwise cor() may return a 1 x 1 matrix, which cannot be compared with
  # the vector of N_PERM null correlations.
  rx <- as.numeric(rx)
  ry <- as.numeric(ry)
  
  rho <- suppressWarnings(
    cor(
      rx,
      ry,
      method = "spearman",
      use = "complete.obs"
    )
  )
  
  list(
    data = d,
    rx = rx,
    ry = ry,
    rho = unname(as.numeric(rho)[1])
  )
}

matching_blocks <- function(d) {
  covars <- c("mean_rank", "sd_rank", "missing", "n_cohorts")
  covars <- covars[vapply(
    covars,
    function(v) sd(d[[v]], na.rm = TRUE) > 0,
    logical(1)
  )]
  
  score <- if (length(covars) == 0L) {
    seq_len(nrow(d))
  } else {
    fitted(lm(
      as.formula(paste("z_kWithin ~", paste(covars, collapse = " + "))),
      data = d
    ))
  }
  
  o <- order(score)
  block <- integer(nrow(d))
  block[o] <- ceiling(seq_along(o) / MATCH_BLOCK)
  block
}

permute_blocks <- function(x, block) {
  out <- x
  for (b in unique(block)) {
    i <- which(block == b)
    if (length(i) > 1L) out[i] <- sample(x[i])
  }
  out
}

rho_p <- function(rho, n) {
  rho <- pmax(pmin(rho, 0.999999), -0.999999)
  t <- rho * sqrt((n - 2) / (1 - rho^2))
  2 * pt(abs(t), n - 2, lower.tail = FALSE)
}

make_strata <- function(d) {
  bin <- function(x) {
    if (sd(x, na.rm = TRUE) == 0) return(rep(1L, length(x)))
    dplyr::ntile(rank(x, na.last = "keep"), N_BINS)
  }
  interaction(bin(d$mean_rank), bin(d$sd_rank), bin(d$missing), drop = TRUE)
}

sample_matched <- function(real, universe) {
  need <- table(real$stratum)
  out <- lapply(names(need), function(s) {
    candidates <- universe[
      universe$stratum == s & !(universe$gene %in% real$gene),
      ,
      drop = FALSE
    ]
    if (nrow(candidates) < need[[s]]) {
      candidates <- universe[universe$stratum == s, , drop = FALSE]
    }
    if (nrow(candidates) == 0L) candidates <- universe
    candidates[
      sample(
        seq_len(nrow(candidates)),
        need[[s]],
        replace = nrow(candidates) < need[[s]]
      ),
      ,
      drop = FALSE
    ]
  })
  dplyr::bind_rows(out)
}

# ------------------------------------------------------------------------------
# 1. Continuous Cox models in independent cohorts
# ------------------------------------------------------------------------------

cox <- dplyr::bind_rows(lapply(
  seq_len(nrow(cohorts)),
  function(i) run_dataset(cohorts[i, , drop = FALSE])
))

gene_meta <- dplyr::bind_rows(lapply(
  split(cox, cox$gene),
  function(d) {
    x <- meta_gene(d)
    if (is.null(x)) return(NULL)
    x$gene <- d$gene[1]
    x
  }
))

write.csv(gene_meta, file.path(OUT, "01_gene_cox_meta.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 2. Join external prognosis to global METABRIC centrality
# ------------------------------------------------------------------------------

load(CENTRALITY_FILE)
if (!exists("centrality_global")) stop("centrality_global not found.")

centrality <- as.data.frame(centrality_global)
centrality$gene <- clean_gene(rownames(centrality))

analysis_data <- centrality |>
  dplyr::filter(!grepl("unconnected|grey", module, ignore.case = TRUE)) |>
  dplyr::distinct(gene, module, .keep_all = TRUE) |>
  dplyr::group_by(module) |>
  dplyr::mutate(z_kWithin = zscore(kWithin)) |>
  dplyr::ungroup() |>
  dplyr::inner_join(gene_meta, by = "gene") |>
  dplyr::mutate(
    direction = dplyr::case_when(meta_beta > 0 ~ "Poor", meta_beta < 0 ~ "Good"),
    y_pval = -log10(pmax(meta_p, .Machine$double.xmin)),
    y_effect = abs(meta_beta)
  ) |>
  dplyr::filter(!is.na(direction))

analysis_data$stratum <- make_strata(analysis_data)

write.csv(
  analysis_data,
  file.path(OUT, "01_centrality_with_external_prognosis.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 3. Module-level Spearman tests, robust adjustment and matched permutation
#    As in Fig. 3: genes split by direction (Poor: meta_beta > 0; Good: meta_beta < 0)
#    and two versions of the association with prognosis:
#      pval:   significance, -log10(meta p)
#      effect: effect size, |meta log HR per SD|
# ------------------------------------------------------------------------------

METRICS <- c(pval = "y_pval", effect = "y_effect")
METRIC_LABELS <- c(pval = "-log10(meta-analytic p-value)", effect = "|meta-analytic log HR per SD|")
COVARS <- c("mean_rank", "sd_rank", "missing", "n_cohorts")

test_module_direction <- function(d) {
  d <- d[complete.cases(d[, c("z_kWithin", METRICS, COVARS)]), , drop = FALSE]
  if (nrow(d) < MIN_GENES) return(NULL)
  
  fits <- lapply(METRICS, function(y) robust_residuals(d, y))
  if (any(vapply(fits, is.null, logical(1)))) return(NULL)
  
  # rx (residual centrality) is the same for both metrics: same genes, same covariates
  rx <- fits[[1]]$rx
  block <- matching_blocks(fits[[1]]$data)
  null <- matrix(NA_real_, N_PERM, length(METRICS), dimnames = list(NULL, names(METRICS)))
  
  for (i in seq_len(N_PERM)) {
    # Permute residual centrality within expression-matched blocks.
    px <- permute_blocks(rx, block)
    for (k in names(METRICS)) {
      null[i, k] <- unname(as.numeric(
        cor(px, fits[[k]]$ry, method = "spearman", use = "complete.obs")
      )[1])
    }
  }
  
  result <- data.frame(module = d$module[1], direction = d$direction[1], n_genes = nrow(d))
  for (k in names(METRICS)) {
    result[[paste0("rho_", k)]] <- fits[[k]]$rho
    result[[paste0("p_", k)]] <- (1 + sum(abs(null[, k]) >= abs(fits[[k]]$rho))) / (N_PERM + 1)
  }
  
  null <- data.frame(
    module = d$module[1],
    direction = d$direction[1],
    permutation = seq_len(N_PERM),
    n_genes = nrow(d),
    rho_pval = null[, "pval"],
    rho_effect = null[, "effect"]
  )
  
  list(result = result, null = null)
}

groups <- split(
  analysis_data,
  list(analysis_data$module, analysis_data$direction),
  drop = TRUE
)

module_objects <- lapply(groups, test_module_direction)
module_objects <- module_objects[!vapply(module_objects, is.null, logical(1))]

# BH within each version across all module x direction tests (as in Fig. 3C)
module_results <- dplyr::bind_rows(lapply(module_objects, function(x) x$result)) |>
  dplyr::mutate(
    FDR_pval = p.adjust(p_pval, method = "BH"),
    FDR_effect = p.adjust(p_effect, method = "BH")
  )

module_null <- dplyr::bind_rows(lapply(module_objects, function(x) x$null))

write.csv(
  module_results,
  file.path(OUT, "02_module_tests.csv"),
  row.names = FALSE
)

# Interpretation (as in Fig. 3C):
# Poor, rho > 0: among genes associated with poor prognosis, central genes are
#                more strongly associated (pval: more significant; effect: larger HR).
# Good, rho > 0: the same among genes associated with good prognosis.
# rho < 0: central genes are less strongly associated than peripheral ones.

# ------------------------------------------------------------------------------
# 4. Size- and expression-matched random gene sets (same direction)
# ------------------------------------------------------------------------------

random_results <- lapply(groups, function(real) {
  if (nrow(real) < MIN_GENES) return(NULL)
  universe <- analysis_data[analysis_data$direction == real$direction[1], , drop = FALSE]
  
  spearman <- function(d, y) {
    unname(as.numeric(
      cor(as.numeric(d$z_kWithin), as.numeric(d[[y]]), method = "spearman", use = "complete.obs")
    )[1])
  }
  
  observed <- vapply(METRICS, function(y) spearman(real, y), numeric(1))
  null <- matrix(NA_real_, N_RANDOM, length(METRICS), dimnames = list(NULL, names(METRICS)))
  
  for (i in seq_len(N_RANDOM)) {
    random <- sample_matched(real, universe)
    null[i, ] <- vapply(METRICS, function(y) spearman(random, y), numeric(1))
  }
  
  out <- data.frame(module = real$module[1], direction = real$direction[1], n_genes = nrow(real))
  for (k in names(METRICS)) {
    out[[paste0("observed_", k)]] <- observed[[k]]
    out[[paste0("p_random_", k)]] <-
      (1 + sum(abs(null[, k]) >= abs(observed[[k]]))) / (N_RANDOM + 1)
  }
  out
})

random_results <- dplyr::bind_rows(random_results) |>
  dplyr::mutate(
    FDR_random_pval = p.adjust(p_random_pval, method = "BH"),
    FDR_random_effect = p.adjust(p_random_effect, method = "BH")
  )

write.csv(
  random_results,
  file.path(OUT, "03_random_module_tests.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. Is the number of modules with a significant correlation (in at least one
#    direction) above the null?
#    Observed and null counts use the same method (parametric Spearman p + BH).
# ------------------------------------------------------------------------------

n_sig_modules <- function(d, rho) {
  sig <- p.adjust(rho_p(rho, d$n_genes), method = "BH") < 0.05
  length(unique(d$module[sig]))
}

null_counts <- lapply(split(module_null, module_null$permutation), function(d) {
  data.frame(
    permutation = d$permutation[1],
    n_pval = n_sig_modules(d, d$rho_pval),
    n_effect = n_sig_modules(d, d$rho_effect)
  )
}) |>
  dplyr::bind_rows()

observed_counts <- c(
  pval = n_sig_modules(module_results, module_results$rho_pval),
  effect = n_sig_modules(module_results, module_results$rho_effect)
)
n_modules <- length(unique(module_results$module))

global_summary <- data.frame(
  analysis = names(METRICS),
  n_modules = n_modules,
  observed_n = unname(observed_counts),
  observed_proportion = unname(observed_counts) / n_modules,
  null_mean_n = c(mean(null_counts$n_pval), mean(null_counts$n_effect)),
  null_95_percentile = c(
    quantile(null_counts$n_pval, 0.95),
    quantile(null_counts$n_effect, 0.95)
  ),
  p_global = c(
    (1 + sum(null_counts$n_pval >= observed_counts[["pval"]])) / (N_PERM + 1),
    (1 + sum(null_counts$n_effect >= observed_counts[["effect"]])) / (N_PERM + 1)
  )
)

write.csv(
  global_summary,
  file.path(OUT, "04_global_null_summary.csv"),
  row.names = FALSE
)

write.csv(
  null_counts,
  file.path(OUT, "04_null_significant_module_counts.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Figures (same layout as Fig. 3A-C)
# ------------------------------------------------------------------------------

plot_gene_module <- function(module_name, dir, metric) {
  d <- analysis_data[
    analysis_data$module == module_name & analysis_data$direction %in% dir, ,
    drop = FALSE
  ]
  if (nrow(d) == 0L) return(NULL)
  
  ggplot(d, aes(z_kWithin, .data[[METRICS[[metric]]]])) +
    geom_point(alpha = 0.55) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_classic() +
    labs(
      x = "Standardized kWithin in METABRIC",
      y = paste0(dir, " prognosis\n", METRIC_LABELS[[metric]]),
      title = module_name
    )
}

examples <- data.frame(
  module = c("E2F_TARGETS", "ESTROGEN_RESPONSE_EARLY"),
  direction = c("Poor", "Good")
)

for (k in names(METRICS)) {
  for (j in seq_len(nrow(examples))) {
    p <- plot_gene_module(examples$module[j], examples$direction[j], k)
    if (!is.null(p)) {
      ggsave(
        file.path(OUT, paste0(examples$module[j], "_", tolower(examples$direction[j]), "_", k, ".pdf")),
        p,
        width = 5,
        height = 5
      )
    }
  }
  
  # Heatmap as in Fig. 3C: only FDR < 0.05, symmetric colour scale
  mods <- unique(module_results$module)
  cor_mat <- matrix(
    NA_real_, length(mods), 2,
    dimnames = list(mods, c("Poor prognosis", "Good prognosis"))
  )
  for (r in seq_len(nrow(module_results))) {
    if (module_results[[paste0("FDR_", k)]][r] < 0.05) {
      cor_mat[module_results$module[r], paste(module_results$direction[r], "prognosis")] <-
        module_results[[paste0("rho_", k)]][r]
    }
  }
  
  if (!all(is.na(cor_mat))) {
    cor_mat <- cor_mat[order(rowSums(cbind(cor_mat[, 1], -cor_mat[, 2]), na.rm = TRUE), decreasing = TRUE), ]
    M <- max(abs(cor_mat), na.rm = TRUE)
    pheatmap::pheatmap(
      cor_mat,
      cellwidth = 15, cellheight = 15,
      breaks = seq(-M, M, length.out = 51),
      color = colorRampPalette(c("#4575B4", "white", "#D73027"))(50),
      cluster_rows = FALSE, cluster_cols = FALSE,
      main = METRIC_LABELS[[k]],
      filename = file.path(OUT, paste0("module_summary_", k, "_pheat.pdf")),
      width = 10, height = 10
    )
  }
  
  p_null <- ggplot(null_counts, aes(.data[[paste0("n_", k)]])) +
    geom_histogram(binwidth = 1, boundary = -0.5) +
    geom_vline(xintercept = observed_counts[[k]], linetype = "dashed") +
    theme_classic() +
    labs(
      x = "Number of modules with FDR < 0.05 under the null",
      y = "Permutations",
      title = paste0(METRIC_LABELS[[k]], " - observed: ", observed_counts[[k]], "/", n_modules, " modules")
    )
  
  ggsave(
    file.path(OUT, paste0("global_null_distribution_", k, ".pdf")),
    p_null,
    width = 7,
    height = 5
  )
}

saveRDS(
  list(
    gene_meta = gene_meta,
    analysis_data = analysis_data,
    module_results = module_results,
    random_results = random_results,
    global_summary = global_summary
  ),
  file.path(OUT, "analysis_summary.rds")
)

message("Completed: ", normalizePath(OUT, winslash = "/", mustWork = FALSE))