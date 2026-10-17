######################
## Load precomputed module preservations (corrected runs)
######################
library(pheatmap)

indir  <- "results/2025_fix/"
outdir <- "results/2025_fix/"

## file stubs as saved by 1_, 2_ and 3_preservation scripts; names are the column labels
datasets <- c("In TCGA"          = "TCGA",
              "In NKI"           = "NKI",
              "In TRANSBIG"      = "TRANSBIG",
              "In UNT"           = "UNT",
              "In UPP"           = "UPP",
              "In MAINZ"         = "MAINZ",
              "In VDX"           = "VDX",
              "In BC cell lines" = "lines")
global_stub <- c(TCGA = "TCGA", NKI = "NKI_global", TRANSBIG = "TRANSBIG", UNT = "UNT", UPP = "UPP",
                 MAINZ = "MAINZ", VDX = "VDX", lines = "lines_global")

## Zsummary.pres of each module in one dataset (selected by column name, not position)
get_Z <- function(file) {
  load(file)                                  # loads object "mp"
  Z <- mp$preservation$Z[[1]][[2]]
  setNames(Z[, "Zsummary.pres"], rownames(Z))
}

## matrix of Zsummary (modules x datasets), matched by module name across datasets
build_matrix <- function(type) {
  stubs <- if (type == "global") global_stub[datasets] else paste0(datasets, "_basal")
  z <- lapply(stubs, function(s) get_Z(paste0(indir, "modulePreservation_metabricVs", s, ".RData")))
  mods <- Reduce(union, lapply(z, names))
  mods <- setdiff(mods, c("grey", "gold", "Unconnected", "b_Unconnected"))   # remove gold and unassigned genes
  m <- sapply(z, function(v) v[mods])
  rownames(m) <- mods
  colnames(m) <- names(datasets)
  m
}

######################
## Color scale on Zsummary itself, with breaks at the standard thresholds
## Z < 2: no evidence of preservation (blue to white); 2 <= Z < 10: moderate (light to mid red);
## Z >= 10: strong evidence of preservation (dark red). Values above 10 are capped for display only.
######################
z_colors <- function(Zmin) {
  low  <- seq(min(Zmin, 0), 2, length.out = 11)
  mid  <- seq(2, 10, length.out = 21)[-1]
  breaks <- c(low, mid, 10.01)
  cols <- c(colorRampPalette(c("#4575B4", "white"))(length(low) - 1),
            colorRampPalette(c("#FDDBC7", "#D6604D"))(length(mid)),
            "#B2182B")
  list(breaks = breaks, colors = cols)
}

plot_preservation_heatmap <- function(Zmat, file_stub) {
  write.csv(Zmat, paste0(outdir, file_stub, "_Zsummary.csv"))      # uncapped values for the supplementary tables
  Zplot <- pmin(Zmat, 10.005)   # capped values fall in the top (Z >= 10) color class
  sc <- z_colors(min(Zmat, na.rm = TRUE))
  numbers <- matrix(ifelse(is.na(Zmat), "", formatC(Zmat, format = "f", digits = 1)),
                    nrow = nrow(Zmat), dimnames = dimnames(Zmat))   # actual (uncapped) Z in each cell
  args <- list(Zplot, cluster_cols = F, cluster_rows = !anyNA(Zplot),
               cellwidth = 18, cellheight = 15, breaks = sc$breaks, color = sc$colors,
               legend_breaks = c(0, 2, 5, 10), legend_labels = c("0", "2", "5", "\u226510"),
               display_numbers = numbers, fontsize_number = 5, number_color = "black",
               na_col = "grey90")
  graphics.off()
  png(paste0(outdir, file_stub, ".png"), res = 300, 2500, 2000)
  do.call(pheatmap, args)
  dev.off()
  graphics.off()
  pdf(paste0(outdir, file_stub, ".pdf"), 7, 10)
  do.call(pheatmap, args)
  dev.off()
}

##############################
#### Figure preservation global (Supplementary Figure 2A)
##############################
Z_global <- build_matrix("global")
plot_preservation_heatmap(Z_global, "preservation_global_modules_alldatasets")

##############################
#### Figure preservation basal (Figure 2E)
##############################
Z_basal <- build_matrix("basal")
plot_preservation_heatmap(Z_basal, "preservation_basal_modules_alldatasets")

## number of datasets with Z >= 2 and Z >= 10 for each module (for the text)
print(cbind(Z_ge_2 = rowSums(Z_basal >= 2, na.rm = TRUE), Z_ge_10 = rowSums(Z_basal >= 10, na.rm = TRUE)))