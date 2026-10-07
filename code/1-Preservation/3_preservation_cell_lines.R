#####Comparison with lines modules
library(WGCNA)
library(ensembldb)
library(EnsDb.Hsapiens.v75)

#change names in gene symbols  (unchanged)
changenames<-function(data, anno){
  annotation_sel=anno[match( rownames(data), anno[,1]),2]
  
  if(length(which(annotation_sel==""))>0){
    data<-data[-which(annotation_sel==""),]
    annotation_sel<-annotation_sel[-which(annotation_sel=="")]
  }
  
  a<-which(duplicated(annotation_sel))
  while(length(a)>0){
    for(i in 1:length(unique(annotation_sel))){
      if(length(which(annotation_sel==unique(annotation_sel)[i]))>1){
        m=which.max(rowMeans(data[which(annotation_sel==unique(annotation_sel)[i]),], na.rm=T))
        data=data[-which(annotation_sel==unique(annotation_sel)[i])[-m],]
        annotation_sel=annotation_sel[-which(annotation_sel==unique(annotation_sel)[i])[-m]]
        print(i)
        print(dim(data))
      }
    }
    
    data=data[which(is.na(annotation_sel)==F),]
    annotation_sel=na.omit(annotation_sel)
    a<-which(duplicated(annotation_sel))
  }
  
  rownames(data)=annotation_sel
  return(data)
}

load(file="data/RData/metabric.RData")
load(file="data/RData/meta.RData")
load(file="data/RData/net_metabric_oneblock.RData")
load(file="data/RData/net_metabric_Basal.RData")
load(file="results/2025/centrality_global.RData")
load(file="results/2025/centrality_basal.RData")

## separate output folder, so that the original results are not overwritten
outdir <- "results/2025_fix/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

################################
## 1. Sample files (GSM only) and series matrix annotation, matched by GSM
################################
gse_dir <- "data/GSE48213/"     # no setwd: explicit paths
files <- list.files(gse_dir, pattern = "^GSM[0-9]+")
gsm_files <- sub("^(GSM[0-9]+).*", "\\1", files)
stopifnot(!anyDuplicated(gsm_files))

meta_lines <- read.csv("data/GSE48213_series_matrix.csv", header = F, stringsAsFactors = FALSE)
# series matrix rows identified by content, not by position
row_with <- function(pattern) {
  hit <- sapply(meta_lines, function(col) grepl(pattern, col))
  which(rowSums(hit) >= length(files))
}
gsm_row <- row_with("^GSM[0-9]+$")[1]
sub_row <- row_with("^subtype:")
stopifnot(length(gsm_row) == 1, !is.na(gsm_row), length(sub_row) == 1)
gsm_meta <- as.character(unlist(meta_lines[gsm_row, ]))
sub_meta <- as.character(unlist(meta_lines[sub_row, ]))
is_sample <- grepl("^GSM[0-9]+$", gsm_meta)       # excludes the row-label column, if present
line_info <- data.frame(gsm = gsm_meta[is_sample],
                        subtype = sub("^subtype: *", "", sub_meta[is_sample]),
                        stringsAsFactors = FALSE)
stopifnot(setequal(gsm_files, line_info$gsm))      # same samples in the files and in the series matrix
print(table(line_info$subtype))

################################
## 2. Read the files, checking that genes are in the same order
################################
x <- NULL
for(i in seq_along(files)){
  y <- read.csv(file.path(gse_dir, files[i]), row.names = 1, header = T, sep = "\t")
  stopifnot(ncol(y) == 1)
  if (!is.null(x)) stopifnot(identical(rownames(y), rownames(x)))   # cbind does not align by rownames
  if (is.null(x)) x <- y else x <- cbind(x, y)
}
x <- as.matrix(x)
colnames(x) <- gsm_files

################################
## 3. log2 transformation
################################
x <- log2(x + 1)

################################
## 4. Ensembl gene ID -> gene symbol with EnsDb.Hsapiens.v75 (Ensembl 75, GRCh37; to be cited in the Methods)
##    chosen for its higher coverage of basal module genes compared with the most recent release
################################
edb <- EnsDb.Hsapiens.v75
print(ensemblVersion(edb))
ids <- sub("\\..*$", "", rownames(x))               # removes the ID version suffix, if present
sym <- mapIds(edb, keys = ids, keytype = "GENEID", column = "SYMBOL")
anno <- data.frame(ensembl_gene_id = rownames(x), hgnc_symbol = as.character(sym),
                   stringsAsFactors = FALSE)
sc_data <- changenames(x, anno)
mod_genes <- rownames(centrality_basal)[centrality_basal$module != "b_Unconnected"]
print(c(mapped_genes = nrow(sc_data),
        in_metabric = length(intersect(rownames(sc_data), rownames(metabric))),
        basal_module_genes_present = sum(mod_genes %in% rownames(sc_data))))

################################
## 5. Basal + claudin-low lines selected by GSM
################################
# basal and claudin-low lines (distinct categories in the GSE48213 annotation)
basal_gsm <- line_info$gsm[line_info$subtype %in% c("Basal", "Claudin-low")]
sc_data_basal <- sc_data[, colnames(sc_data) %in% basal_gsm]
print(files[gsm_files %in% basal_gsm])
# comparison with the selection of the original script (by position)
basal_lines_old <- which(meta_lines[11, ] == "subtype: Basal")
print(list(original = files[basal_lines_old], corrected = files[gsm_files %in% basal_gsm]))

################################
## 6. Module preservation (same parameters as in 1_preservation_bulkBC.R)
################################
plot_preservation <- function(mp, file){
  ref = 1
  test = 2
  statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
  statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])
  print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
               signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
  Z <- mp$preservation$Z[[ref]][[test]]
  plotMods <- !(rownames(Z) %in% c("gold", "Unconnected", "b_Unconnected"))
  moduleSizes <- Z[plotMods, "moduleSize"]
  zs <- Z[plotMods, "Zsummary.pres"]
  min = min(zs, na.rm = TRUE)
  max = max(zs, na.rm = TRUE)
  if (min > -max/10) min = -max/10
  ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  pdf(file)
  plot(moduleSizes, zs, col = 1, bg = "grey", pch = 21,
       main = "Preservation Zsummary", cex = 2.4,
       ylab = "Preservation Zsummary", xlab = "Module size", log = "x",
       ylim = ylim, cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.4)
  labelPoints(moduleSizes, zs, rownames(Z)[plotMods], cex = 1, offs = 0.08)
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
  dev.off()
}

run_preservation <- function(ref_expr, ref_modules, test_expr, test_label, file_stub){
  # remove genes/samples with too many missing values or zero variance in each set (modulePreservation does not remove them itself)
  gsg_ref  <- goodSamplesGenes(t(ref_expr),  verbose = 0)
  gsg_test <- goodSamplesGenes(t(test_expr), verbose = 0)
  message(file_stub, ": removed ", sum(!gsg_ref$goodGenes), " genes / ", sum(!gsg_ref$goodSamples), " samples (reference), ",
          sum(!gsg_test$goodGenes), " genes / ", sum(!gsg_test$goodSamples), " samples (test)")
  ref_expr    <- ref_expr[gsg_ref$goodGenes, gsg_ref$goodSamples]
  ref_modules <- ref_modules[gsg_ref$goodGenes]      # module colors kept aligned with the retained genes
  test_expr   <- test_expr[gsg_test$goodGenes, gsg_test$goodSamples]
  maxSize <- max(table(ref_modules[!ref_modules %in% c("Unconnected", "b_Unconnected")]))
  multiExpr = list(BC = list(data = t(ref_expr)), test = list(data = t(test_expr)))
  names(multiExpr)[2] <- test_label
  multiColor = list(BC = ref_modules)
  print(system.time( {
    mp = modulePreservation(multiExpr, multiColor,
                            referenceNetworks = 1,
                            networkType = "signed",
                            maxModuleSize = maxSize,
                            nPermutations = 200,
                            randomSeed = 1,
                            quickCor = 0,
                            verbose = 3,
                            checkData = TRUE)
  } ))
  save(mp, file = paste0(outdir, "modulePreservation_metabricVs", file_stub, ".RData"))
  plot_preservation(mp, paste0(outdir, "preservation_metabricVs", file_stub, ".pdf"))
  invisible(mp)
}

## basal modules, basal + claudin-low lines (file name unchanged for 4_preservation_pheatmap.R)
run_preservation(metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], centrality_basal$module,
                 sc_data_basal, "lines", "lines_basal")

## global modules, all lines
run_preservation(metabric, centrality_global$module, sc_data, "lines", "lines_global")
