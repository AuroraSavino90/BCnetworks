#####Comparison with lines modules
library(WGCNA)
library(ensembldb)
library(EnsDb.Hsapiens.v75)

#change names in gene symbols  (invariata)
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

## cartella di output separata, per non sovrascrivere i risultati originali
outdir <- "results/2025_fix/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

################################
## 1. File dei campioni (solo GSM) e annotazione dal series matrix, con match per GSM
################################
gse_dir <- "data/GSE48213/"     # niente setwd: percorsi espliciti
files <- list.files(gse_dir, pattern = "^GSM[0-9]+")
gsm_files <- sub("^(GSM[0-9]+).*", "\\1", files)
stopifnot(!anyDuplicated(gsm_files))

meta_lines <- read.csv("data/GSE48213_series_matrix.csv", header = F, stringsAsFactors = FALSE)
# righe del series matrix individuate dal contenuto, non dalla posizione
row_with <- function(pattern) {
  hit <- sapply(meta_lines, function(col) grepl(pattern, col))
  which(rowSums(hit) >= length(files))
}
gsm_row <- row_with("^GSM[0-9]+$")[1]
sub_row <- row_with("^subtype:")
stopifnot(length(gsm_row) == 1, !is.na(gsm_row), length(sub_row) == 1)
gsm_meta <- as.character(unlist(meta_lines[gsm_row, ]))
sub_meta <- as.character(unlist(meta_lines[sub_row, ]))
is_sample <- grepl("^GSM[0-9]+$", gsm_meta)       # esclude l'eventuale colonna di etichette
line_info <- data.frame(gsm = gsm_meta[is_sample],
                        subtype = sub("^subtype: *", "", sub_meta[is_sample]),
                        stringsAsFactors = FALSE)
stopifnot(setequal(gsm_files, line_info$gsm))      # stessi campioni nei file e nel series matrix
print(table(line_info$subtype))

################################
## 2. Lettura dei file, con controllo che i geni siano nello stesso ordine
################################
x <- NULL
for(i in seq_along(files)){
  y <- read.csv(file.path(gse_dir, files[i]), row.names = 1, header = T, sep = "\t")
  stopifnot(ncol(y) == 1)
  if (!is.null(x)) stopifnot(identical(rownames(y), rownames(x)))   # cbind non allinea per rownames
  if (is.null(x)) x <- y else x <- cbind(x, y)
}
x <- as.matrix(x)
colnames(x) <- gsm_files

################################
## 3. Trasformazione log2
################################
x <- log2(x + 1)

################################
## 4. Ensembl gene ID -> gene symbol con EnsDb.Hsapiens.v75 (Ensembl 75, GRCh37; da citare nel Methods)
##    scelta per la maggiore copertura dei geni dei moduli basal rispetto alla release piu' recente
################################
edb <- EnsDb.Hsapiens.v75
print(ensemblVersion(edb))
ids <- sub("\\..*$", "", rownames(x))               # rimuove l'eventuale versione dell'ID
sym <- mapIds(edb, keys = ids, keytype = "GENEID", column = "SYMBOL")
anno <- data.frame(ensembl_gene_id = rownames(x), hgnc_symbol = as.character(sym),
                   stringsAsFactors = FALSE)
sc_data <- changenames(x, anno)
mod_genes <- rownames(centrality_basal)[centrality_basal$module != "b_Unconnected"]
print(c(geni_mappati = nrow(sc_data),
        in_metabric = length(intersect(rownames(sc_data), rownames(metabric))),
        geni_moduli_basal_presenti = sum(mod_genes %in% rownames(sc_data))))

################################
## 5. Linee basal + claudin-low selezionate per GSM
################################
# linee basal e claudin-low (nell'annotazione GSE48213 sono categorie distinte)
basal_gsm <- line_info$gsm[line_info$subtype %in% c("Basal", "Claudin-low")]
sc_data_basal <- sc_data[, colnames(sc_data) %in% basal_gsm]
print(files[gsm_files %in% basal_gsm])
# confronto con la selezione dello script originale (per posizione)
basal_lines_old <- which(meta_lines[11, ] == "subtype: Basal")
print(list(originale = files[basal_lines_old], corretta = files[gsm_files %in% basal_gsm]))

################################
## 6. Module preservation (stessi parametri di 1_preservation_bulkBC.R)
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
  # rimozione di geni/campioni con troppi NA o varianza nulla in ciascun set (modulePreservation non li rimuove da se')
  gsg_ref  <- goodSamplesGenes(t(ref_expr),  verbose = 0)
  gsg_test <- goodSamplesGenes(t(test_expr), verbose = 0)
  message(file_stub, ": rimossi ", sum(!gsg_ref$goodGenes), " geni / ", sum(!gsg_ref$goodSamples), " campioni (riferimento), ",
          sum(!gsg_test$goodGenes), " geni / ", sum(!gsg_test$goodSamples), " campioni (test)")
  ref_expr    <- ref_expr[gsg_ref$goodGenes, gsg_ref$goodSamples]
  ref_modules <- ref_modules[gsg_ref$goodGenes]      # colori allineati ai geni rimasti
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

## moduli basal, linee basal + claudin-low (nome file invariato per 4_preservation_pheatmap.R)
run_preservation(metabric[, meta$NOT_IN_OSLOVAL_Pam50Subtype=="Basal"], centrality_basal$module,
                 sc_data_basal, "lines", "lines_basal")

## moduli globali, tutte le linee
run_preservation(metabric, centrality_global$module, sc_data, "lines", "lines_global")