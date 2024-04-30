.libPaths('/home/4470655/R/x86_64-pc-linux-gnu-library/4.1/')
library(Seurat)
library(dplyr)
library(stringr)
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(org.Hs.eg.db)
library(rtracklayer)
library(foreach)
library(doParallel)

counts1 <- Read10X_h5("/share/lab_wangl/Yijun_Tian/LabNotes/GEX-ATAC-Seq/twistCapture/multiome_library_capture/cellrangerarc_quantify/BPH1e/outs/filtered_feature_bc_matrix.h5",use.names=FALSE)
rna_counts1 <- counts1$`Gene Expression`
BPH1 <- CreateSeuratObject(counts = rna_counts1, min.cells=1)
BPH1 <- subset(BPH1,cells=scan("BPH1e.id.txt",what=""))
BPH1 <- SCTransform(BPH1, verbose = FALSE)
cutmat <- fread("BPH1e.cut.mat.txt", header = T)
humGene <- import("hg38.ensembl.gene.bed")
#registerDoParallel(4)

merge.all <- function(x, ..., by = "row.names") {
  L <- list(...)
  for (i in seq_along(L)) {
    x <- merge(x, L[[i]], by = by)
    rownames(x) <- x$Row.names
    x$Row.names <- NULL
  }
  return(x)
}

res <- data.frame()
#res <- foreach (i = 2:dim(cutmat)[1], .combine = rbind, .verbose = T, .errorhandling = 'remove') %dopar% {
for (i in 2:dim(cutmat)[1]) {
  print(i)
  tmp <- unlist(strsplit(cutmat[i, 1][[1]], "_"))
  snp <- tmp[3]
  cut_cell <- names(cutmat)[-1][!is.na(cutmat[i, -1])]
  region <- GenomicRanges::GRanges(seqnames = Rle(c(tmp[1])), ranges = IRanges(as.numeric(tmp[2]) - 500000, end = as.numeric(tmp[2]) + 500000))
  genelist <- IRanges::subsetByOverlaps(humGene, region)$name
   if (length(genelist) == 0)
  {next }
  else {
  try(s <- subset(BPH1, features = genelist), silent = T)
  s@meta.data$condition <- sapply(rownames(s@meta.data), function(ita) ifelse(ita %in% cut_cell, "cut", "X"))
  s <- Seurat::SetIdent(s, value = "condition")
  slice_t <- Seurat::FindMarkers(s, ident.1 = "X", ident.2 = "cut", slot = "data", logfc.threshold = 0.00, test.use = "t", min.pct = 0.00, assay = "SCT", verbose = F) %>% dplyr::rename(t = p_val, t_adj = p_val_adj)
  slice_roc <- Seurat::FindMarkers(s, ident.1 = "X", ident.2 = "cut", slot = "data", logfc.threshold = 0.00, test.use = "roc", min.pct = 0.00, assay = "SCT", verbose = F) %>%
               dplyr::rename(auc_diff = avg_diff, auc_power = power) %>% dplyr::select(-c(avg_log2FC, pct.1, pct.2))
  slice_wilcox <- Seurat::FindMarkers(s, ident.1 = "X", ident.2 = "cut", slot = "data", logfc.threshold = 0.00, test.use = "wilcox", min.pct = 0.00, assay = "SCT", verbose = F) %>%
               dplyr::rename(wil = p_val, wil_adj = p_val_adj) %>% dplyr::select(-c(avg_log2FC, pct.1, pct.2))
  slice_bm <- Seurat::FindMarkers(s, ident.1 = "X", ident.2 = "cut", slot = "data", logfc.threshold = 0.00, test.use = "bimod", min.pct = 0.00, assay = "SCT", verbose = F) %>%
              dplyr::rename(bm = p_val, bm_adj = p_val_adj) %>% dplyr::select(-c(avg_log2FC, pct.1, pct.2))
  slice_nbn <- Seurat::FindMarkers(s, ident.1 = "X", ident.2 = "cut", slot = "data", logfc.threshold = 0.00, test.use = "negbinom", min.pct = 0.00, assay = "SCT", verbose = F) %>%
              dplyr::rename(nbn = p_val, nbn_adj = p_val_adj) %>% dplyr::select(-c(avg_log2FC, pct.1, pct.2))
  slice_ps <- Seurat::FindMarkers(s, ident.1 = "X", ident.2 = "cut", slot = "data", logfc.threshold = 0.00, test.use = "poisson", min.pct = 0.00, assay = "SCT", verbose = F) %>%
              dplyr::rename(ps = p_val, ps_adj = p_val_adj) %>% dplyr::select(-c(avg_log2FC, pct.1, pct.2))
  slice_lr <- Seurat::FindMarkers(s, ident.1 = "X", ident.2 = "cut", slot = "data", logfc.threshold = 0.00, test.use = "LR", min.pct = 0.00, assay = "SCT", verbose = F) %>%
              dplyr::rename(lr = p_val, lr_adj = p_val_adj) %>% dplyr::select(-c(avg_log2FC, pct.1, pct.2))
  #slice <- dplyr::mutate(slice_t, slice_roc, slice_wilcox, slice_bm, slice_nbn, slice_ps, slice_lr) %>% dplyr::mutate(snps = snp, cut_cell = length(cut_cell), snp_cord = paste(tmp[1], tmp[2], sep="_"))
  #Use merge.all function to avoid less gene reported in the methods
  slice <- merge.all(slice_t, slice_roc, slice_wilcox, slice_bm, slice_nbn, slice_ps, slice_lr) %>% dplyr::mutate(snps = snp, cut_cell = length(cut_cell), snp_cord = paste(tmp[1], tmp[2], sep="_"))
  #return(slice)
  entity <- rbind(res, slice)
  res = entity
  }
}

results <- rownames_to_column(res, "rawid") %>% mutate(genes = substr(rawid, 1, 15)) %>%
              mutate(symbol= mapIds(org.Hs.eg.db, keys= genes, keytype = "ENSEMBL", column="SYMBOL")) %>%
	      mutate(QTLname = paste(snps, genes, sep = "_"))

saveRDS(results, "BPH1e.results.rds")
write.table(results, row.names = F,"BPH1e.results.txt", sep = "\t", quote = F)