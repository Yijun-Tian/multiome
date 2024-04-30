.libPaths('/home/4470655/R/x86_64-pc-linux-gnu-library/4.1/')
library(Seurat)
library(dplyr)
library(stringr)
library(tidyverse)
counts1 <- Read10X_h5("/share/lab_wangl/Yijun_Tian/LabNotes/GEX-ATAC-Seq/twistCapture/multiome_library_capture/cellrangerarc_quantify/BPH1e/outs/filtered_feature_bc_matrix.h5",use.names=FALSE)
rna_counts1 <- counts1$`Gene Expression`
BPH1 <- CreateSeuratObject(counts = rna_counts1, min.cells=1)
BPH1 <- subset(BPH1,cells=scan("BPH1e.id.txt",what=""))
BPH1 <- SCTransform(BPH1, verbose = FALSE)
gex<-as.data.frame(BPH1@assays$SCT@data)
gexMat=gex[!rownames(gex)%in%scan("MT.geneid.txt",what=""),]
eMat=gexMat[rownames(gexMat)%in%scan("e.geneid.txt",what=""),]
write.table(data.frame("geneid"=rownames(eMat), eMat, check.names = FALSE), "GEXmat.eENSEMBL.BPH1e.txt", row.names=F, quote=F, sep="\t")
write.table(data.frame("geneid"=rownames(gexMat), gexMat, check.names = FALSE), "GEXmat.ENSEMBL.BPH1e.txt", row.names=F, quote=F, sep="\t")


