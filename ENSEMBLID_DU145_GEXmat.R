.libPaths('/home/4470655/R/x86_64-pc-linux-gnu-library/4.1/')
library(Seurat)
library(dplyr)
library(stringr)
library(tidyverse)
counts1 <- Read10X_h5("/share/lab_wangl/Yijun_Tian/LabNotes/GEX-ATAC-Seq/twistCapture/multiome_library_capture/cellrangerarc_quantify/DU145e/outs/filtered_feature_bc_matrix.h5",use.names=FALSE)
rna_counts1 <- counts1$`Gene Expression`
DU145 <- CreateSeuratObject(counts = rna_counts1, min.cells=1)
DU145 <- subset(DU145,cells=scan("DU145e.id.txt",what=""))
DU145 <- SCTransform(DU145, verbose = FALSE)
gex<-as.data.frame(DU145@assays$SCT@data)
gexMat=gex[!rownames(gex)%in%scan("MT.geneid.txt",what=""),]
eMat=gexMat[rownames(gexMat)%in%scan("e.geneid.txt",what=""),]
write.table(data.frame("geneid"=rownames(eMat), eMat, check.names = FALSE), "GEXmat.eENSEMBL.DU145e.txt", row.names=F, quote=F, sep="\t")
write.table(data.frame("geneid"=rownames(gexMat), gexMat, check.names = FALSE), "GEXmat.ENSEMBL.DU145e.txt", row.names=F, quote=F, sep="\t")