.libPaths('/home/4470655/R/x86_64-pc-linux-gnu-library/4.1/')
library(Seurat)
library(dplyr)
library(stringr)
library(tidyverse)
counts1 <- Read10X_h5("/share/lab_wangl/Yijun_Tian/LabNotes/GEX-ATAC-Seq/twistCapture/multiome_library_capture/cellrangerarc_quantify/RWPE2R2e/outs/filtered_feature_bc_matrix.h5",use.names=FALSE)
rna_counts1 <- counts1$`Gene Expression`
RWPE2R2 <- CreateSeuratObject(counts = rna_counts1, min.cells=1)
RWPE2R2 <- subset(RWPE2R2,cells=scan("RWPE2R2e.id.txt",what=""))
RWPE2R2 <- SCTransform(RWPE2R2, verbose = FALSE)
gex<-as.data.frame(RWPE2R2@assays$SCT@data)
gexMat=gex[!rownames(gex)%in%scan("MT.geneid.txt",what=""),]
eMat=gexMat[rownames(gexMat)%in%scan("e.geneid.txt",what=""),]
write.table(data.frame("geneid"=rownames(eMat), eMat, check.names = FALSE), "GEXmat.eENSEMBL.RWPE2R2e.txt", row.names=F, quote=F, sep="\t")
write.table(data.frame("geneid"=rownames(gexMat), gexMat, check.names = FALSE), "GEXmat.ENSEMBL.RWPE2R2e.txt", row.names=F, quote=F, sep="\t")


