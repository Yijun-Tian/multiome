#  Identifying Noncoding Regulatory Variation by Multiome Single-Cell Sequencing in Prostate Cells

![gitcover](https://github.com/Yijun-Tian/multiome/assets/56854730/f0b74b55-425f-4d9b-a69f-8c7558a30e8d)


## **ectInsBam.pl used to extract ATAC-seq Tn5 insertion region sequence from BAM files**
e.g.: `perl ectInsBam.pl atac_possorted_bam.bam sth.Tn5ins sth.id.txt`  
Argument 1 is the ATAC-seq BAM file name;  
Argument 2 is the prefix for sam output name;  
Argument 3 is a cell barcode id file used to filter alignments.  
Remember to adjust the SAMtools executable PATH.  

## **cutmat.sh used to generate cut cell matrix per SNP from 10x intermediate fragment files**
e.g." `bash cutmat.sh 0.01`
Argument 1 is a percentage defining how many at least cell should be informatic about each SNP being cut or not.  
Remember to adjust the 10x folder and the BEDtools PATH. The target SNP BED file is also adjustable. The BED file needs to be sorted and looks like below.  
```
$ head -n 3 dbSnp153_1000Genomes_snvGt5Pct.bed
chr1    11007   11008   rs575272151
chr1    11011   11012   rs544419019
chr1    13115   13116   rs62635286
```
## ** ENSEMBLID_*_GEXmat.R used to generate SCT transformed gene expression matrix from 10x intermediate matrix.h5 files**
Make sure to use right version of Seurat (Seurat_4.3.0) and SeuratObject (SeuratObject_4.1.3).
"*.id.txt" can be used for customization of cells by barcode. "*.geneid.txt" can be used to remove mitochondria genes or focus on eQTL genes.

## **cutCompare_*_GEXmat.R used to calculate the differentially expressed genes between cells with and without a cut on each SNP locus**
make sure to use right version of Seurat (Seurat_4.3.0) and SeuratObject (SeuratObject_4.1.3)

