#  Identifying Noncoding Regulatory Variation by Multiome Single-Cell Sequencing in Prostate Cells

![gitcover](https://github.com/Yijun-Tian/multiome/assets/56854730/f0b74b55-425f-4d9b-a69f-8c7558a30e8d)


## **ectInsBam.pl used to extract ATAC-seq Tn5 insertion region sequence from BAM files**
e.g.: `perl ectInsBam.pl atac_possorted_bam.bam sth.Tn5ins sth.id.txt`  
Argument 1 is the ATAC-seq BAM file name;  
Argument 2 is the prefix for sam output name;  
Argument 3 is a cell barcode id file used to filter alignments.  
Remember to adjust the SAMtools executable PATH.  

## **cutmat.sh used to generate cut cell matrix per SNP from 10x intermediate files**
e.g." `bash cutmat.sh 0.01`
Argument 1 is a percentage defining how many at least cell should be informatic about each SNP being cut or not.  
Remember to adjust the 10x folder and the BEDtools PATH. The target SNP BED file is also adjustable. The BED file needs to be sorted and looks like below.  
```
$ head -n 3 dbSnp153_1000Genomes_snvGt5Pct.bed
chr1    11007   11008   rs575272151
chr1    11011   11012   rs544419019
chr1    13115   13116   rs62635286
```

## **cutCompare_*_GEXmat.R used to calculate the differentially expressed genes between cells with and without a cut on each SNP locus**
make sure to use right version of Seurat and SeuratObject:
```
> sessionInfo()
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 8.3 (Ootpa)

Matrix products: default
BLAS/LAPACK: /app/eb/software/FlexiBLAS/3.0.4-GCC-11.2.0/lib64/libflexiblas.so.3.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
 [1] doParallel_1.0.17    iterators_1.0.13     foreach_1.5.2
 [4] rtracklayer_1.54.0   org.Hs.eg.db_3.14.0  AnnotationDbi_1.56.2
 [7] Biobase_2.54.0       GenomicRanges_1.46.1 GenomeInfoDb_1.30.1
[10] IRanges_2.28.0       S4Vectors_0.32.4     BiocGenerics_0.40.0
[13] data.table_1.15.2    lubridate_1.9.3      forcats_1.0.0
[16] purrr_1.0.2          readr_2.1.4          tidyr_1.3.0
[19] tibble_3.2.1         ggplot2_3.4.4        tidyverse_2.0.0
[22] stringr_1.5.1        dplyr_1.1.4          SeuratObject_4.1.3
[25] Seurat_4.3.0
```

