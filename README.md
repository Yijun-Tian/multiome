#  Brief title for this repo
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
