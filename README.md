#  Brief title for this repo
## **ectInsBam.pl used to extract ATAC-seq Tn5 insertion region sequence from BAM files**
Argument 1 is the ATAC-seq BAM file name;  
Argument 2 is the prefix for sam output name;  
Argument 3 is a cell barcode id file used to filter alignments.  
e.g.: `perl ectInsBam.pl atac_possorted_bam.bam sth.Tn5ins sth.id.txt`
