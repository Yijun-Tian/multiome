#  Brief title for this repo
##__ectInsBam.pl: used to extract ATAC-seq Tn5 insertion region sequence from BAM files__
Argument 1 is the ATAC-seq BAM file name
Argument 2 is the prefix for sam output name
Argument 3 is a cell barcode id file for filter alignments
Example: perl ectInsBam.pl atac_possorted_bam.bam sth.Tn5ins sth.id.txt
