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
The result file is a tab deliminated matrix file looked like below:
```
$ head BPH1e.cut.mat.txt | cut -f 1-5 | column -t
id                           AAACAGCCAATAACCT-1  AAACAGCCACAACCTA-1  AAACAGCCAGGTTTGC-1  AAACAGCCAGTTTACG-1
NA                           NA                  NA                  NA
chr10_102628952_rs4917979    NA                  NA                  NA                  NA
chr10_102669988_rs7094325    NA                  NA                  NA                  NA
chr10_102753291_rs11191385   NA                  NA                  NA                  NA
chr10_103528051_rs7080969    NA                  NA                  NA                  NA
chr10_103918331_rs111652908  NA                  NA                  NA                  NA
chr10_124907806_rs7920451    NA                  NA                  NA                  NA
chr10_124907855_rs4478950    1                   NA                  NA                  1
chr10_125000679_rs11598549   NA                  NA                  NA                  NA
```

## **cutCompare_*_GEXmat.R used to calculate the differentially expressed genes between cells with and without a cut on each SNP locus**
make sure to use right version of Seurat (Seurat_4.3.0) and SeuratObject (SeuratObject_4.1.3)
The result looks like below:
```
$ head -n 6 BPH1e.results.txt | column -t
rawid            t                    avg_log2FC            pct.1  pct.2  t_adj               myAUC  auc_diff               auc_power            wil                wil_adj  bm                 bm_adj  nbn                nbn_adj  ps                 ps_adj  lr                 lr_adj  snps       cut_cell  snp_cord         genes            symbol  QTLname
ENSG00000059915  0.645881002467781    -0.00355745538752639  0.003  0.006  1                   0.499  -0.0024658401718317    0.002                0.492810527937918  1        1                  1       0.639292384880882  1        0.576873707322748  1       0.599472268364337  1       rs4917979  175       chr10_102628952  ENSG00000059915  PSD     rs4917979_ENSG00000059915
ENSG00000076685  0.83234715089356     -0.00121380816211492  0.986  0.983  1                   0.503  -0.000841347705310858  0.00600000000000001  0.886767694488973  1        0.926716129619521  1       0.97944499364679   1        0.976797884225056  1       0.825170072484646  1       rs4917979  175       chr10_102628952  ENSG00000076685  NT5C2   rs4917979_ENSG00000076685
ENSG00000077150  0.41986987967329     -0.0190418625596354   0.046  0.063  1                   0.492  -0.0131988133458213    0.016                0.311996360675463  1        0                  0       0.428007967681019  1        0.41200972572952   1       0.391840525639293  1       rs4917979  175       chr10_102628952  ENSG00000077150  NFKB2   rs4917979_ENSG00000077150
ENSG00000107859  0.00269501382870481  0.00155602572073766   0.001  0      0.0673753457176202  0.501  0.00107855484120807    0.002                0.663928462507187  1        NA                 NA      0.994039221539202  1        0.990851581832557  1       0.540842233634519  1       rs4917979  175       chr10_102628952  ENSG00000107859  PITX3   rs4917979_ENSG00000107859
ENSG00000107862  0.281650435524757    -0.0830735536397731   0.674  0.674  1                   0.473  -0.057582199484504     0.054                0.19530304550519   1        0.170445612716625  1       0.162396605970932  1        0.118752334155827  1       0.256574667065844  1       rs4917979  175       chr10_102628952  ENSG00000107862  GBF1    rs4917979_ENSG00000107862
```

