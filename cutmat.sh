for i in $( ls /share/lab_wangl/Yijun_Tian/LabNotes/GEX-ATAC-Seq/twistCapture/multiome_library_capture/cellrangerarc_quantify/RWPE[1-2]e/outs/atac_fragments.tsv.gz| cut -d"/" -f10 | sort | uniq );
do
#Get Tn5 insertion sites within each cell from atac-seq fragment ends: cutsite.*.bed
 zcat /share/lab_wangl/Yijun_Tian/LabNotes/GEX-ATAC-Seq/twistCapture/multiome_library_capture/cellrangerarc_quantify/${i}/outs/atac_fragments.tsv.gz |\
      grep -w -F -f ${i}.id.txt - | grep -v "GL" | grep -v "KI" | grep -v "MT" > atac_fragments_${i}.bed
 awk '{FS=OFS="\t"}{print $1,$2-5,$2+5,$4,$5,"_"NR"_"}' atac_fragments_${i}.bed > cutsite.$i.bed
 awk '{FS=OFS="\t"}{print $1,$3-5,$3+5,$4,$5,"_"NR"_"}' atac_fragments_${i}.bed >> cutsite.$i.bed

#Get Tn5 insertion information from SNP surrounding regions: wider.$i.cut.txt & wider.$i.c3.txt
 bedtools intersect -a cutsite.$i.bed -b <( awk '{FS=OFS="\t"}{if ($2-20 >0) print $1,$2-20,$3+20,$4; else print $1,1,$3+20,$4}' dbSnp153_1000Genomes_snvGt5Pct.bed) -wa -wb |\
      sort -k6,10 -T TEMP/ -S20% --parallel=4 - | uniq --skip-fields=5 > wider.${i}.cut.txt

 awk '{FS=OFS="\t"}{print $4"|"$7"_"$8+20"_"$10}' wider.$i.cut.txt | sort | uniq -c | sed 's/^[[:space:]]*//g' | awk '{FS=" "; OFS="\t"}{print $2,$1}' | sed "s/|/\t/" > wider.$i.c3.txt

#Get contribution of Tn5 insertion for each fragments, so that can remove them from fragment length calculation
 #cut -f6-10 ${i}.cut.txt > $i.cut.contribution

#Experimental: get fragments hovering SNP region from nucleosome-free ones: < 147 bp, may expand to more strengent criteria in the future to infer certain TF binidng
 ## bedtools intersect -a <( awk '{OFS="\t"}{if ($3-$2 <= 147) print $0,"_"NR"_"}' atac_fragments_${i}.bed |\
   ##   grep -v -F -f <(cut -f 1 $i.cut.contribution | sort | uniq ) - ) -b dbSnp153_1000Genomes_snvGt5Pct.bed -wa -wb > $i.temp.bed

#Get uniq snp lists after above two filter
 cell=$(wc -l ${i}.id.txt | cut -d " " -f1)
 thres=$(echo $cell*$1 | bc -l)
 pct=`printf "%.0f\n" $(echo 100*$1 | bc -l)`
 echo "looking for SNP in $cell $i cells with no more than $thres cells empty (detected in at least ${pct}% cells)"
 awk 'BEGIN {FS=OFS="\t"} NR>0 {evidence[$1"|"$2]+=1} END {for (m in evidence) {print m,evidence[m]}}' wider.$i.c3.txt |\
 sed -e 's/|/\t/g' | awk 'BEGIN {FS=OFS="\t"} NR>0 {snp[$2]+=1} END {for (m in snp) {print m,snp[m]}}' > ${i}.cut.profile
 awk -v b=$thres -F"\t" '{if($2 >= b) print $1}' ${i}.cut.profile > ${i}.cut.rsId.txt
 grep -w -F -f ${i}.cut.rsId.txt wider.$i.c3.txt > wider.$i.c3q.txt
 grep -v -F -f <(cut -f1 wider.$i.c3q.txt | sort | uniq ) ${i}.id.txt >> wider.$i.c3q.txt

 awk 'BEGIN {FS=OFS="\t"} NR>0 {total[$2$1]+=$3;BC[$1]++;rsID[$2]++} END {n=asorti(BC,cp1); printf "id"; for (i=1; i<=n; i++) {printf "%s%s",FS,cp1[i]}; print ""; m=asorti(rsID,cp2); for (j=1;j<=m;j++) {printf "%s", cp2[j]; for (k=1; k<=n; k++) {printf "%s%i",FS,total[cp2[j]cp1[k]]}; print ""}}' wider.$i.c3q.txt |\
   sed -e '2,$s/\t0/\tNA/g' > $i.cut.mat.txt
done






