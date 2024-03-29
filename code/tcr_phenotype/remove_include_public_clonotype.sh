#!/bin/sh
id=$1
dd=/data/srlab/kishigaki/data/TCR_robins
    # id=HIP14156

 #1, no public sequence
mkdir -p cdr3/HP_nopub/$id

zcat $dd/$id.tsv.gz |
cut -f 2,8,9,11,17 |
sed -e "1d" |
awk 'BEGIN{FS="\t"}{
   amino_acid=$1;
   productive_frequency=$2;
   v_gene=$4;
   v_gene_org=$4;
    j_gene=$5;
   if( productive_frequency != "null" &&  v_gene != "unresolved" ) { 
      v_gene = substr(v_gene,5,7);
      split(v_gene, v_gene2, "-");
      print  v_gene2[1] "_" amino_acid, v_gene_org, j_gene }
}' |
grep -F -w -v -f cdr3/HP_nopub/pub_cdr3.txt - |
cut -d "_" -f 2  > cdr3/HP_nopub/$id/tmp_all.cdr3

for length in $(seq 12 18);do
   cat cdr3/HP_nopub/$id/tmp_all.cdr3 |
   awk -v L=$length '{if( length($1)==L ){ print }}' |
   gzip -f -c - >  cdr3/HP_nopub/$id/$length.txt.gz #target length cdr3
done

rm -f cdr3/HP_nopub/$id/tmp_all.cdr3


 #2, public sequence
mkdir -p cdr3/HP_pub/$id

zcat $dd/$id.tsv.gz |
cut -f 2,8,9,11,17 |
sed -e "1d" |
awk 'BEGIN{FS="\t"}{
   amino_acid=$1;
   productive_frequency=$2;
   v_gene=$4;
   v_gene_org=$4;
    j_gene=$5;
   if( productive_frequency != "null" &&  v_gene != "unresolved" ) { 
      v_gene = substr(v_gene,5,7);
      split(v_gene, v_gene2, "-");
      print  v_gene2[1] "_" amino_acid, v_gene_org, j_gene }
}' |
grep -F -w -f cdr3/HP_nopub/pub_cdr3.txt - |
cut -d "_" -f 2  > cdr3/HP_pub/$id/tmp_all.cdr3

for length in $(seq 12 18);do
   cat cdr3/HP_pub/$id/tmp_all.cdr3 |
   awk -v L=$length '{if( length($1)==L ){ print }}' |
   gzip -f -c - >  cdr3/HP_pub/$id/$length.txt.gz #target length cdr3
done

rm -f cdr3/HP_pub/$id/tmp_all.cdr3



