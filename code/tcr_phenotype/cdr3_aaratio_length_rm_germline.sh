#!/bin/sh
length=$1
pos=$2
    # length=15; pos=12
dd=/data/srlab/kishigaki/data/TCR_robins
idlist=$( awk '{print $2}' Info/genotypeid_hiptcrid.txt )

odir=phenotype.freq3/L${length}_P${pos}_rm_gl
mkdir -p $odir

ofile=$odir/single_aa_usage_ratio.txt
echo "Sample AA rate" > $ofile

for id in $idlist;do
   echo $id
   
   #reads with functionality, v gene identified, and select cdr3-aa length
   zcat $dd/$id.tsv.gz |
   sed -e "1d" |
   cut -f 2,8,9,11,17 |
   awk -v L=$length 'BEGIN{FS="\t"}{
      amino_acid=$1;
      productive_frequency=$2;
      v_gene=$4;
      j_gene=$5;
      if( productive_frequency != "null" && 
          v_gene != "unresolved" && 
          j_gene != "unresolved" && 
          length(amino_acid) == L ) { 
             print amino_acid, v_gene, j_gene
          }
   }' |
   grep -F -w -v -f data/missing.vgenes - |
   awk -v pos=$pos '{
      amino_acid = $1;
      aa=substr(amino_acid, pos, 1);
      print aa, $2,$3
   }' > $odir/tmp.$id #target position amino acids, vgene, j gene
   
   #exclude V germline encoded amino acids
   if [ -s data/V_P${pos}.exclude.txt ] ; then
      cat $odir/tmp.$id |
      awk '{print $2 "__" $1, $3}' |
      grep -F -v -w -f data/V_P${pos}.exclude.txt - |
      awk '{
         split($1, D, "__");
         AA = D[2]
         print $2 "__" AA
      }' > $odir/tmp.$id.v1
   else
      cat $odir/tmp.$id |
      awk '{print $2 "__" $1, $3}' |
      awk '{
         split($1, D, "__");
         AA = D[2]
         print $2 "__" AA
      }' > $odir/tmp.$id.v1
   fi
   
   if [ -s data/J_L${length}_P${pos}.exclude.txt  ] ; then
      cat $odir/tmp.$id.v1 |
      grep -F -v -w -f data/J_L${length}_P${pos}.exclude.txt - |
      awk '{
         split($1, D, "__");
         AA = D[2]
         print AA
      }' > $odir/tmp2.$id
   else
      cat $odir/tmp.$id.v1 |
      awk '{
         split($1, D, "__");
         AA = D[2]
         print AA
      }' > $odir/tmp2.$id
   fi
   
   #calculate ratio
   Total=$( cat $odir/tmp2.$id | wc -l )
   
   cat $odir/tmp2.$id |
   sort  | uniq -c |
   awk -v Total=$Total -v id=$id '{print id, $2, $1/Total}'  >> $ofile
   
   rm -f $odir/tmp.$id
   rm -f $odir/tmp2.$id
   rm -f $odir/tmp.$id.v1
   
done

gzip -f $ofile





