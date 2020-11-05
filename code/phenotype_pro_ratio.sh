#!/bin/sh
    # length=15
dd=/data/srlab/kishigaki/data/TCR_robins
idlist=$( awk '{print $2}' Info/genotypeid_hiptcrid.txt )

ofile=phenotype.summary/HIP/pro_rate/single_aa_usage_ratio.txt
echo "Sample length N_obs" > $ofile

for id in $idlist;do
   echo $id
   
   #reads with functionality, v gene identified, and select cdr3-aa length
   cat $dd/$id.tsv |
   sed -e "1d" |
   cut -f 2,8,9,11 |
   awk 'BEGIN{ FS="\t" }{
      amino_acid=$1;
      productive_frequency=$2;
      L= length(amino_acid)
      if( productive_frequency != "null" ) { 
         print L
      } else { print "NA" }
   }' |
   sort |
   uniq -c |
   awk -v id=$id '{print id, $2, $1}' >> $ofile #target position amino acids
   
done

gzip -f $ofile



