#!/bin/sh
length=$1
pos=$2
    # length=15; pos=7
dd=/data/srlab/kishigaki/data/TCR_robins
idlist=$( awk '{print $2}' Info/genotypeid_hiptcrid.txt )

odir=phenotype.freq3/L${length}_P${pos}_nopub
mkdir -p $odir

ofile=$odir/single_aa_usage_ratio.txt
echo "Sample AA rate" > $ofile

for id in $idlist;do
   echo $id
   
   #reads with functionality, v gene identified, and select cdr3-aa length
   zcat cdr3/HP_nopub/$id/$length.txt.gz |
   awk -v pos=$pos '{
      amino_acid = $1;
      aa=substr(amino_acid, pos, 1);
      print aa
   }' > $odir/tmp.$id #target position amino acids
   
   #calculate ratio
   Total=$( cat $odir/tmp.$id | wc -l )
   
   cat $odir/tmp.$id |
   sort  | uniq -c |
   awk -v Total=$Total -v id=$id '{print id, $2, $1/Total}'  >> $ofile
   
   rm -f $odir/tmp.$id
   
done

gzip -f $ofile




