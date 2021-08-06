#!/bin/sh
length=$1
pos=$2
    # length=15; pos=7
idlist=$( awk '{print $2}' Info/genotypeid_hiptcrid.txt )

odir=phenotype.freq3/L${length}_P${pos}
mkdir -p $odir

ofile=$odir/single_aa_usage_ratio_downsamples.txt
echo "Sample AA rate" > $ofile

for id in $idlist;do
   echo $id
   
   #reads with functionality, v gene identified, and select cdr3-aa length
   cat mixcr_v4/$id/clones.TRB.txt |
   sed -e "1d" |
   cut -f 33 |
   grep  "\*\|_" |
   awk -v L=$length 'BEGIN{FS="\t"}{
      amino_acid=$1;
      if( length(amino_acid) == L ) { print amino_acid }
   }' |
   awk -v pos=$pos '{
      amino_acid = $1;
      aa=substr(amino_acid, pos, 1);
      print aa
   }' |
   grep -v "\*\|_" > $odir/tmp.$id #target position amino acids
   
   #downsample productive sequence
   Total=$( cat $odir/tmp.$id | wc -l ) #this is the number of sequence with non-productive CDR3
   
   cat mixcr_v4/$id/clones.TRB.txt |
   sed -e "1d" |
   cut -f 33 |
   grep -v "\*\|_" |
   awk -v L=$length 'BEGIN{FS="\t"}{
      amino_acid=$1;
      if( length(amino_acid) == L ) { print amino_acid }
   }' |
   awk -v pos=$pos '{
      amino_acid = $1;
      aa=substr(amino_acid, pos, 1);
      print aa
   }' |
   shuf |
   head -n $Total > $odir/tmp2.$id
   
   #calculate ratio
   cat $odir/tmp2.$id |
   sort  | uniq -c |
   awk -v Total=$Total -v id=$id '{print id, $2, $1/Total}'  >> $ofile
   
   rm -f $odir/tmp.$id
   rm -f $odir/tmp2.$id
   
done

gzip -f $ofile



