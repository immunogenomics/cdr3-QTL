#!/bin/sh
length=$1
pos=$2
    # length=15; pos=7
dd=/data/srlab/kishigaki/data/TCR_robins
idlist=$( awk '{print $2}' Info/genotypeid_hiptcrid.txt )

odir=phenotype.freq3/L${length}_P${pos}
mkdir -p $odir

ofile=$odir/jgene_cond_single_aa_usage_ratio.txt
echo "Sample Vgene AA rate obs" > $ofile

for id in $idlist;do
   echo $id
   
   #reads with functionality, v gene identified, and select cdr3-aa length
   cat $dd/$id.tsv |
   sed -e "1d" |
   cut -f 2,8,9,17 |
   awk -v L=$length 'BEGIN{FS="\t"}{
      amino_acid=$1;
      productive_frequency=$2;
      v_gene=$4;
      if( productive_frequency != "null" && 
          v_gene != "unresolved" && 
          length(amino_acid) == L ) { 
             print amino_acid, v_gene
          }
   }' |
   awk -v pos=$pos '{
      amino_acid = $1;
      v_gene=$2 ;
      aa=substr(amino_acid, pos, 1);
      print aa, v_gene
   }' > $odir/tmp.$id #target position amino acids
   
   #calculate ratio
   vgenelist=$( awk '{print $2}' $odir/tmp.$id | sort | uniq )
   
   for vgene in $vgenelist;do
      cat $odir/tmp.$id |
      awk -v vgene=$vgene '{
         if($2==vgene){ print $1 }
      }' > $odir/tmp.$id.2
      
      Total=$( cat $odir/tmp.$id.2 | wc -l )
      
      cat $odir/tmp.$id.2 |
      sort  | uniq -c |
      awk -v Total=$Total -v id=$id -v vgene=$vgene '{print id, vgene, $2, $1/Total, Total}'  >> $ofile
      
      rm -f $odir/tmp.$id.2
      
   done
   
   rm -f $odir/tmp.$id
   
done

gzip -f $ofile




