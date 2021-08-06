#!/bin/sh
wd=$( pwd )
id=$1

mkdir -p mixcr_v4/$id
cd mixcr_v4/$id

#Step1: Align sequencing reads
gunzip -c $wd/fasta/HIP/$id.fa.gz > tmp.fa

java -Xmx2g -Xms2g -jar /PHShome/ki991/tools/mixcr-2.1.11/mixcr-2.1.11/mixcr.jar \
   align \
   -s hs -f \
   --verbose \
   tmp.fa \
   alignments.vdjca

rm -f tmp.fa

#Step2: assembly without error correction
java -Xmx2g -Xms2g -jar /PHShome/ki991/tools/mixcr-2.1.11/mixcr-2.1.11/mixcr.jar \
   assemble -f \
   alignments.vdjca clones.clns

#Step3: Export
java -Xmx2g -Xms2g -jar /PHShome/ki991/tools/mixcr-2.1.11/mixcr-2.1.11/mixcr.jar \
   exportClones -f \
   -c TRA clones.clns clones.TRA.txt

java -Xmx2g -Xms2g -jar /PHShome/ki991/tools/mixcr-2.1.11/mixcr-2.1.11/mixcr.jar \
   exportClones -f \
   -c TRB clones.clns clones.TRB.txt

rm -f alignments*



