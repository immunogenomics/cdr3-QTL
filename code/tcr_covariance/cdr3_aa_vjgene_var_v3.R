 # remove duplicate count at the center of cdr3
args <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line
dd="/data/srlab/kishigaki/data/TCR_robins"
library(data.table)
LIST=read.table("Info/genotypeid_hiptcrid.txt", stringsAsFactors=F)[,2]

df<-data.frame()

for(Sample in LIST){
   print(Sample)
   ofile <- paste0("phenotype.summary/HIP/aa_vjgene_var3/",Sample,".tcr_aapos_vgene_var.txt.gz")
   
   FILE = paste0(dd,"/",Sample,".tsv")
   d1 <- fread(FILE)
   d1 <- d1[,c(2:8,11,17),with=F]
   d1 <- as.data.frame(d1)
   d1 <- subset(d1, productive_frequency != "null")
   d1 <- subset(d1, v_gene!="unresolved")
   d1 <- subset(d1, j_gene!="unresolved")
   d1$length <- nchar(d1$amino_acid)
   d1 <- subset(d1, length==15) #MODIFY HERE
   
   #prepare matrix
   mat <- sapply(as.character(d1$amino_acid), function(x){
      x <- as.character(x);
      x <- strsplit(x, "");
      x <- unlist( x )
   })
   mat <- t(mat) #cdr3 sequence matirx: 15 col matrix
   
   #entropy
   #cdr3
   res <- data.frame()
   for( i in 1:15 ){
      x1 <- mat[, i]
      TB <- table(x1)
      rate <- c(TB) / sum(TB)
      ent <- sum( - rate * log(rate) )
      norment <- ent / log(length(rate))
      dump <- data.frame(i, ent, norment)
      res  <- rbind(res , dump)
   }
   
   #v gene
   i=0 
   x1 <- d1$v_gene
   TB <- table(x1)
   rate <- c(TB) / sum(TB)
   ent <- sum( - rate * log(rate) )
   norment <- ent / log(length(rate)) #normalize by dividing the max of entropy
   dump <- data.frame(i, ent, norment)
   res  <- rbind(res , dump)
   #j gene
   i=16
   x1 <- d1$j_gene
   TB <- table(x1)
   rate <- c(TB) / sum(TB)
   ent <- sum( - rate * log(rate) )
   norment <- ent / log(length(rate)) #normalize by dividing the max of entropy
   dump <- data.frame(i, ent, norment)
   res  <- rbind(res , dump)
   
   gz1 <- gzfile(ofile, "w")
   write.table(res, gz1, sep = "\t", quote = F, row.names = FALSE, col.names=TRUE)
   close(gz1)
   
}



