 # cdr3 length=15
args <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line
dd="/data/srlab/kishigaki/data/TCR_robins"
library(data.table)
library(entropy)
LIST=read.table("Info/genotypeid_hiptcrid.txt", stringsAsFactors=F)[,2]

df<-data.frame()

for(Sample in LIST){
   print(Sample)
   ofile <- paste0("phenotype.summary/HIP/aa_pos_cor3/",Sample,".tcr_each_aapos_cor.txt.gz")
   
   FILE = paste0(dd,"/",Sample,".tsv")
   d1 <- fread(FILE)
   d1 <- d1[,2:8,with=F]
   d1 <- as.data.frame(d1)
   d1$length <- nchar(d1$amino_acid)
   d1 <- subset(d1, productive_frequency != "null")
   d1 <- subset(d1, length==15) #MODIFY HERE
   
   #prepare matrix
   mat <- sapply(as.character(d1$amino_acid), function(x){
      x <- as.character(x);
      x <- strsplit(x, "");
      x <- unlist( x )
   })
   mat <- t(mat) #cdr3 sequence matirx: 15 col matrix
      
   #mutual information
   res <- data.frame()
   for( i in 1:15 ){
   for( k in 1:15 ){
      print(paste0( i, ":", k))
      d1 <- mat[, c(i,k)]
      d1 <- as.data.frame(d1)
      colnames(d1) <- c("X","Y")
      d1 <- d1[ d1$X!="NA" &  d1$Y!="NA",] #revemo NA
      
      xlist <- as.character( unique(d1$X) )
      ylist <- as.character( unique(d1$Y) )
      
      #make frequency table
      freq <- matrix(0, nrow=length(xlist), ncol=length(ylist))
      row.names(freq) <- xlist
      colnames(freq) <- ylist
      
      for( x in xlist ){
      for( y in ylist ){
         both <- subset(d1, X==x & Y==y)
         freq[ x, y ] <- nrow( both )
      }
      }
      
      MI <- mi.empirical(freq, unit=c("log"))
      Hx <- entropy.empirical( colSums(freq) )
      Hy <- entropy.empirical( rowSums(freq) )
      NMI <- 2 * MI / ( Hx + Hy ) #normalized mutual entropy
      
      dump <- data.frame(i, k, MI, Hx, Hy, NMI)
      res <- rbind(res, dump)
      
   }
   }
   
   gz1 <- gzfile(ofile, "w")
   write.table(res, gz1, sep = "\t", quote = F, row.names = FALSE, col.names=TRUE)
   close(gz1)
   
}



