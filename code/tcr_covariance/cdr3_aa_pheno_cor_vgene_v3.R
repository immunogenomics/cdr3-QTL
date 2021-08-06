 # length = 15
args <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line
dd="/data/srlab/kishigaki/data/TCR_robins"
library(data.table)
library(entropy)
LIST=read.table("Info/genotypeid_hiptcrid.txt", stringsAsFactors=F)[,2]

df<-data.frame()

for(Sample in LIST){
   print(Sample)
   ofile <- paste0("phenotype.summary/HIP/aa_pos_cor_vgene3/",Sample,".tcr_each_aapos_cor_vj.txt.gz")
   
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
   
   #mutual information (cdr3aa vs V or J )
   res <- data.frame()
   for( i in 1:15 ){
      x1 <- mat[, i, drop=F]
      x2 <- d1[,c("v_gene", "j_gene")]
      x3 <- cbind(x1, x2)
      x3 <- as.data.frame(x3)
      colnames(x3) <- c("X","V","J")
      x3 <- x3[ x3$X!="NA", ] #revemo NA
      
      #v gene vs aa
      xlist <- as.character( unique(x3$X) )
      vlist <- as.character( unique(x3$V) )
      freq <- matrix(0, nrow=length(xlist), ncol=length(vlist))
      row.names(freq) <- xlist
      colnames(freq) <- vlist
      for( x in xlist ){
      for( v in vlist ){
         both <- subset(x3, X==x & V==v)
         freq[ x, v ] <- nrow( both )
      }
      }
      
      MI <- mi.empirical(freq, unit=c("log"))
      Hx <- entropy.empirical( colSums(freq) )
      Hy <- entropy.empirical( rowSums(freq) )
      NMI <- 2 * MI / ( Hx + Hy ) #normalized mutual entropy
      dump <- data.frame(i, k=0, MI, Hx, Hy, NMI)
      res  <- rbind(res , dump)
      
      #j gene vs aa
      xlist <- as.character( unique(x3$X) )
      jlist <- as.character( unique(x3$J) )
      freq <- matrix(0, nrow=length(xlist), ncol=length(jlist))
      row.names(freq) <- xlist
      colnames(freq) <- jlist
      for( x in xlist ){
      for( j in jlist ){
         both <- subset(x3, X==x & J==j )
         freq[ x, j ] <- nrow( both )
      }
      }
      
      MI <- mi.empirical(freq, unit=c("log"))
      Hx <- entropy.empirical( colSums(freq) )
      Hy <- entropy.empirical( rowSums(freq) )
      NMI <- 2 * MI / ( Hx + Hy ) #normalized mutual entropy
      dump <- data.frame(i, k=16, MI, Hx, Hy, NMI)
      res  <- rbind(res , dump)
   }
   
   #VJ MI
   x2 <- d1[,c("v_gene", "j_gene")]
   x3 <- as.data.frame(x2)
   colnames(x3) <- c("V","J")

   vlist <- as.character( unique(x3$V) )
   jlist <- as.character( unique(x3$J) )
   freq <- matrix(0, nrow=length(vlist), ncol=length(jlist))
   row.names(freq) <- vlist
   colnames(freq) <- jlist
   for( v in vlist ){
   for( j in jlist ){
      both <- subset(x3, V==v & J==j )
      freq[ v, j ] <- nrow( both )
   }
   }
   
   MI <- mi.empirical(freq, unit=c("log"))
   Hx <- entropy.empirical( colSums(freq) )
   Hy <- entropy.empirical( rowSums(freq) )
   NMI <- 2 * MI / ( Hx + Hy ) #normalized mutual entropy
   dump <- data.frame(i=0, k=16, MI, Hx, Hy, NMI)
   res  <- rbind(res , dump)
   
   gz1 <- gzfile(ofile, "w")
   write.table(res, gz1, sep = "\t", quote = F, row.names = FALSE, col.names=TRUE)
   close(gz1)
   
}




