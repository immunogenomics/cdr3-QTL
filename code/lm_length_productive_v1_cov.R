args   <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line
HLA <- as.character( args[1] )
length <- as.character( args[2] )
pos <- as.character( args[3] )
    # HLA="HLA_DRB1"; length="15"; pos="9"

library(data.table)

 #output file and directory
ofile <- paste0("assoc/HIP/lm/length_productive_v1_cov/", HLA, ".L", length,"_P", pos ,".lm.v2.txt.gz")

 #phenotype
pdata <- read.table(paste0("phenotype.freq3/L",length,"_P",pos,"/single_aa_usage_ratio.txt.gz"), header=T, stringsAsFactors=F)
pdata <- as.data.frame(pdata)
colnames(pdata)[2] <- "pname"
pdata <- pdata[ ,c("Sample", "pname", "rate")]

 #genotype (NEW GENOTYPE): 1% maf qc
gdata <- fread(paste0("genotype_qc/", HLA, ".imgt.aa.postqc.matrix"), header=T, stringsAsFactors=F)
gdata <- as.data.frame(gdata)
row.names(gdata) <- gdata$id
colnames(gdata)[1] <- "Sample"

 #covariate (only use genotype pca)
cdata <- read.table("genotype_qc/HLA_A_B_C_DQB1_DRB1_4digit.pca", header=T, stringsAsFactors=F)
row.names(cdata) <- cdata$id
colnames(cdata)[1] <- "Sample"
cdata <- cdata[,c("Sample", "PC1", "PC2", "PC3")]

 #main calculations
pnamelist <- read.table(paste0("phenotype.freq3/L",length,"_P",pos,"/single_aa_usage_ratio.pass"), stringsAsFactors=F)[,1]
gnamelist <- colnames(gdata)[ -1 ]

res <- data.frame()
for( ptarget in pnamelist ){
for( gtarget in gnamelist ){
        # ptarget = "E"; gtarget="AA_DRB1_11_32552136_VL"
   #prep data
   pdata2 <- subset(pdata, pname==ptarget )
   pdata2 <- pdata2[,c("Sample", "rate")]
   gdata2 <- gdata[ ,c("Sample", gtarget)]
   colnames(gdata2) <- c("Sample", "dose")
   M <- merge(gdata2, pdata2, by="Sample")
   M <- merge(M, cdata, by="Sample")
   x <- M$rate
   M$normrate <- qnorm( (rank(x, na.last="keep") - 0.5) / sum(!is.na(x)) ) #inverse normal normalization
   
   if( var(M$dose) > 0 & var(M$rate) ){
      #non-normalized ratio
      test1 <- lm( rate ~  dose + PC1 + PC2 + PC3, data=M)
      beta1 <- summary(test1)$coefficients["dose","Estimate"]
      pvalue1 <- summary(test1)$coefficients["dose","Pr(>|t|)"]
      
      #normalized ratio
      test2 <- lm( normrate ~  dose + PC1 + PC2 + PC3, data=M)
      beta2 <- summary(test2)$coefficients["dose","Estimate"]
      pvalue2 <- summary(test2)$coefficients["dose","Pr(>|t|)"]
      
      dump <- data.frame( pname=ptarget, gname=gtarget, beta1, pvalue1, beta2, pvalue2)
      res <- rbind(res, dump)
      
   } else {
      dump <- data.frame( pname=ptarget, gname=gtarget, beta1=NA, pvalue1=NA, beta2=NA, pvalue2=NA)
      res <- rbind(res, dump)
   }
}
}

res <- res[order(res$pvalue2), ]

gz1 <- gzfile(ofile, "w")
write.table(res, gz1, sep = "\t", quote = F, row.names = FALSE, col.names=TRUE)
close(gz1)



