args   <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line
length <- as.character( args[1] )
pos <- as.character( args[2] )
    # length="15"; pos="9"

library(data.table)

 #output file and directory
ofile <- paste0("assoc/HIP/ra_score/length_v1/L", length,"_P", pos ,".lm.v1.txt.gz")

 #phenotype
pdata <- fread(paste0("gunzip -c phenotype.freq3/L",length,"_P",pos,"/single_aa_usage_ratio.txt.gz"),
   header=T, stringsAsFactors=F,)
pdata <- as.data.frame(pdata)
colnames(pdata)[2] <- "pname"
pdata <- pdata[ ,c("Sample", "pname", "rate")]

 #genotype (NEW GENOTYPE)
load("RA_haplotype_score_v1.RData") # rascore_mat: N x haplotypee
load("RA_haplotype_v1.RData")
betavec <- log(rascore_info$OR)
names(betavec) <- rascore_info$haplotype
betavec <- betavec[ colnames(rascore_mat) ]
genomat <- as.matrix(rascore_mat) %*% betavec # haplotype geenotype x gwas_beta of haplotype
geno <- data.frame(Sample=row.names(genomat), score=genomat[,1] )

 #covariate (only use genotype pca)
cdata <- fread("genotype_qc/HLA_A_B_C_DQB1_DRB1_4digit.pca", header=T, stringsAsFactors=F)
cdata <- as.data.frame(cdata)
row.names(cdata) <- cdata$id
colnames(cdata)[1] <- "Sample"
cdata <- cdata[,c("Sample", "PC1", "PC2", "PC3")]

 #main calculations
pnamelist <- read.table(paste0("phenotype.freq3/L",length,"_P",pos,"/single_aa_usage_ratio.pass"), stringsAsFactors=F)[,1]

res <- data.frame()

for( ptarget in pnamelist ){
   #prep data
   # ptarget="E"
   print( ptarget)
   
   pdata2 <- subset(pdata, pname==ptarget )
   pdata2 <- pdata2[,c("Sample", "rate")]
   M <- merge(pdata2, geno, by="Sample")
   M <- merge(M, cdata, by="Sample")
   x <- M$rate
   M$normrate <- qnorm( (rank(x, na.last="keep") - 0.5) / sum(!is.na(x)) ) #inverse normal normalization
   
   Nall <- nrow(M)
   phenomean <- mean( M$rate )
   
   #non-normalized data (flip y and x)
   test <- lm( score ~ rate + PC1 + PC2 + PC3 , data=M )
   out <- summary(test)$coefficients[ -1 ,c(1,2,4)] #remove intercept
   out <- data.frame(id=row.names(out), out)
   colnames(out) <- c("id", "beta", "se", "pvalue")
   out$type <- "org_scale"
   out$ptarget <- ptarget
   res <- rbind(res, out)
   
   #normalized data
   test <- lm( normrate ~ score + PC1 + PC2 + PC3 , data=M )
   out <- summary(test)$coefficients[ -1 ,c(1,2,4)] #remove intercept
   out <- data.frame(id=row.names(out), out)
   colnames(out) <- c("id", "beta", "se", "pvalue")
   out$type <- "norm_scale"
   out$ptarget <- ptarget
   res <- rbind(res, out)
   
}

gz1 <- gzfile(ofile, "w")
write.table(res, gz1, sep = "\t", quote = F, row.names = FALSE, col.names=TRUE)
close(gz1)



