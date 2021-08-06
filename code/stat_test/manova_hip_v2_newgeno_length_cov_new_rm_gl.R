#args   <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line
#length <- as.character( args[ 1 ] )
#pos <- as.character( args[ 2 ] )
    # HLA="HLA_DRB1"; length="13"; pos="6"

HLA="HLA_DRB1"; length="13"; pos="6"

.libPaths( "~/bin/R" )
library(data.table)
library(MVLM)

 #for( HLA in c("HLA_A","HLA_B","HLA_C","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1","HLA_DRB1") ) {

 #output file and directory
 #ofile <- paste0("assoc/HIP/manova/v2_newgeno_length_cov_newscript_rm_gl/", HLA, ".", length, ".", pos,".manova.v2_newgeno.txt.gz")

 #phenotype
pdata <- read.table(paste0("../../data/tcr_phenotype/L",length,"_P",pos,"_rm_gl/single_aa_usage_ratio.txt.gz"), 
   header=T, stringsAsFactors=F,)
pdata <- as.data.frame(pdata)
pass <- pnamelist <- read.table(paste0("../../data/tcr_phenotype/L",length,"_P",pos,"_rm_gl/single_aa_usage_ratio.pass"),
   stringsAsFactors=F)[,1]
pdata <- pdata[is.element(pdata$AA, pass), ]
pdata <- subset(pdata, AA!="unresolved")

pnamelist <- unique(pdata$AA)
pall <- data.frame()
for( ptarget in pnamelist ){
    pdata2 <- subset(pdata, AA==ptarget)
    pdata2 <- pdata2[,c("Sample", "rate")]
    x <- pdata2$rate
    pdata2$normrate <- qnorm( (rank(x, na.last="keep") - 0.5) / sum(!is.na(x)) ) #inverse normal normalization
    #head(M)
    #cancel V gene effect by averaging
    df <- pdata2[,c("Sample","normrate")]
    df$ptarget <- ptarget
    pall <- rbind(pall,df)
}

mat <- matrix(0, nrow=length(unique(pall$Sample)), ncol=length(unique(pall$ptarget)) )
row.names(mat) <- unique(pall$Sample)
colnames(mat) <- unique(pall$ptarget)

for( i in unique(pall$Sample)){
for( k in unique(pall$ptarget)){
   out <- subset(pall, Sample==i & ptarget==k)
   if( nrow(out) ==1 ){
      mat[i,k] <- out$normrate
   }else{
       mat[i,k] <- 0
   }
}
}

mat <- data.frame(Sample=row.names(mat),mat)

 #covariate (only use genotype pca)
cdata <- read.table("../../data/genotype/HLA_A_B_C_DQB1_DRB1_4digit.pca", header=T, stringsAsFactors=F)
cdata <- as.data.frame(cdata)
row.names(cdata) <- cdata$id
colnames(cdata)[1] <- "Sample"
cdata <- cdata[,c("Sample", "PC1", "PC2", "PC3")]

 #main calculations
gtargettaglist <- dir("../../data/genotype/v2")
gtargettaglist <- grep( paste0("AA_", unlist(strsplit(HLA,"_"))[2] ,"_"), gtargettaglist, value=T)
gtargettaglist <- gsub(".txt", "", gtargettaglist)

res <- data.frame()
for( gtargettag in gtargettaglist ){
    # gtargettag="AA_DRB1_13_32552130"
   
   gdata2 <- read.table(paste0("../../data/genotype/v2/",gtargettag,".txt"), header=T, stringsAsFactors=F)
      # ref allele has already removed (most common)
   M <- merge(gdata2, mat, by="Sample")
   M <- merge(M, cdata, by="Sample")
   
   #make formula
   x <- paste(pnamelist, collapse=",")
   x <- paste("cbind(",x,") ~")
   y <- paste(colnames(gdata2)[-1], collapse=" + ")
   y <- paste(y,"+ PC1 + PC2 + PC3")
   formula1 <- as.formula(paste(x,y))
   
   x <- paste(pnamelist, collapse=",")
   x <- paste("cbind(",x,") ~")
   y <- paste(" PC1 + PC2 + PC3")
   formula0 <- as.formula(paste(x,y))
   
   #pvalue
   mod1 <- lm( formula1 , data = M)
   mod0 <- lm( formula0, data = M)
   test <- anova(mod0, mod1)
   pvalue_lm <- test$"Pr(>F)"[2]
   
   #STEP2: use MVLM
   #full
   x <- paste(pnamelist, collapse=",")
   x <- paste("cbind(",x,") ~")
   y <- paste(colnames(gdata2)[-1], collapse=" + ")
   y <- paste(y,"+ PC1 + PC2 + PC3")
   formula1 <- as.formula(paste(x,y))
   mvlm.res <- mvlm( formula1 , data = M)
      # summary(mvlm.res)
   pseudo_rsq_full <- mvlm.res$pseudo.rsq["Omnibus Effect",1] 
   
   #null
   x <- paste(pnamelist, collapse=",")
   x <- paste("cbind(",x,") ~")
   y <- paste(" PC1 + PC2 + PC3")
   formula0 <- as.formula(paste(x,y))
   mvlm.res <- mvlm( formula0 , data = M)
      # summary(mvlm.res)
   pseudo_rsq_null <- mvlm.res$pseudo.rsq["Omnibus Effect",1] 
   
   dump <- data.frame(HLA, length, pos, gtargettag, pvalue_lm,  
      pseudo_rsq_full, pseudo_rsq_null)
   
   res <- rbind(res, dump)
}

res <- res[order(res$pvalue_lm),]

 #gz1 <- gzfile(ofile, "w")
 #write.table(res, gz1, sep = "\t", quote = F, row.names = FALSE, col.names=TRUE)
 #close(gz1)
 #}

show(head(res))


