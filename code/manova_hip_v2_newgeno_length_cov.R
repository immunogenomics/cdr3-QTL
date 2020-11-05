args   <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line
length <- as.character( args[ 1 ] )
pos <- as.character( args[ 2 ] )
    # HLA="HLA_DRB1"; length="14"; pos="8"

.libPaths( "~/bin/R" )
library(data.table)
library(MVLM)

for( HLA in c("HLA_A","HLA_B","HLA_C","HLA_DPA1","HLA_DPB1","HLA_DQA1","HLA_DQB1","HLA_DRB1") ) {

 #output file and directory
ofile <- paste0("assoc/HIP/manova/v2_newgeno_length_cov/", HLA, ".", length, ".", pos,".manova.v2_newgeno.txt.gz")

 #phenotype
pdata <- fread(paste0("gunzip -c phenotype.freq3/L",length,"_P",pos,"/single_aa_usage_ratio.txt.gz"), 
   header=T, stringsAsFactors=F,)
pdata <- as.data.frame(pdata)
pass <- pnamelist <- read.table(paste0("phenotype.freq3/L",length,"_P",pos,"/single_aa_usage_ratio.pass"),
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

 #genotype
 #gdata <- fread(paste0("genotype_qc/", HLA, ".imgt.aa.postqc.maf5.matrix"), header=T, stringsAsFactors=F)
gdata <- fread(paste0("genotype_qc/", HLA, ".imgt.aa.postqc.matrix"), header=T, stringsAsFactors=F)
gdata <- as.data.frame(gdata)
row.names(gdata) <- gdata$id
colnames(gdata)[1] <- "Sample"

 #covariate (only use genotype pca)
cdata <- fread("genotype_qc/HLA_A_B_C_DQB1_DRB1_4digit.pca", header=T, stringsAsFactors=F)
cdata <- as.data.frame(cdata)
row.names(cdata) <- cdata$id
colnames(cdata)[1] <- "Sample"
cdata <- cdata[,c("Sample", "PC1", "PC2", "PC3")]

 #main calculations
gtargettaglist <- unique( sapply(strsplit( colnames(gdata)[-1], "_"), 
   function(x){paste0(x[1],"_",x[2],"_",x[3],"_",x[4])} ))
pnamelist <- gsub("-",".",pnamelist)

res <- data.frame()
for( gtargettag in gtargettaglist ){
    # gtargettag="AA_DRB1_13_32552130"
   gnamelist <- grep(gtargettag, colnames(gdata),value=T)
   
   if(  length(gnamelist)==1 ){
      #bi-allelic
      gnamelist <- gnamelist
   } else {
      gnamelistlen <- sapply(gnamelist,function(x){
         x <- unlist(strsplit(x,"_"));
         aa <- x[5];
         nchar(aa) })
      gnamelist <- gnamelist[ gnamelistlen==1 ]
      if(  length(gnamelist)==1 ){
         gnamelist <- gnamelist
      } else {
         gnamelist <- gnamelist[-1] #use one aa as reference
      }
   }
   
   gdata2 <- gdata[ ,c("Sample", gnamelist)]
   colnames(gdata2) <- c("Sample", paste0("dose", 1:length(gnamelist)))
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
   #no pc cov
   x <- paste(pnamelist, collapse=",")
   x <- paste("cbind(",x,") ~")
   y <- paste(colnames(gdata2)[-1], collapse=" + ")
   formula1 <- as.formula(paste(x,y))
   mvlm.res <- mvlm( formula1 , data = M)
      #summary(mvlm.res)
   pvalue_omni <- mvlm.res$pv["Omnibus Effect", 1]
   pseudo_rsq <- mvlm.res$pseudo.rsq["Omnibus Effect",1]
   pseudo_rsq2 <- sum( mvlm.res$pseudo.rsq[,1] ) - pseudo_rsq #only dose term (without intercept term)
   
   #full
   x <- paste(pnamelist, collapse=",")
   x <- paste("cbind(",x,") ~")
   y <- paste(colnames(gdata2)[-1], collapse=" + ")
   y <- paste(y,"+ PC1 + PC2 + PC3")
   formula1 <- as.formula(paste(x,y))
   mvlm.res <- mvlm( formula1 , data = M)
      # summary(mvlm.res)
   pseudo_rsq_full <- mvlm.res$pseudo.rsq["Omnibus Effect",1] 
   pseudo_rsq_pc <- sum( mvlm.res$pseudo.rsq[c("PC1","PC2","PC3"),1]  ) #pc
   pseudo_rsq_dose <- sum( mvlm.res$pseudo.rsq[,1] ) - pseudo_rsq_full - pseudo_rsq_pc
   
   #null
   x <- paste(pnamelist, collapse=",")
   x <- paste("cbind(",x,") ~")
   y <- paste(" PC1 + PC2 + PC3")
   formula0 <- as.formula(paste(x,y))
   mvlm.res <- mvlm( formula0 , data = M)
      # summary(mvlm.res)
   pseudo_rsq_null <- mvlm.res$pseudo.rsq["Omnibus Effect",1] 
   
   dump <- data.frame(HLA, length, pos, gtargettag, pvalue_lm, pvalue_omni, 
      pseudo_rsq, pseudo_rsq2, pseudo_rsq_full, pseudo_rsq_pc, pseudo_rsq_dose, pseudo_rsq_null)
   
     #dump <- data.frame(HLA, length, pos, gtargettag, pvalue_lm, pvalue_omni, pseudo_rsq)
   
   res <- rbind(res, dump)
}

gz1 <- gzfile(ofile, "w")
write.table(res, gz1, sep = "\t", quote = F, row.names = FALSE, col.names=TRUE)
close(gz1)

}




