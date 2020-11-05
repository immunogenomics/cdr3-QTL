args   <- commandArgs(trailingOnly = T)   # enable reading factors from linux command line
length <- as.character( args[1] )
pos <- as.character( args[2] )
    #  length="15"; pos="7"

library(data.table)
library(lme4)

 #read info
info <- read.table("assoc/HIP/lmm/v3_length_cov/lead_info.txt",header=T,stringsAsFactors=F)
info <- info[ info$length==length & info$pos==pos, ]

 #output file and directory
ofile <- paste0("assoc/HIP/lmm/v3_length_vfixed_cov/L", length, ".", pos, ".lmm.v3_length.txt.gz")

 #phenotype
pdata <- read.table(paste0("phenotype.freq3/L", length, "_P", pos, "/vgene_cond_single_aa_usage_ratio.txt.gz"),
   header=T, stringsAsFactors=F)
pass <- read.table( paste0("phenotype.freq3/L", length, "_P", pos, "/single_aa_usage_ratio.pass"), stringsAsFactors=F)[,1]
pdata <- pdata[is.element(pdata$AA, pass), ]

 #covariate (only use genotype pca)
cdata <- read.table("genotype_qc/HLA_A_B_C_DQB1_DRB1_4digit.pca", header=T, stringsAsFactors=F)
row.names(cdata) <- cdata$id
colnames(cdata)[1] <- "Sample"
cdata <- cdata[,c("Sample", "PC1", "PC2", "PC3")]

res <- data.frame()

for( i in 1:nrow(info) ){
      # i = 1
   print(i)
   
   ptarget = info[i, "pname"]
   gtarget = info[i, "gname"]
   beta1_lm = info[i, "beta1"]
   pvalue1_lm = info[i, "pvalue1"]
   beta2_lm = info[i, "beta2"]
   pvalue2_lm = info[i, "pvalue2"]
   HLA = info[i, "gdata"]
   
   #genotype (NEW GENOTYPE)
   gdata <- read.table(paste0("genotype_qc/HLA_", HLA, ".imgt.aa.postqc.matrix"), header=T, stringsAsFactors=F)
   row.names(gdata) <- gdata$id
   colnames(gdata)[1] <- "Sample"
   
   #prep data
   pdata2 <- subset(pdata, AA==ptarget & obs > 1000 )
   pdata2 <- pdata2[,c("Sample", "rate", "Vgene")]
   
   TB <- table(pdata2$Vgene)
   keepv <- names( TB[ TB > ( length(unique(pdata2$Sample)) /2  ) ] )
   
   pdata2 <- pdata2[ pdata2$Vgene %in% keepv, ]
   
   gdata2 <- gdata[ ,c("Sample", gtarget)]
   colnames(gdata2) <- c("Sample", "dose")
   
   M <- merge(gdata2, pdata2, by="Sample")
   M <- merge(M, cdata, by="Sample")
   x <- M$rate
   M$normrate <- qnorm( (rank(x, na.last="keep") - 0.5) / sum(!is.na(x)) ) #inverse normal normalization
   
   if( var(M$rate) >0 & var(M$dose) > 0 & length(unique(M$Vgene)) > 1 & length(unique(M$Sample)) > 1 ){
      #non-normalized ratio
      test <- lmer( rate ~ dose + PC1 + PC2 + PC3 + Vgene + (1|Sample), data=M)
      beta1 <- summary(test)$coefficients["dose","Estimate"]
      se1 <- summary(test)$coefficients["dose","Std. Error"]
      
      H0 <- lmer( rate ~             PC1 + PC2 + PC3 + Vgene + (1|Sample), data=M, REML = F )
      H1 <- lmer( rate ~ dose + PC1 + PC2 + PC3 + Vgene + (1|Sample), data=M, REML = F )
      ANNO <- anova(H0, H1)
      lrpval1 <- ANNO[["Pr(>Chisq)"]][2]
      
      #normalized ratio
      test <- lmer( normrate ~ dose + PC1 + PC2 + PC3 + Vgene + (1|Sample), data=M)
      beta2 <- summary(test)$coefficients["dose","Estimate"]
      se2 <- summary(test)$coefficients["dose","Std. Error"]
      
      H0 <- lmer( normrate ~             PC1 + PC2 + PC3 +  Vgene + (1|Sample), data=M, REML = F )
      H1 <- lmer( normrate ~ dose + PC1 + PC2 + PC3 +  Vgene + (1|Sample), data=M, REML = F )
      ANNO <- anova(H0, H1)
      lrpval2 <- ANNO[["Pr(>Chisq)"]][2]
      
      dump <- data.frame( pname=ptarget, gname=gtarget, 
         beta1_lm, pvalue1_lm,
         beta2_lm, pvalue2_lm, HLA,
         beta1, se1, lrpval1,
         beta2, se2, lrpval2
         )
      
      res <- rbind(res, dump)
      
   }
}

gz1 <- gzfile(ofile, "w")
write.table(res, gz1, sep = "\t", quote = F, row.names = FALSE, col.names=TRUE)
close(gz1)



