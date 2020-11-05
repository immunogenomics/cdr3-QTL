#demo of calculating CDR3 risk score

x <- c(
   "CASGKQQGEAF","CASSVGQGLGYTF","CSQAPSREERGEDTQYF",
   "CASSQDRVGGTDTQYF","CASSGSLTGTGGAEAF","CSASGTGEGYTF","CASSLRGFIQPNEQF","CASFLGPVFPGGYTF"
)  #simulated CDR3 sequences

#example of celiac disease score
load("HIP_celiac_score_aausage_normbetas_v1.RData")
pos_beta_all <- as.data.frame(out)

score <- sapply( x, function( x ){
   len_cdr3 <- nchar(x);
   if( len_cdr3 >=12 & len_cdr3 <= 18 ){
        df <- data.frame( ptarget=unlist(strsplit(x, "")), pos=1:len_cdr3 );
        df$tag <- paste0(df$ptarget, ":", df$pos);
       
        pos_beta <- subset(pos_beta_all, length==len_cdr3 ) #len_cdr3=15
        pos_beta$tag <- paste0(pos_beta$ptarget, ":", pos_beta$pos)
        pos_beta <- pos_beta[,c("tag","beta")]
        
        df <- merge( pos_beta, df, by="tag" );
        if ( nrow(df) == 0 ){
            score = 0
        } else {
            score = sum(df$beta)
        };
    } else {
        score = 0
    };
    return( score )
})



