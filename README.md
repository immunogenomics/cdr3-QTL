# Overview
We tested associations between HLA genotypes and TCR-CDR3 amino acid compositions. We treated the amino acid composition of CDR3 as a quantitative trait, and tested its association with HLA genotype; we call this CDR3 quantitative trait loci analysis (cdr3-QTL). This webpage contains the codes and the summary statistics from cdr3-QTL analysis. 
- Reference: our manuscript is under submission (doi: in progress).

# Main dataset
We utilized publicly available TCR dataset. All raw TCR sequence data and genotype data of the main dataset are available at Adaptive Biotechnologies immuneACCESS site: https://clients.adaptivebiotech.com/pub/emerson-2017-natgen
- Reference: Emerson, R. O. et al. Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA-mediated effects on the T cell repertoire. Nat. Genet. 1–10 (2017). doi:10.1038/ng.3822

## Results 1: the main results from cdr3-QTL analysis (./summary_stats/)
- Our strategy to prepare quantitative phenotype of CDR3. We utilized amino acid frequencies (panel a). In this example, alanine (A) usage ratio at CDR3 position 110 is calculated for each individual.
- We utilized two different linear models (panel b).
- (i) multivariate multiple linear regression: a vector of frequency of 20 amino acids at a given position of CDR3 is the response variable; all amino acid alleles except one at a site of HLA are the explanatory variables. 
- (ii) linear regression model: the frequency of a single amino acid at a position of CDR3 is the response variable; a single amino acid allele at a site of HLA is the explanatory variable.

![image](./figure/Fig1_1.png)

## Results 2: CDR3 risk score
From results 1, we found many CDR3 amino acid patterns associated wit HLA risk of autoimmune diseases. We developed a scoring system that quantifies the enrichment of these patterns in a given CDR3 sequence (we refer to this as the CDR3 risk score). This is a schematic explanation of our strategy to calculate CDR3 risk score. The table shows the effect size estimate of cdr3-QTL analysis based on HLA risk scores (explained in our manuscript). Effect sizes for corresponding amino acids are summed when such amino acids exist in a target CDR3 sequence. This sum is defined as the CDR3 risk score. Effect sizes which passed a P value threshold were utilized. Although the threshold of P < 0.05 is used in this figure, we utlized the Bonferroni corrected P < 0.05 in the final analysis.

![image](./figure/Fig2.png)
