# Overview
We tested associations between HLA genotypes and TCR-CDR3 amino acid compositions. We treated the amino acid composition of CDR3 as a quantitative trait, and tested its association with HLA genotype; we call this CDR3 quantitative trait loci analysis (cdr3-QTL). This webpage contains the codes and the summary statistics from cdr3-QTL analysis. 
Reference: the manuscript is in submission (doi: in progress).

# Main dataset
We utilized publicly available TCR dataset. All raw TCR sequence data and genotype data of the discovery dataset and the replication dataset are available at Adaptive Biotechnologies immuneACCESS site: https://clients.adaptivebiotech.com/pub/emerson-2017-natgen
- Reference: Emerson, R. O. et al. Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA-mediated effects on the T cell repertoire. Nat. Genet. 1–10 (2017). doi:10.1038/ng.3822

## Results 1: the main results from cdr3-QTL analysis (./summary_stats/)
We utilized two different linear models as explained in this figure.
This is a schematic explanation of the linear model used in this study.
(a) Our strategy to calculate amino acid frequencies which we utilized for the main analysis. In this example, alanine (A) usage ratio at CDR3 position 110 is calculated for each individual.
(b) Two main linear models utilized in this study. In a multivariate multiple linear regression, a vector of frequency of 20 amino acids at a given position of CDR3 is the response variable; all amino acid alleles except one at a site of HLA are the explanatory variables. In a linear regression model, the frequency of
a single amino acid at a position of CDR3 is the response variable; a single amino acid allele at a site of HLA is the explanatory variable.

![image](./figure/Fig1_1.png)

## Results 2: CDR3 risk score
We developed a scoring system that quantifies the enrichment of these patterns in a given CDR3 sequence (we refer to this as the CDR3 risk score).
This is a schematic explanation of our strategy to calculate CDR3 risk score. The table shows the effect size estimate of cdr3-QTL analysis based on HLA risk scores (explained in our manuscript). Effect sizes for corresponding amino acid positions are summed up to calculate the CDR3 risk score. Effect sizes which passed a P value threshold were utilized. In the final analysis, we utlized Bonferroni corrected P < 0.05 as a threshold.

![image](./figure/Fig2.png)

