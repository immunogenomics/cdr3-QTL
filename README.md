# Overview
We tested associations between HLA genotypes and TCR-CDR3 amino acid compositions. We treated the amino acid composition of CDR3 as a quantitative trait, and tested its association with HLA genotypes; we call this CDR3 quantitative trait loci analysis (cdr3-QTL). This webpage contains the codes and the summary statistics from our cdr3-QTL analysis. 
- Reference: Ishigaki K et al. HLA autoimmune risk alleles restrict the hypervariable region of T cell receptors. MedRxiv 2020 (https://www.medrxiv.org/content/10.1101/2020.11.08.20227983v1

## Main dataset
We analyzed publicly available TCR dataset. All raw TCR sequence data and genotype data of the main dataset are available at Adaptive Biotechnologies immuneACCESS site: https://clients.adaptivebiotech.com/pub/emerson-2017-natgen
- Reference: Emerson, R. O. et al. Immunosequencing identifies signatures of cytomegalovirus exposure history and HLA-mediated effects on the T cell repertoire. Nat. Genet. 1–10 (2017). doi:10.1038/ng.3822

## Why we conducted this study?
- Polymorphisms in the human leukocyte antigen (HLA) genes within the major histocompatibility complex (MHC) locus strongly influence autoimmune disease risk.
- Two non-exclusive hypotheses exist about the pathogenic role of HLA alleles:
- i) the central hypothesis, where HLA risk alleles influence thymic selection so that the probability of T cell receptors (TCRs) reactive to pathogenic antigens is increased
- ii) the peripheral hypothesis, where HLA risk alleles increase the affinity for pathogenic antigens. 
- The peripheral hypothesis has been the main research focus in autoimmunity, while human data on the central hypothesis are lacking.
- Therefore, to explore genetic evidence supporting the central hypothesis, we investigated the influence of HLA alleles on TCR composition at the highly diverse complementarity determining region 3 (CDR3), where TCR recognizes antigens.

## How we generated quantitative phenotypes from TCR sequence data?
- Our strategy to prepare quantitative phenotype of CDR3. We utilized amino acid frequencies (panel a). In this example, alanine (A) usage ratio at CDR3 position 110 is calculated for each individual.

## Main models in our analysis (panel b)
- *We prepared a vignette and example input datasets to explan our models (./vignettes/example_cdr3QTL_v2.ipynb and ./vignettes/DRB1_site13_L13CDR3_p109.RData)*
- (i) multivariate multiple linear regression: a vector of frequency of 20 amino acids at a given position of CDR3 is the response variable; all amino acid alleles except one at a site of HLA are the explanatory variables.
- (ii) linear regression model: the frequency of a single amino acid at a position of CDR3 is the response variable; a single amino acid allele at a site of HLA is the explanatory variable.

![image](./figure/Fig1_1.png)

## Results from our cdr3-QTL analysis

### Results 1: the main results from cdr3-QTL analysis (./summary_stats/). 
- When a file name includes "rm_gl", it indicates that the results were based on CDR3 phenotypes excluding germline-encoded sequences.

### Results 2: CDR3 risk score
- From results 1, we found many CDR3 amino acid patterns associated wit HLA risk of autoimmune diseases. We developed a scoring system that quantifies the enrichment of these patterns in a given CDR3 sequence (we refer to this as the CDR3 risk score). 
- This is a schematic explanation of our strategy to calculate CDR3 risk score. The table shows the effect size estimate of cdr3-QTL analysis based on HLA risk scores (explained in our manuscript). Effect sizes for corresponding amino acids are summed when such amino acids exist in a target CDR3 sequence. This sum is defined as the CDR3 risk score. Effect sizes which passed a P value threshold were utilized. Although the threshold of P < 0.05 is used in this example figure, we used the Bonferroni corrected P < 0.05 in our final analysis.

![image](./figure/Fig2.png)

