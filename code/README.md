## The code used in this project.

### Explanation about multivariate multple linear regressino
- manova_cdr3QTL_demo_v1.pdf: explains how P values (based on MANOVA model comparison) and variance explained were calculated in this project.
- example_data_HLA_DRB1_site13_L13P109.RData: data used in this explanation.

### Scripts for CDR3 data structure (entropy and mutual information)
- cdr3_aa_pheno_cor_v3.R
- cdr3_aa_pheno_cor_vgene_v3.R
- cdr3_aa_vjgene_var_v3.R
- phenotype_pro_ratio.sh

### Scripts for generating CDR3 phenotypes (=amino acid usage frequency)
- cdr3_aaratio_length.sh: generating productive CDR3 phenotypes (used for the main analysis)
- cdr3_aaratio_length_ce.sh: generating productive CDR3 phenotypes including clonal expansion effects
- cdr3_aaratio_length_cond_vgene.sh: generating productive CDR3 phenotypes conditioning on V gene usage
- cdr3_aaratio_length_cond_jgene.sh: generating productive CDR3 phenotypes conditioning on J gene usage
- cdr3_aaratio_length_nopub.sh: generating productive CDR3 phenotypes excluding public clonotypes
- cdr3_aaratio_nonpro_length.sh: generating non-productive CDR3 phenotypes
- phenotype_pro_ratio_downsample.sh: generating downsampled productive CDR3 phenotypes

### Scripts for cdr3-QTL analysis using HLA amino acid genotype (linear regression)
- lm_length_productive_v1_cov.R

### Scripts for cdr3-QTL analysis using HLA amino acid genotype conditioning on V/J usage (linear mixed model)
- lmm_hip_v3_length_vfixed_cov.R
- lmm_hip_v3_length_jfixed_cov.R

### Scripts for cdr3-QTL analysis using HLA amino acid genotype (multivariate multiple regression)
- manova_hip_v2_newgeno_length_cov.R

### Scripts for cdr3-QTL analysis using HLA risk score (linear regression)
- celiac_score_length_v1.R
- ra_score_length_v1.R
- t1d_score_length_v1.R

### Scripts for MIXCR (used to infer non-productive CDR3 sequence)
- mixcr_tcr_hip.sh

### Scripts for removing public clonotypes
- remove_public_clonotype.sh
