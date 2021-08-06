## The code used in this project.

### Script of genotype QC and PCA
- ./genotype/genotype_4digit_qc_pca.ipynb


### Scripts for cdr3-QTL analysis using HLA amino acid genotype (multivariate multiple regression)
- ./stat_test/manova_hip_v2_newgeno_length_cov_new_rm_gl.R


### Scripts for cdr3-QTL analysis using HLA amino acid genotype (linear regression)
- ./stat_test/lm_length_productive_v1_cov_rm_gl.R


### Scripts for cdr3-QTL analysis using HLA amino acid genotype conditioning on V/J usage (linear mixed model)
- ./stat_test/lmm_hip_v3_length_vfixed_cov_rm_gl.R
- ./stat_test/lmm_hip_v3_length_jfixed_cov_rm_gl.R


### Scripts for cdr3-QTL analysis using HLA risk score (linear regression)
- ./stat_test/celiac_score_length_v1_rm_gl.R
- ./stat_test/ra_score_length_v1_rm_gl.R
- ./stat_test/t1d_score_length_v1_rm_gl.R


### Scripts for CDR3 data structure (entropy and mutual information)
- ./tcr_covariance/cdr3_aa_pheno_cor_v3.R
- ./tcr_covariance/cdr3_aa_pheno_cor_vgene_v3.R
- ./tcr_covariance/cdr3_aa_vjgene_var_v3.R
- ./tcr_covariance/phenotype_pro_ratio.sh


### Scripts for generating CDR3 phenotypes (=amino acid usage frequency)
- ./tcr_phenotype/cdr3_aaratio_length_rm_germline.sh: generating productive CDR3 phenotypes (used for the main analysis)
- ./tcr_phenotype/cdr3_aaratio_length_ce_rm_gl.sh: generating productive CDR3 phenotypes including clonal expansion effects
- ./tcr_phenotype/cdr3_aaratio_length_rm_germline_vgene.sh: generating productive CDR3 phenotypes conditioning on V gene usage
- ./tcr_phenotype/cdr3_aaratio_length_rm_germline_jgene.sh: generating productive CDR3 phenotypes conditioning on J gene usage
- ./tcr_phenotype/cdr3_aaratio_length_nopub_rm_gl.sh: generating productive CDR3 phenotypes excluding public clonotypes
- ./tcr_phenotype/cdr3_aaratio_nonpro_length.sh: generating non-productive CDR3 phenotypes
- ./tcr_phenotype/phenotype_pro_ratio_downsample.sh: generating downsampled productive CDR3 phenotypes


### Scripts for MIXCR (used to infer non-productive CDR3 sequence)
- ./tcr_phenotype/mixcr_tcr_hip.sh


### Scripts for removing public clonotypes
- ./tcr_phenotype/remove_include_public_clonotype.sh


