# Buitrago et al., 2022 RAD-Seq data analysis

This directory contains the scripts used to run analyses and produce figures for the RAD-Seq data.

## Workflow
1. Processing of raw sequencing data (00_preprocessing.sh)

### Stylophora pistillata
2. Reference genome RADseq genotyping and filtering of VCF files (01_Spis_RADfiltering.sh)
3. Preliminary relatedness analysis using Identity-By-Descent and visualization (02_Spis_SNPRelate.R)
4. Analysis of Molecular Variance - AMOVA (03_Spis_AMOVA.R)
5. Pairwise FST analysis and visualization (04_Spis_PairwiseFST.R)
6. Isolation by distance analysis and visualization (05_Spis_IBD.R)
7. Principal component analysis and visualization (06_Spis_PCAstructure.R)
8. Admixture analysis and visualization (07a_Spis_ADMIXTURE.sh; 07b_Spis_ADMIXTUREvisualization.R)
9. Linkage Disequilibrium (LD) analysis and visualization of LD decay (08a_Spis_LDanalysis.sh; 08b_Spis_LDvisualization.R)
10. Candidate SNPs for positive selection analyses and visualization (09a_Spis_CandidateSNPs_FormatingFiles.R; 09b_Spis_CandidateSNPs_BAYPASS.R; 09c_Spis_CandidateSNPs_BAYESCAN.R; 09d_Spis_CandidateSNPs_visualization.R)

### Pocillopora verrucosa
11. Reference genome RADseq genotyping and filtering of VCF files (01_Pver_RADfiltering.sh)
12. Preliminary relatedness analysis using Identity-By-Descent and visualization (02_Pver_SNPRelate.R)
13. Analysis of Molecular Variance - AMOVA (03_Pver_AMOVA.R)
14. Pairwise FST analysis and visualization (04_Pver_PairwiseFST.R)
15. Isolation by distance analysis and visualization (05_Pver_IBD.R)
16. Principal component analysis and visualization (06_Pver_PCAstructure.R)
17. Admixture analysis and visualization (07a_Pver_ADMIXTURE.sh; 07b_Pver_ADMIXTUREvisualization.R)
18. Linkage Disequilibrium (LD) analysis and visualization of LD decay (08a_Pver_LDanalysis.sh; 08b_Pver_LDvisualization.R)
19. Candidate SNPs for positive selection analyses and visualization (09a_Pver_CandidateSNPs_FormatingFiles.R; 09b_Pver_CandidateSNPs_BAYPASS.R; 09c_Pver_CandidateSNPs_BAYESCAN.R; 09d_Pver_CandidateSNPs_visualization.R)

