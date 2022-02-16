
# Buitrago et al., 2022 16S rRNA amplicon data analysis

This directory contains the scripts used to produce figures and run analyses for the 16S rRNA data.

## Workflow

### inference of ASV, QCs and plotting overall bacterial diversity 

1. Amplicon Sequence variance (ASV) were inferred using [dada2](https://github.com/benjjneb/dada2) using the script `Spis_Pver_dada2.R`
2. Quality checks (i.g., removal of putatively contaminant ASVs and removal of samples with < 1000 reads) were done   using the script `Spis_Pver_QC.R`
3. barplots of most abundant bacterial genera and families were done using the script `Spis_Pver_barplots.R`

### Beta diversity plotting and stats
4. Ordination plots were done using the script `Spis_Pver_ordination.R`
5. Beta dispersion was evaluated using the script `Spis_Pver_BetaDispersion.R`
6. PERMANOVAs were done using the script `Spis_Pver_permanova.R`

### Alpha diversity plotting and stats
7. Alpha diversity estimates and statistical comparisons were done using the script `Spis_Pver_alpha.R`
8. Indicator bacterial ASVs across host population clusters and temperature clusters were identified using the script  `Spis_Pver_indicSpecies.R`
