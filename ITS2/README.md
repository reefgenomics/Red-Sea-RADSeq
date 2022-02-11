
# Buitrago et al., 2022 ITS2 analyses

This directory contains the scripts used to produce figures and run analyses for the ITS2 data.

Input files required for the scripts are detailed below. If the files are not present in this repo. please
check the 

# Bar figure plotting and general stats

The barplot ITS2 figure and general stats were generated using `./buitrago.py`.
This script was also used to make additional figures that were not used in the MS. As such, there is redundant code in the script.
The two classes relevant to the plotting and stats generation are `BuitragoBars_clustered_profiles` and `CalculateAverageProfDistances`.
Instances of both of these will be instantiated and the relevant outputs will be produced when the script is run.

Use the `./buitrago_env.yml` to generate a conda envronment containing the required dependencies.

The following files are required as input for this script:

- `./pver.ind.ordered.byclusters.txt`: list of the *P. verrucosa* samples to use in plotting and ordinations

- `./spis.ind.ordered.byclusters.txt`: list of the *S. pistillata* samples to use in plotting and ordinations

- `./sp_output/post_med_seqs/131_20201203_DBV_20201207T095144.seqs.absolute.abund_and_meta.txt`: post-MED ITS2 sequence absolute count table

- `./sp_output/its2_type_profiles/131_20201203_DBV_20201207T095144.profiles.absolute.abund_and_meta.txt`: ITS2 type profile absolute count table

- `./sp_output/between_sample_distances/A/20201207T095144_braycurtis_sample_distances_A_sqrt.dist`: between sample distances based on square root transformed ITS2 sequence abundances (post-MED).

- [only required if running with dist_type="uf" (UniFrac)] `./sp_output/between_sample_distances/A/20201207T095144_unifrac_sample_distances_A_sqrt.dist`

- `./131_20201203_DBV_20201207T095144.profiles.absolute.abund_and_meta.clustered.tsv`: ITS2 type profile absolute count table where the profiles have been clustered by having 3 or more DIVs in common.

- `./sp_output/between_profile_distances/A/20201207T095144_braycurtis_profile_distances_A_sqrt.dist`: between ITS2 type profile distances for *Symbiodinium* profiles (based on BrayCrutis and square root transformed counts)
- `./sp_output/between_profile_distances/C/20201207T095144_braycurtis_profile_distances_C_sqrt.dist`: between ITS2 type profile distances for *Cladocopium* profiles (based on BrayCrutis and square root transformed counts)
- `./sp_output/between_profile_distances/D/20201207T095144_braycurtis_profile_distances_D_sqrt.dist`: between ITS2 type profile distances for *Durusdinium* profiles (based on BrayCrutis and square root transformed counts)

# Ordination plots

Ordination plots were generated using the R script `plot_buitrago.R`.

The following files are required as input for this script:

- `./pver.genclust.strata.K2.csv`: genetic clusters assigned to *P. verrucosa*

- `./spis.genclust.strata.K6.csv`: genetic clusters assigned to *S. pistillata*

- `./sp_output/post_med_seqs/131_20201203_DBV_20201207T095144.seqs.absolute.abund_and_meta.txt`: post-MED ITS2 sequence absolute count table

# PERMANOVA and Beta diversity testing

The PERMANOVA and Beta were conducted using the permanove_buitrago.r script.

The following files are required as input:

- `./pver.genclust.strata.K2.csv`: genetic clusters assigned to *P. verrucosa*

- `./spis.genclust.strata.K6.csv`: genetic clusters assigned to *S. pistillata*

- `./sp_output/between_sample_distances/A/20201207T095144_braycurtis_sample_distances_A_sqrt.dist`: between sample distances based on square root transformed ITS2 sequence abundances (post-MED).

- `./reef_temp.csv`: file containing the sst_clim_mmm in deg C.

The PERMANOVAs were also computed using UniFrac-derived between sample distances. The computation of these PERMANOVAs has been commented out in the R script. To run these, the following file must be present:

- `./sp_output/between_sample_distances/A/20201207T095144_unifrac_sample_distances_A_sqrt.dist`
