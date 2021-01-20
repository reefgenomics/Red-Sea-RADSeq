# PATHS
library(vegan)
library(glue)


PerformPermanova <- function(
  meta_location_path, meta_genotype_path, dist_path, title
){
  # A function for formatting the meta info and distance matrices
  # and ouputting the results of a PERMANOVA
  # SYMBIODINIUM dist matrix processing
  dist_df <- read.table(dist_path, sep="\t", header = FALSE)
  dist_names = dist_df$V1
  # drop the uids and the hostnames
  dist_df <- dist_df[,-c(1,2)]
  # set the row names to the host names
  rownames(dist_df) <- dist_names
  colnames(dist_df) <- dist_names
  
  # location meta info
  location_meta_df <- read.table(meta_location_path, sep='\t', header=TRUE)
  location_names<- location_meta_df$INDIVIDUALS
  rownames(location_meta_df) <- location_names
  # Get the intersect of the two set of names
  spis_names <- intersect(dist_names, location_names)
  # Drop the rows and columns not in the intersect
  dist_df <- dist_df[spis_names, spis_names]
  location_meta_df <- location_meta_df[spis_names,]
  
  # Now we have the matrix we want to work with for the PERMANOVA
  dist_matrix <- as.dist(data.matrix(dist_df))
  
  # Now we need to prepare the meta information
  # We want levels of REEF REGION and HOST_GENOTYPE
  # REEF and REGION levels are already contained in the spis_host_df
  # Get HOST_GENOTYPE
  genotype_meta_df <- read.table(meta_genotype_path, head=TRUE, sep='\t')
  rownames(genotype_meta_df) <- genotype_meta_df$INDIVIDUALS
  genotype_meta_df <- genotype_meta_df[spis_names,]
  meta_df <- merge(genotype_meta_df, location_meta_df,by="INDIVIDUALS")
  rownames(meta_df) <- meta_df$INDIVIDUALS
  meta_df <- meta_df[spis_names,]
  
  # PERMANOVA for spis_symbiodinium
  print(glue("PERMANOVA for {title}"))
  print("dist_matrix ~ REEF*GEN_CLUSTER")
  a_result <- adonis(formula = dist_matrix ~ REEF*GEN_CLUSTER, data=meta_df)
  print(a_result)
  cat("\n\n\n")
  print("dist_matrix ~ REGION*REEF*GEN_CLUSTER")
  a_result <- adonis(formula = dist_matrix ~ REGION*REEF*GEN_CLUSTER, data=meta_df)
  print(a_result)
  cat("\n\n\n")
  print("dist_matrix ~ GEN_CLUSTER*REGION*REEF")
  a_result <- adonis(formula = dist_matrix ~ GEN_CLUSTER*REGION*REEF, data=meta_df)
  print(a_result)
  cat("\n\n\n")
}


# PERMANOVAS
print("Performing PVER-BrayCurtis PERMANOVAs")
PerformPermanova(
  meta_location_path = "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/pver.14reef.277ind.ordered.strata.txt",
  meta_genotype_path = "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/Pver.major.host.genetic.cluster.txt",
  dist_path = "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/sp_output/between_sample_distances/A/20201207T095144_braycurtis_sample_distances_A_sqrt.dist",
  title="SPIS - BrayCurtis"
)

print("Performing PVER-UniFrac PERMANOVAs")
PerformPermanova(
  meta_location_path = "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/pver.14reef.277ind.ordered.strata.txt",
  meta_genotype_path = "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/Pver.major.host.genetic.cluster.txt",
  dist_path = "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/sp_output/between_sample_distances/A/20201207T095144_unifrac_sample_distances_A_sqrt.dist",
  title="SPIS - BrayCurtis"
)
print("Performing SPIS-BrayCurtis PERMANOVAs")
PerformPermanova(
  meta_location_path = "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/spis.18reef.368ind.ordered.strata.txt",
  meta_genotype_path = "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/Spis.major.host.genetic.cluster.txt",
  dist_path = "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/sp_output/between_sample_distances/A/20201207T095144_braycurtis_sample_distances_A_sqrt.dist",
  title="SPIS - BrayCurtis"
)

print("Performing SPIS-UniFrac PERMANOVAs")
PerformPermanova(
  meta_location_path = "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/spis.18reef.368ind.ordered.strata.txt",
  meta_genotype_path = "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/Spis.major.host.genetic.cluster.txt",
  dist_path = "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/sp_output/between_sample_distances/A/20201207T095144_unifrac_sample_distances_A_sqrt.dist",
  title="SPIS - BrayCurtis"
)




# # The path to the BrayCurtis distance output
# bc_distance_path <- file.path(".", "sp_output", "between_sample_distances", "A", "20201207T095144_unifrac_sample_distances_A_sqrt.dist")
# spis_host_sample_names_path = "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/spis.18reef.368ind.ordered.strata.txt"
# spis_host_genotype_path = "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/Spis.major.host.genetic.cluster.txt"
# pver_host_genotype_path = "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/Pver.major.host.genetic.cluster.txt"
# # First get the lists of samples that we will be working with
# # We will need to read in Carol's list of sample names that she is using
# # in the host analysis and then compare these to the samples that
# # we have in the distance SymPortal output
# # There will be samples in each of these that aren't in the other
# # Because: Carol did not analyse all samples submitted to SymPortal,
# # and the samples that contained very little Symbiodinium will not be in the
# # Symbiodinum distances matrix. As such we want to work with an intersaction
# # of the two lists
# 
# # SYMBIODINIUM dist matrix processing
# # Spis bc dist matrix
# symbiodinium_dist_df <- read.table(bc_distance_path, sep="\t", header = FALSE)
# symbiodinium_dist_names = symbiodinium_dist_df$V1
# # drop the uids and the hostnames
# symbiodinium_dist_df <- symbiodinium_dist_df[,-c(1,2)]
# # set the row names to the host names
# rownames(symbiodinium_dist_df) <- symbiodinium_dist_names
# colnames(symbiodinium_dist_df) <- symbiodinium_dist_names
# 
# 
# # Spis host samples
# spis_host_df <- read.table(spis_host_sample_names_path, sep='\t', header=TRUE)
# spis_host_names <- spis_host_df$INDIVIDUALS
# rownames(spis_host_df) <- spis_host_names 
# # Get the intersect of the two set of names
# spis_names <- intersect(symbiodinium_dist_names, spis_host_names)
# # Drop the rows and columns not in the intersect
# spis_dist_df <- symbiodinium_dist_df[spis_names, spis_names]
# spis_host_df <- spis_host_df[spis_names,]
# 
# # Now we have the matrix we want to work with for the PERMANOVA
# spis_matrix <- as.dist(data.matrix(spis_dist_df))
# 
# # Now we need to prepare the meta information
# # We want levels of REEF REGION and HOST_GENOTYPE
# # REEF and REGION levels are already contained in the spis_host_df
# # Get HOST_GENOTYPE
# spis_host_genotype_df <- read.table(spis_host_genotype_path, head=TRUE, sep='\t')
# rownames(spis_host_genotype_df) <- spis_host_genotype_df$INDIVIDUALS
# spis_host_genotype_df <- spis_host_genotype_df[spis_names,]
# spis_meta_df <- merge(spis_host_genotype_df, spis_host_df,by="INDIVIDUALS")
# rownames(spis_meta_df) <- spis_meta_df$INDIVIDUALS
# spis_meta_df <- spis_meta_df[spis_names,]
# 
# # PERMANOVA for spis_symbiodinium
# adonis2(formula = spis_matrix ~ REGION*REEF*GEN_CLUSTER, data=spis_meta_df)











