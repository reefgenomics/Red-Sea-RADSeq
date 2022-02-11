# PATHS
library(vegan)
library(glue)
library(stringr)
library(devtools)
# Install instructions here: https://github.com/pmartinezarbizu/pairwiseAdonis
library(pairwiseAdonis)

PerformPermanova <- function(
  genotype.meta.pth, sym.dist.path, title
){
  # A function for formatting the reef and host genotype 
  # meta info and Symbiodiniaceae distance matrices
  # and ouputting the results of a PERMANOVA
  # SYMBIODINIUM dist matrix processing
  
  # The samples we are working with are defined by the samples in the
  # genotype.meta.pth file. As such, process this first,
  # then find the intersect with the Symbiodiniaceae distance
  # and use this intersect to contduct the analyses.
  
  # Get HOST_GENOTYPE
  genotype.meta.df <- read.table(genotype.meta.pth, head=TRUE, sep=',')
  sample_names_to_work_with <- genotype.meta.df$INDIVIDUALS
  rownames(genotype.meta.df) <- sample_names_to_work_with
  
  # Need to check that all of these are in the Symbiodinium distance matrix
  sym.dist.df <- read.table(sym.dist.path, sep="\t", header = FALSE)
  dist_names = sym.dist.df$V1
  # drop the uids and the hostnames
  sym.dist.df <- sym.dist.df[,-c(1,2)]
  # set the row names to the host names
  rownames(sym.dist.df) <- dist_names
  colnames(sym.dist.df) <- dist_names
  
  # Check to see if all of the sample_names_to_work_with names
  # are found in the sym.dist.df
  all(sample_names_to_work_with %in% rownames(sym.dist.df))
  table(sample_names_to_work_with %in% rownames(sym.dist.df))
  # FALSE  TRUE 
  # 1   202 
  # NB one of the samples in the host_genotype file is missing from the
  # ITS2 sample_names
  
  # Get a list of the names in common and work with this
  sample_names_to_work_with <- intersect(sample_names_to_work_with, rownames(sym.dist.df))
  
  # Then subset the genotype and zooxs matrix to the above intersect list
  genotype.meta.df <- subset(genotype.meta.df, rownames(genotype.meta.df) %in% sample_names_to_work_with)
  sym.dist.df <- sym.dist.df[sample_names_to_work_with, sample_names_to_work_with]

  # Generate reef info by using the name
  genotype.meta.df <- transform(genotype.meta.df, REEF=str_extract(INDIVIDUALS, "[A-Z]{3}-R\\d+"))
  # TODO here we should join the temperature data to the genotype.meta.df
  reef.temp.df = read.csv("reef_temp.csv")
  rownames(reef.temp.df) = reef.temp.df$reef
  colnames(reef.temp.df)[1] = "REEF"
  genotype.meta.df = transform(genotype.meta.df, temp=reef.temp.df[REEF, "temp"])
  genotype.meta.df$tempf<-
    with(genotype.meta.df,
         ifelse(temp >= 27 & temp < 28, "cool",
                ifelse(temp >= 29 & temp < 30, "med-cool",
                       ifelse(temp >= 30 & temp < 31, "med-warm",
                              ifelse(temp >= 31 & temp < 32, "warm", "other")))))
  genotype.meta.df$tempf <- as.factor(genotype.meta.df$tempf)
    
  # Change ambiguous STRATA name to GEN_CLUSTER
  colnames(genotype.meta.df)[2] <- "GEN_CLUSTER"

  # Now we have the matrix we want to work with for the PERMANOVA
  sym.dist.df <- as.dist(data.matrix(sym.dist.df))
  
  # PERMANOVA
  sink(glue("{title}.permanova.results.main.txt"))
  print(glue("PERMANOVA for {title}"))
  print("dist_matrix ~ temp*GEN_CLUSTER")
  a_result <- adonis(formula = sym.dist.df ~ tempf*GEN_CLUSTER, data=genotype.meta.df)
  print(a_result)
  cat("\n\n\n")
  
  print("dist_matrix ~ GEN_CLUSTER*temp")
  a_result <- adonis(formula = sym.dist.df ~ GEN_CLUSTER*tempf, data=genotype.meta.df)
  print(a_result)
  cat("\n\n\n")
  sink()
  
  sink(glue("{title}.permanova.results.pairwise.txt"))
  print("Pairwsie comparisons")
  print(pairwise.adonis(sym.dist.df, genotype.meta.df$GEN_CLUSTER, p.adjust.m="fdr"))
  print(pairwise.adonis(sym.dist.df, genotype.meta.df$tempf, p.adjust.m="fdr"))
  sink()
  
  sink(glue("{title}.permanova.results.strata.txt"))
  print(adonis(formula = sym.dist.df ~ GEN_CLUSTER, strata=genotype.meta.df$tempf, data=genotype.meta.df))
  sink()
  
  # BETADISPER
  # REEF
  mod.tempf = betadisper(sym.dist.df, genotype.meta.df$tempf)
  sink(glue("{title}.reef.betadisper.anova.txt"))
  print(anova(mod.tempf))
  print("\n\n\n")
  # mod_reef_type.HSD = TukeyHSD(mod_reef_type)
  print(permutest(mod.tempf, pairwise=TRUE))
  sink()
  
  pdf(file=glue("{title}.reef.betadisper.plot.pdf"), height=10, width=10)
  plot(mod.tempf, segments = FALSE, ellipse=TRUE, label = TRUE)
  dev.off()
  
  # GEN_CLUSTER
  mod.gen.cluster = betadisper(sym.dist.df, genotype.meta.df$GEN_CLUSTER)
  sink(glue("{title}.gen_cluster.betadisper.anova.txt"))
  print(anova(mod.gen.cluster))
  # mod_reef_type.HSD = TukeyHSD(mod_reef_type)
  print(permutest(mod.gen.cluster, pairwise=TRUE))
  sink()
  
  pdf(file=glue("{title}.gen_cluster.betadisper.plot.pdf"), height=10, width=10)
  plot(mod.gen.cluster, segments = FALSE, ellipse=TRUE, label=TRUE)
  dev.off()
}


# PERMANOVAS
# NB the braycurtis gave smaller residuals in both species so we will work
# with those.
print("Performing PVER-BrayCurtis PERMANOVAs")
PerformPermanova(
  genotype.meta.pth = "pver.genclust.strata.K2.csv",
  sym.dist.path = "sp_output/between_sample_distances/A/20201207T095144_braycurtis_sample_distances_A_sqrt.dist",
  title="pver.bc"
)

# print("Performing PVER-UniFrac PERMANOVAs")
# PerformPermanova(
#   genotype.meta.pth = "pver.genclust.strata.K2.csv",
#   sym.dist.path = "sp_output/between_sample_distances/A/20201207T095144_unifrac_sample_distances_A_sqrt.dist",
#   title="PVER - UniFrac"
# )

print("Performing SPIS-BrayCurtis PERMANOVAs")
PerformPermanova(
  genotype.meta.pth = "spis.genclust.strata.K6.csv",
  sym.dist.path = "sp_output/between_sample_distances/A/20201207T095144_braycurtis_sample_distances_A_sqrt.dist",
  title="spis.bc"
)

# print("Performing SPIS-UniFrac PERMANOVAs")
# PerformPermanova(
#   genotype.meta.pth = "spis.genclust.strata.K6.csv",
#   sym.dist.path = "sp_output/between_sample_distances/A/20201207T095144_unifrac_sample_distances_A_sqrt.dist",
#   title="SPIS - UniFrac"
# )
