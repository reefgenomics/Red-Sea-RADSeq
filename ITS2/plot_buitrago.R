library(stringr)
library(phyloseq)
library(patchwork)
library(gridExtra)

# Plot up ordinations
pver.gentype.path <- "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/pver.genclust.strata.K2.csv"
spis.gentype.path <- "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/spis.genclust.strata.K6.csv"
symportal.post.med.count.table.path <- "/Users/benjaminhume/Documents/projects/20210113_buitrago/ITS2/sp_output/post_med_seqs/131_20201203_DBV_20201207T095144.seqs.absolute.abund_and_meta.txt"

# PVER meta df
pver.meta.df <- read.table(pver.gentype.path, head=TRUE, sep=',')
pver.meta.df <- base::transform(pver.meta.df, REEF=str_extract(pver.meta.df$INDIVIDUALS, "[A-Z]{3}-R\\d+"))
colnames(pver.meta.df)[2] <- "GEN_CLUSTER"
pver.meta.df = base::transform(pver.meta.df, SPECIES="PVER")

# SPIS meta df
spis.meta.df <- read.table(spis.gentype.path, head=TRUE, sep=',')
spis.meta.df <- base::transform(spis.meta.df, REEF=str_extract(spis.meta.df$INDIVIDUALS, "[A-Z]{3}-R\\d+"))
colnames(spis.meta.df)[2] <- "GEN_CLUSTER"
# Drop SWAJ-R1-43
spis.meta.df = spis.meta.df[!(spis.meta.df$INDIVIDUALS %in% c("SWAJ-R1-43")),]
spis.meta.df = base::transform(spis.meta.df, SPECIES="SPIS")

all.meta.df = rbind(pver.meta.df, spis.meta.df)
rownames(all.meta.df) = all.meta.df$INDIVIDUALS

# We will compute a braycurtis on the post_MED table
post.med.df = read.table(symportal.post.med.count.table.path, sep="\t", header=TRUE, row.names = 2, check.names=FALSE)
post.med.df = head(post.med.df, -1)
post.med.df = post.med.df[,-c(1:33)]
# Drop the samples not to be plotted
post.med.df = post.med.df[(rownames(post.med.df) %in% all.meta.df$INDIVIDUALS),]
# Drop samples with no seq info
post.med.df = post.med.df[which(rowSums(post.med.df) != 0),]
# Drop 0 seqs
post.med.df = post.med.df[,which(colSums(post.med.df) != 0)]
# Make compositional and square root transform
post.med.df.sqrt = sqrt(post.med.df)

otu.t= otu_table(post.med.df.sqrt, taxa_are_rows=FALSE)
sam.t= sample_data(all.meta.df)

phyloseq.all = phyloseq(otu.t, sam.t)
# NB microbiome package masked the base::transform method used above so
# we will load and then detach
library(microbiome)
summarize_phyloseq(phyloseq.all)


phyloseq.pver=subset_samples(phyloseq.all, SPECIES=="PVER")
phyloseq.spis=subset_samples(phyloseq.all, SPECIES=="SPIS")

P6=c("#222E50", "#007991", "#BCD8C1", "#E9D985", "#F29469", "#BE3F23") #27, 29, 32, 34
P10=c("#669900","#99cc33","#ccee66","#006699","#3399cc","#990066","#cc3399","#ff6600","#ff9900","#ffcc00", "#ADADAD")
P10=c("#669900","#99cc33","#ccee66","#006699","#3399cc","#990066","#cc3399","#ff6600","#ff9900","#ffcc00",  "#2077b4","#aec7e8","#ff7f0e","#ffbb78","#2ca02c","#98df89","#17becf","#9edae5","#e377c2","#f7b6d2","#ADADAD")
P20=c("#fad390", "#f6b93b", "#fa983a", "#e58e26", "#f8c291", "#e55039", "#eb2f06", "#b71540", "#6a89cc", "#4a69bd","#1e3799", "#0c2461", "#82ccdd", "#60a3bc", "#3c6382", "#0a3d62", "#b8e994", "#78e08f", "#38ada9", "#079992", "#C0C0C0")

pver.ord = ordinate(phyloseq.pver, method = "PCoA", distance = "bray")
pver.plot.reef = plot_ordination(physeq=phyloseq.pver, ordination=pver.ord, color="REEF") + geom_point(size = 3, alpha = 1) + theme_bw() + ggtitle("PVER.REEF") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P20)
pver.plot.gen_cluster = plot_ordination(physeq=phyloseq.pver, ordination=pver.ord, color="GEN_CLUSTER") + geom_point(size = 3, alpha = 1) + theme_bw() + ggtitle("PVER.GEN_CLUSTER") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P6)

spis.ord = ordinate(phyloseq.spis, method = "PCoA", distance = "bray")
spis.plot.reef = plot_ordination(physeq=phyloseq.spis, ordination=spis.ord, color="REEF") + geom_point(size = 3, alpha = 1) + theme_bw() + ggtitle("SPIS.REEF") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P20)
spis.plot.gen_cluster = plot_ordination(physeq=phyloseq.spis, ordination=spis.ord, color="GEN_CLUSTER") + geom_point(size = 3, alpha = 1) + theme_bw() + ggtitle("SPIS.GEN_CLUSTER") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P6)

pdf("plots/ordinations.pdf", width=10,height=10, pointsize = 10)
gridExtra::grid.arrange( pver.plot.reef, spis.plot.reef, pver.plot.gen_cluster, spis.plot.gen_cluster, ncol=2, nrow=2)
dev.off()
detach("package:microbiome", unload=TRUE)
