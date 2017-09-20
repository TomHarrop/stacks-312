#!/usr/bin/env Rscript

library(ape)
library(data.table)
library(ggplot2)
library(ggtree)
library(SNPRelate)


###########
# GLOBALS #
###########

set.seed(42)

# input VCF
vcf <- "output/stacks/batch_1.vcf"

# snprelate data
gds_file <- "output/snprelate/batch_1.gds"

# order of populations
population_order <- c("Ruakura", "Lincoln", "Fortrose")

##################
# PREPARE SNPSET #
##################

# convert to gds
snpgdsVCF2GDS(vcf.fn = vcf,
              out.fn = gds_file,
              method = "biallelic.only",
              verbose = TRUE)

# read in the GDS
gds <- snpgdsOpen(gds_file)

# prune the snps for LD
snpset <- snpgdsLDpruning(gds,
                          ld.threshold=0.2,
                          autosome.only = FALSE)
snpset_ids <- unlist(snpset)

#################################
# PRINCIPAL COMPONENTS ANALYSIS #
#################################

# run the pca
pca <- snpgdsPCA(gds,
                 snp.id = snpset_ids,
                 autosome.only = FALSE)

# set up plotting data
plotdata <- data.table(sample_id = pca$sample.id,
                       PCA1 = pca$eigenvect[, 1],
                       PCA2 = pca$eigenvect[, 2])
plotdata[, population := factor(gsub("[[:digit:]]", "", sample_id),
                                levels = population_order)]

# set up axis labels
pv_labs <- paste0(
    c("PCA1", "PCA2"),
    " (",
    signif(pca$varprop[c(1, 2)] * 100, 3),
    "%)")

# draw the plot
ggplot(plotdata, aes(x = PCA1, y = PCA2, colour = population)) +
    xlab(pv_labs[1]) + ylab(pv_labs[2]) +        
    scale_color_brewer(palette = "Set1",
                       guide = guide_legend(title = NULL)) +
    geom_point(size = 4, alpha = 0.8)

####################
# IDENTITY HEATMAP #
####################

# calculate identity by state for SNPs
ibs <- snpgdsIBS(gds,
                 snp.id = snpset_ids,
                 autosome.only = FALSE)

# hclust the IBS matrix (essentially an NJ-tree)
hc <- snpgdsHCluster(ibs, hang = 0.1)

# prepare plotting data.table
ibs_dt <- data.table(ibs$ibs)
setnames(ibs_dt, names(ibs_dt), ibs$sample.id)
ibs_dt[, id1 := ibs$sample.id]
ibs_pd <- melt(ibs_dt,
               id.vars = "id1",
               variable.name = "id2",
               value.name = "ibs")

# set the x- and y-axis order using the hclust results
sample_order <- hc$hclust$labels[hc$hclust$order]
ibs_pd[, id1 := factor(id1, levels = sample_order)]
ibs_pd[, id2 := factor(id2, levels = sample_order)]
ibs_pd[, pop1 := factor(gsub("[[:digit:]]", "", id1),
                        levels = population_order)]

# set up axis text colour list
Set1 <- RColorBrewer::brewer.pal(9, "Set1")
colour_list <- Set1[c(1:3)]
names(colour_list) <- population_order
setkey(ibs_pd, id1)
col_order <- unique(ibs_pd[, colour_list[pop1], by = id1])$V1

# draw the plot
heatscale <- RColorBrewer::brewer.pal(6, "YlOrRd")
ggplot(ibs_pd, aes(x = id1, y = id2, fill = ibs)) +
    theme_minimal() +
    theme(axis.text.y = element_text(colour = col_order),
          axis.text.x = element_text(colour = col_order,
                                     angle = 90,
                                     hjust = 1)) +
    xlab(NULL) + ylab(NULL) +
    scale_fill_gradientn(colours = heatscale,
                         guide = guide_colourbar(title = "Identity")) +
    geom_raster()

###############################
# CLADOGRAMS, NOT THAT USEFUL #
###############################

# with SNPRelate
ct <- snpgdsCutTree(hc)
snpgdsDrawTree(ct)

# draw cladogram with ggtree
hc_phylo <- nj(as.dist(hc$dist))
gt <- ggtree(hc_phylo, layout = "fan", branch.length = "none")

# add labels
gt %<+% data.frame(plotdata) +
    geom_tiplab2(aes(colour = population)) +
    scale_color_brewer(palette = "Set1",
                       guide = guide_legend(title = NULL))

###################
# SNPRelate FSTs? #
###################

# wrap snpgdsFst
RunFstCalc <- function(s1, s2, gds, snpset_ids){
    my_samples <- unique(c(s1, s2))
    my_populations <- factor(gsub("[[:digit:]]", "", my_samples))
    snpgdsFst(gds,
              snp.id = snpset_ids,
              autosome.only = FALSE,
              sample.id = my_samples,
              population = my_populations,
              method = "W&C84")
}

# set up samples for pairwise comparison
all_samples <- read.gdsn(index.gdsn(gds, "sample.id"))
rua_samples <- grep("Ruakura", all_samples, value = TRUE)
lin_samples <- grep("Lincoln", all_samples, value = TRUE)
for_samples <- grep("Fortrose", all_samples, value = TRUE)

# check all
fst_all <- RunFstCalc(all_samples,
                      all_samples,
                      gds = gds,
                      snpset_ids = snpset_ids)

# run pairwise
pw_fst <- list(
    rua_lin = RunFstCalc(rua_samples, lin_samples, gds, snpset_ids),
    rua_for = RunFstCalc(rua_samples, for_samples, gds, snpset_ids),
    for_lin = RunFstCalc(for_samples, lin_samples, gds, snpset_ids)
)

# make a table
sapply(pw_fst, function(x) x$Fst)



