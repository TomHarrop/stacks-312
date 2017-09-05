#!/usr/bin/env Rscript

log_file <- snakemake@log[[1]]
sink(log_file, type = "message")

library(SNPRelate)

set.seed(42)

# read the gds file
gds_file = snakemake@input[["gds"]]
gds <- snpgdsOpen(gds_file)

# generate snpset
snpset <- snpgdsLDpruning(gds,
                          ld.threshold=0.2,
                          missing.rate = 0.9,
                          autosome.only = FALSE,
                          num.thread = snakemake@threads)
snpset_ids <- unlist(snpset)

# process pca
pca <- snpgdsPCA(gds,
                 snp.id = snpset_ids,
                 num.thread = snakemake@threads)

# save output
saveRDS(pca,
        snakemake@output[["rds"]])

# close log
sink(type="message")
close(log_file)
