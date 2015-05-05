source("intSiteRetriever/intSiteRetriever.R")
source("GCcontent/GCcontent.R")
source("CancerGeneList/onco_genes.R")
source("utils.R")

library(colorspace)
library(hiAnnotator)
library(pipeUtils)

referenceGenome <- "hg18"
heat_map_result_dir <- "./heatmap"
sampleName <- c("pool1-1", "E9")

# should have at least two samples
stopifnot(length(sampleName) != 1)
# check that all samples processed with the same reference genome
stopifnot(unique(getRefGenome(sampleName)$refGenome) == referenceGenome)
stopifnot(all(setNameExists(sampleName)))

reference_genome_sequence <- get_reference_genome(referenceGenome)
sites_mrcs <- get_integration_sites_with_mrcs(sampleName, reference_genome_sequence)

# TODO: populate from local database, at present pulled from UCSC web-site
refSeq_genes <- getRefSeq_genes(referenceGenome)
CpG_islands <- getCpG_islands(referenceGenome)
DNaseI <- getDNaseI(referenceGenome)


oncogene_file <- "CancerGeneList/allonco_no_pipes.csv"
oncogenes <- get_oncogene_from_file(oncogene_file)
# END annotation loading

sites_mrcs <- getSitesInFeature(
            sites_mrcs, refSeq_genes, "within_refSeq_gene", asBool=TRUE)

# is there oncogene closer than 50k
refSeq_gene_symbols <- refSeq_genes$name2
is_refSeq_oncogene <- is_onco_gene(refSeq_gene_symbols, oncogenes)
refSeq_oncogene <- refSeq_genes[is_refSeq_oncogene]
sites_mrcs <- getNearestFeature(
    sites_mrcs, refSeq_oncogene, dists.only=TRUE, colnam="onco")
    #sites_mrcs, refSeq_oncogene, dists.only=TRUE, colnam="onco.100k")
sites_mrcs$onco.100k <- abs(sites_mrcs$oncoDist) <= 50000
sites_mrcs$oncoDist <- NULL
# end oncogene

sites_mrcs <- getPositionalValuesOfFeature(sites_mrcs, refSeq_genes)

window_size_refSeq <- c("10k"=1e4, "100k"=1e5, "1M"=1e6)
sites_mrcs <- getFeatureCounts(sites_mrcs, refSeq_genes, "refSeq_counts", 
                          width=window_size_refSeq)

window_size_GC <- c("50"=50, "100"=100, "250"=250,
    "500"=500, "1k"=1000, "2k"=2000, "5k"=5000,
    "10k"=1e4, "25k"=2.5e4, "50k"=5e4, "100k"=1e5, 
    "250k"=2.5e5, "500k"=5e5, "1M"=1e6)
sites_mrcs <- getGCpercentage(
    sites_mrcs, "GC", window_size_GC, reference_genome_sequence)

window_size_CpG_counts <- c("2k"=2e3, "10k"=1e4)
sites_mrcs <- getFeatureCounts(sites_mrcs, CpG_islands, "CpG_counts", 
                          width=window_size_CpG_counts)

window_size_CpG_density <- c("10k"=1e4, "100k"=1e5, "1M"=1e6)
sites_mrcs <- getFeatureCounts(sites_mrcs, CpG_islands, "CpG_density", 
                          width=window_size_CpG_density)
sites_mrcs <- from_counts_to_density(sites_mrcs, 
    "CpG_density", window_size_CpG_density)

window_size_DNaseI <- c("1k"=1e3, "10k"=1e4, "100k"=1e5, "1M"=1e6)
sites_mrcs <- getFeatureCounts(sites_mrcs, DNaseI, "DNaseI_count", 
                          width=window_size_DNaseI)

sites_mrcs <- as.data.frame(sites_mrcs)

annotation_columns <- get_annotation_columns(sites_mrcs)

rset <- with(sites_mrcs, ROC.setup(
    rep(TRUE, nrow(sites_mrcs)), type, siteID, sampleName))
roc.res <- ROC.strata(annotation_columns, rset, add.var=TRUE, sites_mrcs)
ROCSVG(roc.res, heat_map_result_dir)


