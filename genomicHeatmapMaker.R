source("intSiteRetriever/intSiteRetriever.R")
source("utils.R")

library(colorspace)
library(hiAnnotator)
library(pipeUtils)

referenceGenome <- "hg18"
heat_map_result_dir <- "./heatmap"

# should have at least two samples
sampleName <- c("pool1-1", "HIV_CTRL_noLig-1")
stopifnot(length(sampleName) != 1)
# check that all samples processed with the same reference genome
stopifnot(unique(getRefGenome(sampleName)$refGenome) == referenceGenome)
stopifnot(all(setNameExists(sampleName)))

sites_mrcs <- get_integration_sites_with_mrcs(sampleName)

# TODO: populate from local database, at present pulled from UCSC web-site
refSeq_genes <- getRefSeq_genes(referenceGenome)
# END annotation loading

sites_mrcs <- getSitesInFeature(
            sites_mrcs, refSeq_genes, "within_refSeq_gene", asBool=TRUE)

window_size_refSeq <- c("10k"=1e4, "100k"=1e5, "1M"=1e6)
sites_mrcs <- getFeatureCounts(sites_mrcs, refSeq_genes, "refSeq.density", 
                          width=window_size_refSeq)

granges_column_names <- c("seqnames", "start", "end", "width", "strand")
int_site_column_names <- c("siteID", "sampleName", "chr", "strand", "position")
required_columns <- unique(c(
    granges_column_names, int_site_column_names, "type"))

sites_mrcs <- as.data.frame(sites_mrcs)
stopifnot(all(required_columns %in% names(sites_mrcs)))
annotation_columns <- setdiff(names(sites_mrcs), required_columns)

rset <- with(sites_mrcs, ROC.setup(
    rep(TRUE, nrow(sites_mrcs)), type, siteID, sampleName))
roc.res <- ROC.strata(annotation_columns, rset, add.var=TRUE, sites_mrcs)
ROCSVG(roc.res, heat_map_result_dir)


