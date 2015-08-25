source("utils.R")

library(colorspace)
library(hiAnnotator)
library(pipeUtils)
library(intSiteRetriever)
library(GCcontent)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg18)
library(BSgenome.Hsapiens.UCSC.hg19)

make_heatmap <- function(sampleName_GTSP, referenceGenome, output_dir, connection) {
    if ( ! "label" %in% colnames(sampleName_GTSP)) {
        sampleName_GTSP$label <- sampleName_GTSP$GTSP
    }
    sampleName_GTSP <- select(sampleName_GTSP, sampleName, GTSP, label)

    # should have at least two samples
    stopifnot(length(unique(sampleName_GTSP$GTSP)) != 1)

    sampleName_GTSP$refGenome <- rep(referenceGenome, nrow(sampleName_GTSP))

    # check that all samples processed with the same reference genome
    stopifnot(all(setNameExists(sampleName_GTSP, connection)))

    reference_genome_sequence <- get_reference_genome(referenceGenome)
    sites_mrcs <- get_integration_sites_with_mrcs(
        sampleName_GTSP, reference_genome_sequence, connection)

    # TODO: populate from local database, at present pulled from UCSC web-site
    refSeq_genes <- getRefSeq_genes(referenceGenome)
    CpG_islands <- getCpG_islands(referenceGenome)
    DNaseI <- getDNaseI(referenceGenome)

    oncogene_file <- "allonco_no_pipes.csv"
    # @return vector of gene symbols
    get_oncogene_from_file <- function(filename) {
        onco <- read.csv(filename, header=FALSE, stringsAsFactors=FALSE)
        as.character(onco$V1)
    }
    oncogenes <- get_oncogene_from_file(oncogene_file)
    # END annotation loading

    sites_mrcs <- getSitesInFeature(
      sites_mrcs, refSeq_genes, "within_refSeq_gene", asBool=TRUE)

    # is there oncogene closer than 50k
    refSeq_gene_symbols <- refSeq_genes$name2
    #' check if gene is onco gene list(curated by Bushman's lab)
    #' @return TRUE if onco-gene FALSE if not
    is_onco_gene <- function(gene_symbol_sites, oncogenes) {
        toupper(gene_symbol_sites) %in% toupper(oncogenes)
    }
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
    # temporary smaller one to save CPU time:
    # for 50K sites and 150K mrcs takes about several minutes
    window_size_GC <- c("50"=50, "100"=100, "250"=250,
        "500"=500, "1k"=1000, "2k"=2000, "5k"=5000, 
        "10k"=1e4, "25k"=2.5e4, "50k"=5e4)
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

    #restore ordering of values in sites_mrcs$sampleName column so that heatmap
    #order reflects sample input order
    #sites_mrcs$sampleName <- factor(sites_mrcs$sampleName, levels=names(sampleNameInput))

    rset <- with(sites_mrcs, ROC.setup(
      rep(TRUE, nrow(sites_mrcs)), type, siteID, sampleName))
    roc.res <- ROC.strata(annotation_columns, rset, add.var=TRUE, sites_mrcs)
    ROCSVG(roc.res, heat_map_result_dir)

}
