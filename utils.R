getRefSeq_genes <- function(reference_genome) {
    makeGRanges(
        getUCSCtable("refGene", "RefSeq Genes", freeze=reference_genome),
        freeze=reference_genome
    )
}

get_integration_sites_with_mrcs <- function(sampleName) {
    sites <- getUniqueSites(sampleName)
    sites$type <- "insertion"

    mrcs <- getMRCs(sampleName)
    mrcs$type <- "match"

    sites_mrcs <- rbind(sites, mrcs)

    sites_mrcs <- makeGRanges(sites_mrcs, soloStart=TRUE,
        chromCol='chr', strandCol='strand', startCol='position') 
    sites_mrcs
}
