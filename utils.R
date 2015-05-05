getRefSeq_genes <- function(reference_genome) {
    refSeq <- makeGRanges(
        getUCSCtable("refGene", "RefSeq Genes", freeze=reference_genome),
        freeze=reference_genome
    )
}

getCpG_islands <- function(reference_genome) {
    cpg <- getUCSCtable("cpgIslandExt", "CpG Islands", freeze=reference_genome)
    cpg$strand <- "*" # either strand
    makeGRanges(cpg, freeze=reference_genome, chromCol='chrom')
}

getDNaseI <- function(reference_genome) {
    DNaseI <- getUCSCtable("wgEncodeRegDnaseClustered", 
        "DNase Clusters", freeze=reference_genome)
    DNaseI$strand <- "*" # either strand
    makeGRanges(DNaseI, freeze=reference_genome, chromCol='chrom')
}

get_integration_sites_with_mrcs <- function(sampleName, refGenomeSeq) {
    sites <- getUniqueSites(sampleName)
    sites$type <- "insertion"

    mrcs <- getMRCs(sampleName)
    mrcs$type <- "match"

    sites_mrcs <- rbind(sites, mrcs)

    sites_mrcs <- makeGRanges(sites_mrcs, soloStart=TRUE,
        chromCol='chr', strandCol='strand', startCol='position')

    #seqinfo needs to be exact here or trimming will be wrong
    newSeqInfo <- seqinfo(refGenomeSeq)
    seqInfo.new2old <- match(seqnames(newSeqInfo),
                             seqnames(seqinfo(sites_mrcs)))
    seqinfo(sites_mrcs, new2old=seqInfo.new2old) <- newSeqInfo

    sites_mrcs
}

#' return genome seq for human readable UCSC format
#'
#' format is: hg18, ...
get_reference_genome <- function(reference_genome) {
  pattern <- paste0("\\.", reference_genome, "$")
  match_index <- which(grepl(pattern, installed.genomes()))
  stopifnot(length(match_index) == 1)
  BS_genome_full_name <- installed.genomes()[match_index]
  get(BS_genome_full_name)
}

get_annotation_columns <- function(sites) {
    granges_column_names <- c("seqnames", "start", "end", "width", "strand")
    int_site_column_names <- c("siteID", "sampleName", "chr", "strand", "position")
    required_columns <- unique(c(
        granges_column_names, int_site_column_names, "type"))
    stopifnot(all(required_columns %in% names(sites)))
    setdiff(names(sites), required_columns)
}

from_counts_to_density <- function(sites, column_prefix, window_size) {
    metadata <- mcols(sites)
    sapply(seq(window_size), function(i) {
        val <- window_size[i]
        name <- names(window_size)[i]
        column_name <- paste0(column_prefix, ".", name)
        metadata[[column_name]] <<- metadata[[column_name]]/val
    })
    mcols(sites) <- metadata
    sites
}

getPositionalValuesOfFeature <- function(sites, genomicData) {
    #### Boundary Distances #### Nirav Malani code TODO: refactor into several functions
    ## (refSeq boundary.dist), Start (refSeq start.dist), non-width (), General (general.width)
    ## when inGene is FALSE then set following: ref.left.pos, ref.right.pos, ref.left.strand, ref.right.strand
    ## when inGene is TRUE then set following: ref.start.pos, ref.end.pos, ref.gene.strand

    ## prepare the new columns ##
    colnam <- paste("ref", c("left.pos", "right.pos", "left.strand", "right.strand", 
                             "start.pos", "end.pos", "gene.strand"), sep=".") 
    mcols(sites)[colnam] <- NA

    ## add the respective columns as needed ##
    ## beware: precede returns range which is following the query and
    ## follow returns the range which is preceding the query!
    ## so do a switcheroo in terms of extracting the start & stop ##
    left <- follow(sites, genomicData, ignore.strand=TRUE)
    left[is.na(left) | sites$within_refSeq_gene] <- NA
    rows <- na.omit(left)
    sites$ref.left.pos[!is.na(left)] <- end(genomicData[rows])
    sites$ref.left.strand[!is.na(left)] <- as.character(strand(genomicData[rows]))

    right <- precede(sites, genomicData, ignore.strand=TRUE)
    right[is.na(right) | sites$within_refSeq_gene] <- NA
    rows <- na.omit(right)
    sites$ref.right.pos[!is.na(right)] <- start(genomicData[rows])
    sites$ref.right.strand[!is.na(right)] <- as.character(strand(genomicData[rows]))

    inIt <- findOverlaps(sites, genomicData, ignore.strand=TRUE, select="arbitrary")
    inIt[is.na(inIt) | !sites$within_refSeq_gene] <- NA
    rows <- na.omit(inIt)
    sites$ref.start.pos[!is.na(inIt)] <- start(genomicData[rows])
    sites$ref.end.pos[!is.na(inIt)] <- end(genomicData[rows])
    sites$ref.gene.strand[!is.na(inIt)] <- as.character(strand(genomicData[rows]))

    sites$boundary.dist <-
        eval(expression(pmin((ref.end.pos-position)/(ref.end.pos-ref.start.pos),
                             (position-ref.start.pos)/(ref.end.pos-ref.start.pos),
                             (ref.right.pos-position)/(ref.right.pos-ref.left.pos),
                             (position-ref.left.pos)/(ref.right.pos-ref.left.pos),
                             na.rm=T)), mcols(sites))

    sites$start.dist <-
        eval(expression(pmin(ifelse(ref.gene.strand=="-",
                                    (ref.end.pos-position)/(ref.end.pos-ref.start.pos),
                                    (position-ref.start.pos)/(ref.end.pos-ref.start.pos)),
                             ifelse(ref.right.strand=="-",
                                    (ref.right.pos-position)/(ref.right.pos-ref.left.pos),
                                    NA),
                             ifelse(ref.left.strand=="+",
                                    (position-ref.left.pos)/(ref.right.pos-ref.left.pos),
                                    NA),na.rm=T)), mcols(sites))

    sites$general.width <- eval(expression(pmin(ref.end.pos-ref.start.pos, 
                                                ref.right.pos-ref.left.pos,na.rm=T)),
                                mcols(sites))
    sites$gene.width <- eval(expression(ref.end.pos-ref.start.pos ), mcols(sites))

    meta <- mcols(sites)
    meta <- meta[ , ! (names(meta) %in% colnam)]
    mcols(sites) <- meta

    sites 
}
