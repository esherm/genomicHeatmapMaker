source("genomicHeatmapMaker.R")

library(argparse)
library(DBI)
library(RMySQL)

parser <- ArgumentParser(description="Make genomic heatmap for sites from database")
parser$add_argument("-r", "--reference_genome", type="character", required=TRUE, 
    help="reference genome used for all samples")
parser$add_argument("-o", "--output_dir", type="character", required=TRUE, 
    help="output folder where genomic heat maps files will be saved")
parser$add_argument("-g", "--gtsp_label", type="character", required=TRUE, 
    help="tab separated file: GTSP# LABEL_FOR_HEAT_MAP")
parser$add_argument("-c", "--connection_group", type="character", required=TRUE, 
    help="which group to use for connection")

args <- parser$parse_args()

connection <- dbConnect(MySQL(), group=args$connection_group)

referenceGenome <- args$reference_genome
heat_map_result_dir <- args$output_dir 
gtsp_label <- read.delim(args$gtsp_label, stringsAsFactors=FALSE)

sampleNameInput <- gtsp_label[, 1]
# for compatibility with intSiteRetriever
sampleNameInput <- sapply(sampleNameInput, function(x) {paste0(x, '%')})
labels <- gtsp_label[, 2]

make_heatmap(sampleNameInput, labels, 
    referenceGenome, heat_map_result_dir, connection)
