source("genomicHeatmapMaker.R")
source("utils.R")

library(argparse, quietly=TRUE)
library(DBI, quietly=TRUE)
library(RMySQL, quietly=TRUE)
library(dplyr)

parser <- ArgumentParser(description="Make genomic heatmap for sites from database")
parser$add_argument("sample_gtsp", nargs='?', default='sampleName_GTSP.csv')
parser$add_argument("-o", "--output_dir", type="character", default="heatmap_output",
    help="output folder where genomic heat maps files will be saved")
parser$add_argument("-r", "--ref_genome", type="character", default="hg18", 
    help="reference genome used for all samples")
parser$add_argument("-s", "--sites_group", type="character", default="intsites_miseq.read", 
    help="which group to use for connection")

args <- parser$parse_args()
args

referenceGenome <- args$ref_genome
heat_map_result_dir <- args$output_dir 

csvfile <- args$sample_gtsp
if( ! file.exists(csvfile) ) stop(csvfile, " not found")
sampleName_GTSP <- read.csv(csvfile)
stopifnot(all(c("sampleName", "GTSP") %in% colnames(sampleName_GTSP)))
message("\nGenerating report from the following sets")
print(sampleName_GTSP)

connection <- dbConnect(MySQL(), group=args$sites_group)
info <- dbGetInfo(connection)
connection <- src_sql("mysql", connection, info = info)

make_heatmap(sampleName_GTSP, referenceGenome, heat_map_result_dir, connection)
