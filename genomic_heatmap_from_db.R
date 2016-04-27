source("genomicHeatmapMaker.R")
source("utils.R")

libs <- c("argparse", "DBI", "RMySQL", "dplyr")
invisible(sapply(libs, library, character.only=TRUE))

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

loaded_ref_genomes <- c ("hg18", "mm9")
if ( ! referenceGenome %in% loaded_ref_genomes) {
    message("Only following genomes are loaded:")
    message(paste(loaded_ref_genomes, collapse=" "))
    message("Install and add new genomes to genomicHeatmapMaker.R")
    message("and add it to loaded_ref_genomes vector in genomic_heatmap_from_db.R")
    stop(0)
}

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
