# Genomic Heatmap Maker for Integration Sites


# Script to run genomic heatmap maker from sites from DB

For list of GTSP and replicates generate heatmap
for a given reference genome and given integration site database:
```
Rscript genomic_heatmap_from_db.R sampleName_GTSP.csv -o heatmap --ref_genome hg18  --sites_group intSitesDev237
```
Group should be present in ~/.my.cnf.

File `sampleName_GTSP.csv` should have at least 2 columns: sampleName and GTSP.
If column `label` is given it will be used as a label in heatmap otherrwise GTSP will be used.
All replicates of given GTSP merged together. See geneTherapyPatientReportMaker and intSiteCaller
for instruction how to generate `sampleName_GTSP.csv`.

Defaults: reference genome is "hg18" and sites group is "intsites_miseq".


# Database configuration file 

Configuration file should be in home directory and called .my.cnf,
(~/.my.cnf).

The .my.cnf format is:

```
[group_name]
user=YYYYYYY
password=XXXXXX
host=microbZZZ.med.upenn.edu
port=330X
database=XXX
```

# Dependencies

intSiteRetriever
hiAnnotator
pipeUtils
colorspace
GCcontent

List of cancer genes are copied from:
CancerGeneList

# implementation details: pipeUtils requirement

needs a column type which is hard-coded as "insertion" for integration site
and "match" for match random contol(mrc).

pipeUtils can only generate figures for 2 or more samples.


