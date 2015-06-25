# Database config file location

config file should be in home directory and called .my.cnf,
e.g. ~/.my.cnf

The .my.cnf format is as follows:

```
[intSitesDEV-dev]
user=YYYYYYY
password=XXXXXX
host=microb98.med.upenn.edu
port=3309
database=intsitesdev
```

# Script to run maker from sites from DB

```
Rscript map_from_db.R -r hg18 -o heatmap -g gtsp_label.tsv -c intSitesDev237
```

# Dependencies

intSiteRetriever see https://github.com/BushmanLab/intSiteRetriever.git
for installation instructions.

hiAnnotator
pipeUtils
colorspace

GCcontent: https://github.com/anatolydryga/GCcontent
CancerGeneList: https://github.com/BushmanLab/CancerGeneList.git

last 2 dependenices should be cloned in genomicHeatmapMaker directory.

# pipeUtils requirement

needs a column type which is hard-coded as "insertion" for integration site
and "match" for match random contol(mrc).


