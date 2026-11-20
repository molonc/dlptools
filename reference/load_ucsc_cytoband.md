# load cytoband file

file loader for other functions. Targets USCS cytoband files. Cytoband
files originally obtained from UCSC golden path genomes and parsed to
collect the total centromere spans.

## Usage

``` r
load_ucsc_cytoband(version = c("hg19", "hg38"))
```

## Arguments

- version:

  str. hg19 or hg38. Default: hg19

## Value

tibble

## Details

hg19:
https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
hg38:
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
