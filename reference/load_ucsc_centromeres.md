# load and prep default centromerefiles from UCSC

Cytoband files originally obtained from UCSC golden path genomes and
parsed to collect the total centromere spans.

## Usage

``` r
load_ucsc_centromeres(version = c("hg19", "hg38"), centro_file = NULL)
```

## Arguments

- centro_file:

  NULL or string to the path if you download yourself

- hg19:

  boolean to target hg19 for loading

- hg38:

  boolean to target hg38 for loading

## Value

tibble of parsed centromere spans

## Details

hg19:
https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
hg38:
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz

TODO: can I stash the file if it's already loaded once?
