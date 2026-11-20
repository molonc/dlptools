# load and prep default centromerefiles from UCSC

load and prep default centromerefiles from UCSC

## Usage

``` r
load_ucsc_centromeres(
  version = c("hg19", "hg38"),
  source = c("cytoband", "gap")
)
```

## Arguments

- version:

  string. genome version hg19 or hg38. Default: hg19

- source:

  string. cytoband or gap. Where do you want the centromere information
  to come from, UCSC cytband or gap tables?

## Value

tibble of parsed centromere spans
