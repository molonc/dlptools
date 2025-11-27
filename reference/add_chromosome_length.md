# add chromosome lengths to a cn dataframe

Handles chromosomes specified with "chr#" or just \#/X/Y.

## Usage

``` r
add_chromosome_length(cn_df, version = c("hg19", "hg38"), chrom_col = "chr")
```

## Arguments

- cn_df:

  dataframe of copy number information

- version:

  string. hg19 (default) or hg38

- chrom_col:

  string. name of column with chromosome information.
