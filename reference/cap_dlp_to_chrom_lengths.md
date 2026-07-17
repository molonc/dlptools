# cap DLP data at the lengths of chromosomes

DLP can output bins/segments longer than the actual chromosomes as it
chops to 500kb bins. This function sets the `end` column to be no longer
than the chromosome it occurs on.

## Usage

``` r
cap_dlp_to_chrom_lengths(
  rs_df,
  chrom_col = "chr",
  genome_version = c("hg19", "hg38")
)
```

## Arguments

- rs_df:

  tibble/dataframe. Reads or segments dataframe

- chrom_col:

  string. Name of chromosome column

- genome_version:

  string. Name of genome version to use. Default: hg19

## Value

input dataframe with final bins/segments not exceeding total chromosome
lengths
