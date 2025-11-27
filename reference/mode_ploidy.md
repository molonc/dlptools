# find the mode CN per chromosome, then mode across the chromosomes

Done per chromosome first so that large chromosomes don't dominate the
result. It is critical that the input be even bin level data, or this
makes no sense to do. Measures of the copy number states need to be done
in evenly sized bins.

## Usage

``` r
mode_ploidy(
  bin_df,
  sample_col = "cell_id",
  cn_col = "state",
  chrom_col = "chr"
)
```

## Arguments

- bin_df:

  dataframe of bin level data

- sample_col:

  string column name identifying samples

- cn_col:

  string column name for copy number states

- chrom_col:

  string column name for chromosomes

## Value

tibble/dataframe of results by cell_id/sample
