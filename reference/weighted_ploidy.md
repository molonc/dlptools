# calculate sample ploidy with a weighted CN mean

CN is first weighted by the length of its segment, a weighted mean of
these segments for each chromosome is made, then the mean across all
chrosomomes is made for each sample.

## Usage

``` r
weighted_ploidy(
  cn_df,
  sample_col = "cell_id",
  start_col = "start",
  end_col = "end",
  cn_col = "state",
  chrom_col = "chr"
)
```

## Arguments

- cn_df:

  dataframe input of segmented (or not) sample copy numbers

- sample_col:

  string column name identifying samples

- start_col:

  string column name identifying start of segments

- end_col:

  string column name identifying end of segments

- cn_col:

  string column name for copy number states

- chrom_col:

  string column name for chromosomes

## Value

tibble of per-sample ploidy estimates

## Details

In this way, the CN contribution to the ploidy is based on it's length,
and each chromosome contributes equally to the overall sample ploidy.
