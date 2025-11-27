# sizes of the segments

This function is as simple as it sounds, end - start.

## Usage

``` r
extract_segment_sizes(
  segs_df,
  sample_col = "cell_id",
  return = c("values", "counts")
)
```

## Arguments

- segs_df:

  dataframe. copy number segments for samples.

- sample_col:

  string. Name of column with cell/sample names

- return:

  string. "values" (default) or "counts". Values are the observed values
  for cells, counts are the counts of these values in pre-determined
  categories.

## Value

dataframe. sample ids and all observed segment sizes.

## Details

Used as a setup for extracting process based features *á la*:

- Macintyre et al. 2018

- Drews et al. 2018

Really, something like this should be used to generate values that then
you define categories for to count occurrences.

Related functions include: extract_extract_changepoint and
[extract_breakpoints](https://molonc.github.io/dlptools/reference/extract_breakpoints.md)

Can also summaries counts of pre-defined categories. With categories
being:

\< 5mb, 5-10 Mb, 20-50 Mb, 50-100 Mb, and 100+Mb.
