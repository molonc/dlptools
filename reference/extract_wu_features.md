# Extract CN features following Wu et al

This function extracts copy number features in the style of the paper:

## Usage

``` r
extract_wu_features(
  segs_df,
  sample_col = "cell_id",
  state_bin_max = 5,
  bin_breaks = NA,
  annotate_input = FALSE,
  return_matrix = FALSE,
  ...
)
```

## Arguments

- segs_df:

  dataframe. CN segments

- sample_col:

  string. Name of the column with cell_id/other sample name

- state_bin_max:

  int. Maximum CN to consider for bins. All CNs of this value and higher
  are grouped together. Default of 5 follows paper.

- bin_breaks:

  floats, how to break up segment sizes. Bins will be one more than
  breaks. Defaults follow paper. Default is \< 5Mb, 5–10Mb, \> 10Mb
  specified as c(5e6, 10e6 + 1). Internally, base::cut() is used, so 2
  splits produces 3 bins.

- annotate_input:

  boolean. return input dataframe annotating each segment with the
  feature categories it falls into.

- return_matrix:

  boolean. Return a cell-by-feature matrix of counts.

- ...:

  can pass change_split_val to alter critical value for AA/BB split

## Value

default return is a tibble of feature counts for each cell id.

## Details

Wu et al. 2025. Single-cell copy number alteration signature analysis
reveals masked patterns and potential biomarkers for cancer. bioRxiv.

https://www.biorxiv.org/content/10.1101/2025.03.02.641098v1

They employ 4 base features, that they cross for 90 categories:

1.  CN states: 5 bins, 0-1, 2, 3, 4, 5+

2.  segment size: 3 bins, \<5 MB, 5-10Mb, 10Mb,

3.  segment shape: 3 bins, LL (low left, low right segment), HH (high
    left, high right segment), OT (other)

4.  segment change: 2 bins: AA (difference between surrounding segments
    \<= some critical value), or BB

Some issues:

- whole chromosome amplifications/losses not captured

- chromosomes need to have at least 3 changes to be truly reflected in
  these categories. Those with fewer are backfilled based on hardcoded
  rules

- for segment change, the paper says they considered changes \> 2 as BB,
  but in their actual code is set to 1. This function follows their
  code, but can be altered.
