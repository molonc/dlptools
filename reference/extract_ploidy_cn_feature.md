# Count changes of state relative to ploidy.

Marks CN segments as a gain or loss, relative to the mode ploidy of the
sample. Internally using
[mark_cn_relative_to_ploidy](https://molonc.github.io/dlptools/reference/mark_cn_relative_to_ploidy.md).
See that function for argument details.

## Usage

``` r
extract_ploidy_cn_feature(
  segs_df = NA,
  sample_col = "cell_id",
  annotate_input = FALSE,
  return_matrix = FALSE
)
```

## Arguments

- segs_df:

  dataframe. CN segments

- sample_col:

  string. Name of the column with cell_id/other sample name

- annotate_input:

  boolean. return input dataframe annotating each

- return_matrix:

  boolean. Return a cell-by-feature matrix of counts.

## Value

tibble/dataframe of counts
