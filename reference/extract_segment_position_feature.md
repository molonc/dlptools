# count the segment-span-on-chromosome event types.

Critical to this function is
[`mark_segs_chromosome_span()`](https://molonc.github.io/dlptools/reference/mark_segs_chromosome_span.md).
It is important to read and understand that function and its arguments.

## Usage

``` r
extract_segment_position_feature(
  segs_df,
  sample_col = "cell_id",
  annotate_input = FALSE,
  return_matrix = FALSE,
  ...
)
```

## Arguments

- segs_df:

  dataframe of CN segments

- sample_col:

  string. Name of the column with cell_id/other sample name

- annotate_input:

  boolean. return input dataframe annotating each segment.

- return_matrix:

  boolean. Return a cell-by-feature matrix of counts.

## Value

tibble/dataframe of counts

## Details

This function basically just calls
[`mark_segs_chromosome_span()`](https://molonc.github.io/dlptools/reference/mark_segs_chromosome_span.md)
and summarizes the results. Arguments can be passed to that underlying
function. Passing no arguments means you are happy with the defaults.
See
[`mark_segs_chromosome_span()`](https://molonc.github.io/dlptools/reference/mark_segs_chromosome_span.md)
to understand what the defaults are.
