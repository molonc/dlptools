# extract copy state values.

This function is really only here to offer the counting in pre-set
categories, and provide alignment with other feature extraction types.
It is mostly a pointless function and you should probably do this
yourself.

## Usage

``` r
extract_cn(
  segs_df,
  sample_col = "cell_id",
  state_col = "state",
  return = c("values", "counts")
)
```

## Arguments

- segs_df:

  dataframe. copy number segments for samples.

- sample_col:

  string. Name of column with CN state values

- return:

  string. "values" (default) or "counts". Values are the observed values
  for cells, counts are the counts of these values in pre-determined
  categories.

## Details

Also, it is a question if copy state should even be used in signature
fitting, as this can occassionally lead to the same process being
artifically split by ploidy of samples. See [Drews et
al.](https://www.nature.com/articles/s41586-022-04789-9) for a
discussion of this idea.

Predefined categories for summarizing counts are: CN = 1, 2, 3, 4, 5, 6+
