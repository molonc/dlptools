# Add back chr, start, end of bins missing from cells

Often we filter read bins based on various criteria, or tools we use
might drop bins of read/state data. This function will put missing bins
back into a dataframe for cells.

## Usage

``` r
add_missing_bins_for_cells(
  state_df,
  version = c("hg19", "hg38"),
  bin_size = 5e+05,
  cell_metadata_cols = c(),
  input_chroms_only = TRUE
)
```

## Arguments

- state_df:

  the dataframe/tibble to insert the missing bins into

- version:

  which genome version you are using. Default "hg19", or select "hg38"

- bin_size:

  size of bins that should be there.

- cell_metadata_cols:

  vector of columns to carry through to inferred bins

- input_chroms_only:

  bool. Only insert bins for chromosomes in the input state dataframe.
  I.e., if you've filtered X or Y, don't add bins for X or Y. Setting to
  FALSE adds bins for entire genome. Default: TRUE

## Value

input dataframe with extra bins for each cells that were missing.

## Details

Basic operation will result in NAs for every column except cell_id, chr,
start, end for the newly added bins. But you can also carry over cell
metadata by specifying which columns to keep for the inferred bins. This
should be cell level data.
