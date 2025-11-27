# internal workhorse of pairwise_bin_difference function

internal workhorse of pairwise_bin_difference function

## Usage

``` r
compare_two_cells(cell_1_df, cell_2_df, min_seg_length)
```

## Arguments

- cell_1_df:

  dataframe of bin state data for one cell

- cell_2_df:

  dataframe of bin state data for the other cell

- min_seg_length:

  int/double. minimum length of matching/non-matching runs of bins to
  consider in the calculation of the differences between the cells.

## Value

1 row tibble summarizing differences.
