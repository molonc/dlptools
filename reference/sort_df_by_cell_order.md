# sort a table given a vector of cell_ids

Typically used to sort a dataframe based on the plotted tip order to
align states heatmap/annotations/clone IDs to the plotted tree

## Usage

``` r
sort_df_by_cell_order(targ_df, cell_order)
```

## Arguments

- targ_df:

  a table with cell_ids to sort

- cell_order:

  a vector of cell_ids in the desired order (e.g., pulled from a ggplot
  of a tree)

## Value

table that has been sorted
