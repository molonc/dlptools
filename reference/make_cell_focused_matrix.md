# mostly internal for rearranging the pairwise DF to focus on a cell.

This is useful because with non-redundant pair comparisons, the "last
cell" in the list is never the "index_cell", as it has already been
compared ot all others. But sometimes useful to still see it.

## Usage

``` r
make_cell_focused_matrix(in_df, cell)
```

## Arguments

- in_df:

  dataframe most likely from
  [`pairwise_bin_difference()`](https://molonc.github.io/dlptools/reference/pairwise_bin_difference.md)

- cell:

  string. Target cell ID to put as index cell

## Value

tibble

## Details

e.g., cells: A, B, C comparisons would be: A-B, A-C, B-C generate C
focus with this function: C-A, C-B or B focus: B-A, B-C
