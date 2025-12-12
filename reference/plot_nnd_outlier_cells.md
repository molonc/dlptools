# visualize cell nearest neighbour distance values

Generates a simple plot of the nearest neighbour distance for each cell,
which is highlighted with the point. Outlier cells are indicated by
colour. A violin plot is used to show distance to all other cells.

## Usage

``` r
plot_nnd_outlier_cells(pairwise_diffs, nn_cells, outlier_cells)
```

## Arguments

- pairwise_diffs:

  dataframe of
  [`pairwise_bin_difference()`](https://molonc.github.io/dlptools/reference/pairwise_bin_difference.md)

- nn_cells:

  dataframe of
  [`find_nearest_neighbours()`](https://molonc.github.io/dlptools/reference/find_nearest_neighbours.md)

- outlier_Cells:

  dataframe of
  [`find_outlier_cells()`](https://molonc.github.io/dlptools/reference/find_outlier_cells.md)

## Value

ggplot object
