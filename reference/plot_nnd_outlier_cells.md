# visualize cell nearest neighbour distance values

Generates a simple plot of each cell and highlights the min NN-distance
for each. Also highlights cells that are considered outliers.

## Usage

``` r
plot_nnd_outlier_cells(cell_diffs, outlier_cells)
```

## Arguments

- cell_diffs:

  output dataframe of dlptools::pairwise_bin_difference()

- outlier_Cells:

  output dataframe of dlptools::find_outlier_cells()

## Value

ggplot object
