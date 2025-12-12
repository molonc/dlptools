# Find outlier cells using a beta distribution

Also inspired by MSKCC SPECTRUM paper. Using the measured pairwise bin
distances between cells `pairwise_bin_differences()`, this function
takes the nearest neighbour of each cell and fits a beta distribution.
Then using this distribution, it finds outlier cells based on a selected
percentile of the distribution (default 99th).

## Usage

``` r
find_outlier_cells(pairwise_diffs, nn_cells, outlier_percentile = 0.99)
```

## Arguments

- nn_cells:

  dataframe of
  [`find_nearest_neighbours()`](https://molonc.github.io/dlptools/reference/find_nearest_neighbours.md)

- outlier_percentile:

  double. Default 0.99. What percentile of the distribution to consider
  an outlier cell.

- cell_diffs:

  dataframe of differences from `pairwise_bin_differences()`

## Value

NA or tibble of information on cells considered outliers. NA if no
outliers found.
