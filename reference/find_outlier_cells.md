# Find outlier cells using a beta distribution

Also inspired by MSKCC SPECTRUM paper. Using the measured pairwise bin
distances between cells (dlptools::pairwise_bin_differences()), this
function takes the nearest neighbour of each cell and fits a beta
distribution. Then using this distribution, it finds outlier cells based
on a selected percentile of the distribution (default 99th).

## Usage

``` r
find_outlier_cells(cell_diffs, outlier_percentile = 0.99)
```

## Arguments

- cell_diffs:

  the dataframe of differences from dlptools::pairwise_bin_differences()

- outlier_percentile:

  double. Default 0.99. What percentile of the distribution to consider
  an outlier cell.

## Value

NA or tibble of information on cells considered outliers. NA if no
outliers found.
