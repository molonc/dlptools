# convert cell distances to a pairwise matrix

Takes the output of dlptools::pairwise_bin_difference() and returns a
pairwise matrix of the distances. Useful to then pass to functions like
stats::hclust() and related ideas.

## Usage

``` r
convert_dists_to_pairwise(cell_dists)
```

## Arguments

- cell_dists:

  dataframe from dlptools::pairwise_bin_difference()

## Value

matrix of index vs comp cell distances.
