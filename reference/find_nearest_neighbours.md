# find the nearest neighbour of each cell

After running
[`pairwise_bin_difference()`](https://molonc.github.io/dlptools/reference/pairwise_bin_difference.md),
we can feed that output dataframe to this function to find the nearest
neighbour of each cell. The pairwise function does the comparisons in a
non-redundant way (A-B, A-C, B-C), leaving out some redundant
comparisons (i.e., C-A, C-B). If A's nearest neighbour is C, it's not
guaranteed that C's nearest neighbour is A. This function takes a cell
specific focus and finds the nearest neighbour of each.

## Usage

``` r
find_nearest_neighbours(pairwise_diffs)
```

## Arguments

- pairwise_diffs:

  tibble output of
  [`pairwise_bin_difference()`](https://molonc.github.io/dlptools/reference/pairwise_bin_difference.md)

## Value

tibble/dataframe

## Details

As with
[`pairwise_bin_difference()`](https://molonc.github.io/dlptools/reference/pairwise_bin_difference.md),
using
[`furrr::future_map()`](https://furrr.futureverse.org/reference/future_map.html)
internally, so speed improvements with a
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
being set.
