# metric of pairwise differences between two cells.

Inspired by MSKCC SPECTRUM paper. Bins are aligned between two cells and
marked for if they have the same state Then segments of matching and
non-matching runs of bins are found (filtering those smaller than a
specified minimum). These segments are then re-split into 500kb bins and
the difference becomes the number of matching bins divided by the number
of considered bins.

## Usage

``` r
pairwise_bin_difference(
  bin_df,
  cells = c(),
  min_seg_length = 2500000,
  return_pairs_matrix = FALSE
)
```

## Arguments

- bin_df:

  a dataframe of read bins with states. Expected columns of: cell_id,
  chr, start, end, state

- cells:

  optional vector specifying cells to compare. If it's blank, all cells
  are compared. If it's 1 cell, then that one cell is compared to all
  others. If it's 2 or more, then just the specified cells are compared
  to each other.

- min_seg_length:

  double. This is the minium length of matching segment bins to use when
  measuring similarity.

- return_pairs_matrix:

  boolean. If TRUE, returns a pairwise matrix object of distances. This
  is useful to then pass to functions like hclust() and so forth. Can
  also do afterwards with dlptools::convert_dists_to_pairwise()

## Value

tibble of cell pairs and metrics about their differences.

## Details

This function is slow, and the number of pairwise comparisons grows
quickly. Dramatic speed improvements can be had by setting up a parallel
plan for furrr like so:

future::plan(future::multicore, workers=N_CORES_YOU_WANT)

Example for 100 cells, which is 4950 pairs, this function will take 4
minutes with 4 cores.

The returned DF is organized by each cell and the distances to each
other cell (so there are some redundant comparisons, like cell 1 vs cell
2 and cell 2 vs cell 1). There is also a column "nearest_neighbour"
which is a boolean identifying which comparison is the minimum distance
for each cell.
