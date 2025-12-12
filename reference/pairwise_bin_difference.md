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
  targ_cells = c(),
  min_seg_length = 2500000,
  return_pairs_matrix = FALSE
)
```

## Arguments

- bin_df:

  a dataframe of read bins with states. Expected columns of: cell_id,
  chr, start, end, state

- targ_cells:

  optional vector specifying which cells to compare to all other cells.

- min_seg_length:

  double. This is the minium length of matching segment bins to use when
  measuring similarity.

## Value

tibble of cell pairs and metrics about their differences.

## Details

This function is slow, and the number of pairwise comparisons grows
quickly. Dramatic speed improvements can be had by setting up a parallel
plan for furrr like so:

future::plan(future::multisession, workers=N_CORES_YOU_WANT)

Example for 100 cells, which is 4950 pairs, this function will take 4
minutes with 4 cores.

These comparisons are unique pairs! So if 3 input cells: A, B, C,
comparisons made are A-B, A-C, B-C. So one input cell will always be
"missing" from the index column. See [the
vignette](https://molonc.github.io/dlptools/articles/pairwise-differences.html)
for how you can turn this around, but it can create a very large DF with
redundant comparisons to rearrange to have each cell against all others.
