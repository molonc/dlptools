# calculate jaccard similarity using breakpoint matrices

Jaccard is the number of shared breakpoints (intersection), over the
total number of unique breakpoints (union) between two cells.

## Usage

``` r
calc_jaccard_similarity(bps_mtx, shared_bps)
```

## Arguments

- bps_mtx:

  sparse matrix. Breakpoint presence in cells

- shared_bps:

  sparse matrix. the numbers of shared breakpoints between cells.

## Value

sparse matrix of jaccard similarity
