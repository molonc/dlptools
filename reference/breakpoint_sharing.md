# summarize pairwise breakpoint sharing across cells

Will summarize the number of shared breakpoints (changes in copy number)
for all pairs of cells. Alternatively, will return jaccard similarity
between all pairs (intersection of shared breakpoints over the union of
breakpoints).

Returns a pairwise matrix that can be used for various downstream
things, like clustering.

## Usage

``` r
breakpoint_sharing(seg_df, jaccard = FALSE)
```

## Arguments

- seg_df:

  dataframe. Segmented copy number profiles per cell.

- jaccard:

  bool. True will calculate and return jaccard similarity.

## Value

pairwise matrix of the number of shared breakpoints between each pair of
cells, or the jaccard similarity.
