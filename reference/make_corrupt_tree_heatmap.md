# create a complex heatmap object with a tree.

Basically creates an empty heatmap with a tree added as a left side
annotation.

## Usage

``` r
make_corrupt_tree_heatmap(
  phylo_tree,
  clones_df = NULL,
  color_clones = FALSE,
  clone_palette = NULL
)
```

## Arguments

- phylo_tree:

  a phylo object of a tree.

- clones:

  (optional) adds clones as OTUs, will eventually allow for color

## Value

a ComlexHeatmap heatmap object with a tree.
