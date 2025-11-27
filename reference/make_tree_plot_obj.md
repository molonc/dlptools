# create a ggplot object of a phylogenetic tree

Will make a plot of a tree that can be added to a heatmap object. The
plot itself is just a simple tree with no tip labels.

## Usage

``` r
make_tree_plot_obj(
  phylo_tree,
  clones_df = NULL,
  clone_palette = NULL,
  color_clones = FALSE
)
```

## Arguments

- phylo_tree:

  a phylo object

- clones:

  (optional) adds clones as OTUs (TODO is to allow this to then make
  colors)

## Value

ggplot object
