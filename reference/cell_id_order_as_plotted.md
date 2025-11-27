# grab cell ids in the order that they are plotted

This is to align state calls and other things. It will make a ggplot of
the tree the same as the heatmap code does, and then pull out the cell
ID in the order that they are plotted

## Usage

``` r
cell_id_order_as_plotted(phylo_tree)
```

## Arguments

- phylo_tree:

  a phylo object of the tree being used in the heatmap

## Value

vector of cell_ids (or really, whatever are the tip labels)
