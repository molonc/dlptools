# add parent node labels to a dataframe of cell ids with states

add parent node labels to a dataframe of cell ids with states

## Usage

``` r
add_tip_ancestors_to_df(states_df, phylogeny)
```

## Arguments

- states_df:

  bin based dataframe of cell ids and their states

- phylogeny:

  phylo object

## Value

input dataframe, but with parent node information added. Creates a
'parent_node' column.
