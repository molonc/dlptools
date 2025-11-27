# get distance of states to siblings for a tree tip

Computes a basic string distance between a tip and all siblings in a
tree. A tip can have one or more siblings, and the mean distance to them
will be returned.

## Usage

``` r
get_dist_to_sibs(
  tip_label,
  states_string,
  tree,
  state_str_df,
  tip_label_col = "tip_label"
)
```

## Arguments

- states_string:

  string of the states for a tree tip

- tree:

  the full tree of data

- state_str_df:

  dataframe of all tips and the states strings

- tip_label_col:

  name of the column in the state_str_df containing the names of the
  tips on the tree

- tip:

  label of a tree tip

## Details

e.g., in tree (A, (B, C)), B has one sibling C, but for A both B & C are
siblings. So the mean distance of A and C and B is returned.

This function is intended to be used through it's wrapper
dlptools::compute_all_tip_sibling_distances()
