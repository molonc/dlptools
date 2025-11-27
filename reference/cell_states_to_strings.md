# collapse cell states to strings

Converts a long format states dataframe into 1 row per cell with all
states converted to a letter string.

## Usage

``` r
cell_states_to_strings(
  states_df,
  states_col = "state",
  cell_id_col = "cell_id"
)
```

## Arguments

- states_df:

  a dataframe with states for each bin for cells (or whatever the tree
  tips are)

- states_col:

  name of the column with state values

- cell_id_col:

  name of the column with the cell ids (tree tip names)
