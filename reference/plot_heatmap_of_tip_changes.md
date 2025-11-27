# make a plot of just the foreground events

Colors highlight the change in state between the tip and parent node.

## Usage

``` r
plot_heatmap_of_tip_changes(
  states_df,
  phylogeny,
  file_name,
  changes_col = "fg_change"
)
```

## Arguments

- states_df:

  bin based dataframe of state changes. Typically the output of
  characterize_foreground_total_state_changes() or from
  medicc_profiles_to_foreground(). But really, any bin dataframe with a
  column indicating change between tip and parent will work.

- phylogeny:

  phylo object of the tree being used for analysis.

- file_name:

  what to call the resulting plot png.

- fg_change:

  the column to plot with the changes between parent and tips.
