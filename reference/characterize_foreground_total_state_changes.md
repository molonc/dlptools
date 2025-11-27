# mark the types of change between a parent node state and the tip.

This compares the tip and parent state columns to evaluate if a
foreground change has occurred, what type of change, and how much of a
change.

## Usage

``` r
characterize_foreground_total_state_changes(
  bin_states_df,
  state_col = "state",
  parent_state_col = "parent_state"
)
```

## Arguments

- bin_states_df:

  bin based dataframe with parent and tip state columns.

- state_col:

  name of tip state column

- parent_state_col:

  name of parent state column.

## Value

input dataframe with new columns

## Details

e.g., if parent state is 5, tip is 4, it will get marked as a
'foreground' event, fg-gain, with a change of 1.

4 coulmns are added: foreground: boolean, if parent and tip states are
different fg_type: character, fg-gain, fg-loss, or NA fg_change:
numeric, tip - parent state. Positive indicates gain in tip. background:
boolean, inverse of foreground
