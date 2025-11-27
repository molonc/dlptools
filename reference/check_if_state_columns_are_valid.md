# internal check for using the characterization function

To characterize forground changes, need wither numeric columns or
factors with a defined order to them.

## Usage

``` r
check_if_state_columns_are_valid(
  states_df,
  state_col = "state",
  parent_state_col = "parent_state"
)
```
