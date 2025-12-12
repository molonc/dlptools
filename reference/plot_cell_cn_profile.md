# plot copy number profile plot of a cell.

Basically a geom_point of the float copy number, colored by the integer
state number. singals::plotCNprofile() is similar, possibly with more
features. But this function works in a pinch.

## Usage

``` r
plot_cell_cn_profile(
  reads_df,
  cell_ids = c(),
  pseduo_log_y = FALSE,
  yaxis = "copy"
)
```

## Arguments

- reads_df:

  dataframe of read bins like information.

- cell_ids:

  string/vector of cell ids to plot. Remember, too many and the plot
  will be useless. If blank, will plot all cells in dataframe, but only
  if less than 10, because honestly at that point, you probably need a
  different plot.

- yaxis:

  string. Defaults to copy, which is typically what is wanted for
  plotting. But you can change to other bin values, like state.

- pseudo_log_y:

  boolean. Pseudo-log transform the y axis, as sometimes high copy
  numbers obscure the plot. Or you can always set a
  ggplot2::coord_cartesian() yourself on the plot returned by this
  function.

## Value

ggplot2 object
