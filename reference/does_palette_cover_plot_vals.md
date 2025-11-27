# internal to confirm color choice was ok

internal to confirm color choice was ok

## Usage

``` r
does_palette_cover_plot_vals(p_vals, palette, exclude_na = TRUE)
```

## Arguments

- p_vals:

  vector of values to be plotted in the heatmap with the palette

- palette:

  palette being used. A named vector with names corresponding to the
  plot values.

- exclude_na:

  bool. Whether to drop NAs from the check. ComplexHeatmap handles NAs
  with a default color regardless of palette.
