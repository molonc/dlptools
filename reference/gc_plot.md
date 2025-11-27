# visualize GC correction for a cell

Seeing how GC correction worked for a cell can give an indication if the
multiplier chosen by HMMCopy was appropriate to form integer copy
numbers. Also can just be informative to see what's going on behind the
float copy number to integer CN value.

## Usage

``` r
gc_plot(reads_df, cellid, plot_choice = c("both", "raw", "corrected"))
```

## Arguments

- reads_df:

  dataframe of reads like information.

- cellid:

  string of cell to plot

- plot_choice:

  string of plots to show. Default is both, alternatively " 'raw' for
  just the raw GC values or 'corrected' for just the corrected values.

## Value

ggplot2 object
