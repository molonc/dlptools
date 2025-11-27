# Plot metrics across the chip

Map of a metric as observed in the layout of the DLP chip.

## Usage

``` r
chip_plot(metrics_df, targ_val)
```

## Arguments

- metrics_df:

  dataframe of metrics for a run.

- targ_val:

  the value of interest for plotting. E.g., experimental_condition,
  quality, total_reads.

## Value

ggplot2 object
