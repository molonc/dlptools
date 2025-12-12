# create a tile plot of read state calls.

builds a basic ggplot with geom_tile on a reads df. Only expects columns
of cell_id, start, state, chr

## Usage

``` r
plot_read_bins_basic(reads_df)
```

## Arguments

- reads_df:

  a table of reads data

## Value

ggplot object
