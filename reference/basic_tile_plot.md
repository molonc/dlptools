# create a tile plot of read state calls.

builds a basic ggplot with geom_tile on a reads df. Only expects columns
of cell_id, start, state, chr

## Usage

``` r
basic_tile_plot(reads_df)
```

## Arguments

- reads_df:

  a table of reads data (e.g., could load with
  [`import_dlp_files()`](https://molonc.github.io/dlptools/reference/import_dlp_files.md))

## Value

ggplot object
