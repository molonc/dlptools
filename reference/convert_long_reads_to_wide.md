# convert long format reads to wide format

A common manipulation with reads files for various analyses is to
reshape long format reads data (each row is a 500kb bin with state
values for each cell) to wide format, with chr_start_end rows and
cell_id columns.

## Usage

``` r
convert_long_reads_to_wide(reads_df, state_col = "state")
```

## Arguments

- reads_df:

  is the reads table to convert.

## Value

wide format table

## Details

minimal required columns for input are: chr,start,end,cell_id,state
