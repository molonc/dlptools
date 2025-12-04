# convert reads table to segments table

Often the segments file given by the DLP pipeline is not what you want.
What you likely want is to do various filtering at the reads level, and
then make a segments file from that. This function will take a reads
table and covert it to a segments table.

## Usage

``` r
reads_to_segs(reads_df)
```

## Arguments

- reads_df:

  a table with standard reads data

## Value

tibble/dataframe with read bins organized into segment blocks.
