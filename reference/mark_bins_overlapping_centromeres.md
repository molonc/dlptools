# add boolean if bin overlaps with a centromere.

Can optionally specify a padding to mark locations of bins as being
close enough to a centromere. Often bins near centromeres are corrupt in
their state calls. In the past, we have filtered within 3 Mb of
centromeres.

## Usage

``` r
mark_bins_overlapping_centromeres(
  reads_df,
  padding = 0,
  bin_start_col = "start",
  bin_end_col = "end",
  version = c("hg19", "hg38")
)
```

## Arguments

- reads_df:

  tibble of read data

- padding:

  int of number of BP to add to each side of the centromere

- bin_start_col:

  Default: start. column name of the start of bins.

- bin_end_col:

  Default: end column name of the end of bins.

- version:

  default 'hg19', or choose 'hg38' for locations of centromeres.

## Value

input table, but with a boolean 'within_centro' column added (and
potentially other centromere information columns, if needed)

## Details

IMPORTANT! This adds the padding to each side of the centromere. So if
you specify a padding of 3 Mb, it will be within 3 Mb of the start and 3
Mb of the end of the centromere.
