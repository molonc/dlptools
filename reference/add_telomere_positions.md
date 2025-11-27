# add telomere positions to a CN df

See
[`import_telos_file()`](https://molonc.github.io/dlptools/reference/import_telos_file.md)
but really just a convenience function. Telomere positions in the UCSC
gap file are all just 10Kb, except chr17, which I set to 3Kb and 5Kb,
again see linked function for details.

## Usage

``` r
add_telomere_positions(cn_df, version = c("hg19", "hg38"))
```

## Arguments

- cn_df:

  dataframe of CN information, segments or read bins. Just needs a 'chr'
  column.

- version:

  string. hg19 (default) or hg38, not that it makes a difference.

## Value

input dataframe with new columns of telomere information
