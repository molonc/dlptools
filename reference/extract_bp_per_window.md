# convenience function to extracting breakpoints per window size

just calling the generic
[extract_breakpoints](https://molonc.github.io/dlptools/reference/extract_breakpoints.md),
with some pre-loaded options. See that function for details.

## Usage

``` r
extract_bp_per_window(
  segs_df,
  window_size = 1e+07,
  sample_col = "cell_id",
  chrom_col = "chr",
  genome_version = c("hg19", "hg38"),
  return = c("values", "counts")
)
```

## Value

dataframe. Sample IDs and the observed breakpoint counts observed in the
specified window size.
