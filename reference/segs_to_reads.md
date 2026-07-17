# Chop Genomic Segments into Fixed-Size Genomic Bins

Splits variable-length copy number segments into fixed-size genomic bins
(windows) of a specified width (e.g., 500 Kb).

## Usage

``` r
segs_to_reads(
  segs_df,
  bin_size = 5e+05,
  genome_version = c("hg19", "hg38"),
  return_type = c("tibble", "granges"),
  seg_start_col = "start",
  seg_end_col = "end",
  sample_id_col = "cell_id",
  state_col = "state",
  other_meta_cols = c(),
  chrom_col = "chr"
)
```

## Arguments

- segs_df:

  A data.frame or tibble containing genomic segments with start, end,
  and chromosome name columns.

- bin_size:

  Numeric. The width of the genomic windows in base pairs. Defaults to
  `5e5` (500 Kb).

- genome_version:

  Character. The assembly build to use for creating genomic windows.
  Must be one of `"hg19"` or `"hg38"`.

- return_type:

  Character. The format of the returned object. Must be one of
  `"granges"` or `"tibble"`. Defaults to `"tibble"`.

- seg_start_col:

  Character. The name of the column in `segs_df` representing segment
  start coordinates. Defaults to `"start"`.

- seg_end_col:

  Character. The name of the column in `segs_df` representing segment
  end coordinates. Defaults to `"end"`.

- sample_id_col:

  Character. The name of the column in `segs_df` identifying the sample
  or cell. Defaults to `"cell_id"`.

- state_col:

  Character. The name of the column in `segs_df` representing the copy
  number state. Defaults to `"state"`.

- other_meta_cols:

  Character vector. Optional additional metadata column names in
  `segs_df` to carry over to the output. Defaults to an empty vector.

- chrom_col:

  Character. The name of the column in `segs_df` representing
  chromosomes. Defaults to `"chr"`.

## Value

Depending on `return_type`:

- `"granges"`: A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object containing the chopped bins with associated metadata columns.

- `"tibble"`: A
  [`tibble`](https://tibble.tidyverse.org/reference/tibble.html) version
  of the GRanges object, with the chromosome column renamed back to
  `chrom_col`.

In both formats, the output retains the original segment start and end
positions in the `seg_start` and `seg_end` metadata columns.
