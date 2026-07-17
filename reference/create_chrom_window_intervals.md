# break a chromosome up into intervals of a defined window size

break a chromosome up into intervals of a defined window size

## Usage

``` r
create_chrom_window_intervals(
  window_size = 1e+07,
  genome_version = c("hg19", "hg38"),
  return_type = c("granges", "tibble")
)
```

## Arguments

- window_size:

  integer. The size of window to split the chromosome into.

- genome_version:

  string. "hg19" (default) or "hg38"

## Value

list. Named by chromosome, vectors of window starts.
