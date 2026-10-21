# break a chromosome up into intervals of a defined window size

break a chromosome up into intervals of a defined window size

## Usage

``` r
create_chrom_window_intervals(
  window_size = 5e+05,
  genome_version = c("hg19", "hg38"),
  return_type = c("granges", "tibble"),
  cap_at_chromosome_lengths = FALSE
)
```

## Arguments

- window_size:

  integer. The size of window to split the chromosome into.

- genome_version:

  string. "hg19" (default) or "hg38"

- cap_at_chromosome_lengths:

  bool. TRUE will make sure the last bin on a chromsome ends at the real
  lenght of a chromosome, and not exceed it. Default is FALSE to produce
  evenly sized bins.

## Value

list. Named by chromosome, vectors of window starts.
