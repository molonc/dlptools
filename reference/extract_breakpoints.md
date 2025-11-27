# extracting counts of breakpoints per user-defined scope

Counting the number of breakpoints (i.e., transitions between copy
number segments) per arm or per window (typically 10Mb).

## Usage

``` r
extract_breakpoints(
  segs_df,
  scope = c("windows", "arms"),
  genome_version = c("hg19", "hg38"),
  window_size = 1e+07,
  sample_col = "cell_id",
  chrom_col = "chr",
  return = c("values", "counts")
)
```

## Arguments

- segs_df:

  dataframe. Copy number segments for samples

- scope:

  string. "windows" (default) or "arms", i.e., what to target for
  counting

- genome_version:

  string. "hg19" (default) or "hg38"

- window_size:

  integer. How big of a window to use, if extracting counts per a
  "windows" scope. Most publications use 10Mb, which is the default
  (10e6)

- sample_col:

  string. Name of column with cell/sample names

- chrom_col:

  string. Name of column with chromosome names

- return:

  string. "values" (default) or "counts". Values are the observed values
  for cells, counts are the counts of these values in pre-determined
  categories.

## Value

dataframe. Sample IDs and the observed breakpoint counts per scope.

## Details

Optional summarizing in pre-defined categories of:

Arm breakpoints: 0, 1, 2, 3, 4, 5, 6+ Window breakpoints: 0, 1, 2, 3+
