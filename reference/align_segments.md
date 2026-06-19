# cut segments to equal widths across cells/samples.

Segments are continuous stretches of a chromosome with the same copy
number. But, generally, these will be different across cells/samples.
This function will create segments in the smallest consistent intervals
possible across the samples. See example, but it's basically equalizing
the segment widths among all of the cells, so copy number is measured in
the same "bin" across cells.

## Usage

``` r
align_segments(segs_df, chrom_col = "chr", cell_col = "cell_id")
```

## Arguments

- segs_df:

  dataframe. Copy number segments of cells/samples. Need chromosome,
  start, end.

- chrom_col:

  string. name of column containing chromosome information

- cell_col:

  string. name of column with cell/sample IDs

## Value

tibble of broken down segments with same input columns + seg_width of
the new segments.

## Details

The main value of this is that some tools (medicc, or plotting) want the
same intervals across all samples (i.e., bins of equal width). But going
all the way back to 500 kb bins (i.e., typical DLP) can be excessive
when there are longer segments within cells. Cutting segments into their
minimal spans will generally be a smaller set of data than a full 500kb
bin breakdown (unless CN data is highly noisy among samples).

## Examples

``` r
sgs <- tibble::tibble(
  cell_id = c("A", "A", "B", "B"),
  chr = rep("chr1", 4),
  start = c(1, 11, 1, 8),
  end = c(10, 25, 7, 25),
  state = c(2, 4, 3, 8)
)

sgs
#> # A tibble: 4 × 5
#>   cell_id chr   start   end state
#>   <chr>   <chr> <dbl> <dbl> <dbl>
#> 1 A       chr1      1    10     2
#> 2 A       chr1     11    25     4
#> 3 B       chr1      1     7     3
#> 4 B       chr1      8    25     8

align_segments(sgs)
#> # A tibble: 6 × 6
#>   cell_id chr   start   end state seg_width
#>   <chr>   <fct> <int> <int> <dbl>     <int>
#> 1 A       chr1      1     7     2         6
#> 2 A       chr1      8    10     2         2
#> 3 A       chr1     11    25     4        14
#> 4 B       chr1      1     7     3         6
#> 5 B       chr1      8    10     8         2
#> 6 B       chr1     11    25     8        14
```
