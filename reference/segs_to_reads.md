# split segments into bins of a requested size.

Takes a dataframe with segment start and end columns and returns a
dataframe with those segments split into individual bins. All segment
information applied to the generated bins.

## Usage

``` r
segs_to_reads(
  segs_df,
  bin_size = 5e+05,
  seg_start_col = "start",
  seg_end_col = "end"
)
```

## Arguments

- segs_df:

  dataframe (or similar) of segment dat

- bin_size:

  width of bins to split into. Default is standard 500kb

- seg_start_col:

  name of the column that indicates the start of a segment

- seg_end_col:

  name of the column that indicates the end of a segment

## Value

input frame split into bins

## Details

Warnings issues if requested bin size is bigger than some segments or if
segments can't be split evenly into the bins.

Bin end will not exceed a segment end. Depending on what is input and
requested, this can lead to one smaller bin at the end of the segment.
