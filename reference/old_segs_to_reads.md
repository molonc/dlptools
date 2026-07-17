# original function to split segments back to read bins

kept for posterity. This splits into exact 500kb bins, which may extend
beyond end of chromosome, as original DLP did. New function is wildly
faster and more simple.

## Usage

``` r
old_segs_to_reads(
  segs_df,
  bin_size = 5e+05,
  seg_start_col = "start",
  seg_end_col = "end"
)
```
