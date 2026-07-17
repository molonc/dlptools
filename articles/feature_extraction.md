# Feature Extraction for CN Signatures

``` r

library(ggplot2)
library(dlptools)
```

## Feature Extraction

`dlptools` has a series of functions to extract a variety of features.
These could be used for subsequent signature extraction, or be use to
quantify differences between samples.

On offer are the “classic” process based features, inspired by

- [Macintyre et
  al. 2018](https://www.nature.com/articles/s41588-018-0179-8)
- [Maclachlan et
  al. 2021](https://www.nature.com/articles/s41467-021-25469-8)
- [Drews et
  al. 2022](https://www.nature.com/articles/s41586-022-04789-9)
- and others

These include:

1.  segment size
2.  breakpoints per arm
3.  breakpoints per window
4.  changepoints: change in CN state between adjacent segments
5.  length of oscillating chains
6.  CN state (though there is debate if this should be used at all)

  

The functions for these (outlined below) will extract the observed
values for each of these features, and offers to summarise counts in
pre-defined categories. These categories are informed by the
publications above, but designed to be accommodating for DLP data, which
has a more coarse resolution than most of those publications.

  

Additionally, there are functions to extract two externally defined
feature extraction types:

1.  Wu et
    al. <https://www.biorxiv.org/content/10.1101/2025.03.02.641098v1>
2.  Wang et
    al. <https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009557>

The Wang et al features are just a thin wrapper around
[sigminer](https://github.com/ShixiangWang/sigminer?tab=readme-ov-file)
functions. Check out that package for more tools and actual signature
extraction.

  
  

A likely starting point of this exercise might be some read bin data, or
some segment data. It’s probably good that it’s been cleaned up to your
liking.

Perhaps the workflow would look something like this:

``` r

reads_df <- vroom::vroom("data/example_full_reads.tsv.gz")

segs_df <- reads_df |>
  dplyr::filter(!chr %in% c("X", "Y")) |>
  dlptools::mark_mask_regions() |>
  dlptools::mark_bins_overlapping_centromeres(padding = 3e6) |>
  dplyr::filter(
    gc != -1 & map > 0.99 & !overlaps_centro & !mask
  ) |>
  dlptools::add_missing_bins_for_cells() |>
  dlptools::fill_state_from_neighbours() |>
  dlptools::reads_to_segs() |>
  dplyr::filter(
    seg_width > 1.5e6
  ) |>
  dlptools::segs_to_reads() |>
  dlptools::add_missing_bins_for_cells() |>
  dlptools::fill_state_from_neighbours() |>
  dlptools::reads_to_segs()


dplyr::slice_head(segs_df, n = 5)
#> # A tibble: 5 × 6
#>   cell_id                  chr      start      end state seg_width
#>   <chr>                    <chr>    <dbl>    <dbl> <dbl>     <dbl>
#> 1 AT21352-A144173A-R03-C12 1            1 10000000     4   9999999
#> 2 AT21352-A144173A-R03-C12 1     10000001 21500000     5  11499999
#> 3 AT21352-A144173A-R03-C12 1     21500001 39000000     4  17499999
#> 4 AT21352-A144173A-R03-C12 1     39000001 45000000     5   5999999
#> 5 AT21352-A144173A-R03-C12 1     45000001 57500000     4  12499999
```

  
  
  
  

## Process-Based Features

For the 6 functions in this section, each by default return the observed
values. But, each function also has an option to `return = counts`,
which will use pre-defined categories.

  

1.  segment lengths.

This is as simple as it sounds, and the function is unnecessary. These
lengths would need to be subsequently categorized somehow.

``` r

dlptools::extract_segment_sizes(
  segs_df = segs_df
) |>
  dplyr::slice_head(n = 10)
#> # A tibble: 10 × 2
#>    cell_id                  seg_size
#>    <chr>                       <dbl>
#>  1 AT21352-A144173A-R03-C12 10000000
#>  2 AT21352-A144173A-R03-C12 11500000
#>  3 AT21352-A144173A-R03-C12 17500000
#>  4 AT21352-A144173A-R03-C12  6000000
#>  5 AT21352-A144173A-R03-C12 12500000
#>  6 AT21352-A144173A-R03-C12 18500000
#>  7 AT21352-A144173A-R03-C12  3000000
#>  8 AT21352-A144173A-R03-C12  7500000
#>  9 AT21352-A144173A-R03-C12  9500000
#> 10 AT21352-A144173A-R03-C12  8500000
```

  

alternatively you can specify the `return = counts` to get counts in
pre-defined categories. This is an argument available for all of these
process based functions.

``` r

dlptools::extract_segment_sizes(
  segs_df = segs_df,
  return = "counts"
) |>
  dplyr::slice_head(n = 10)
#> # A tibble: 120 × 3
#> # Groups:   cell_id, feat_cat [120]
#>    cell_id                  feat_cat        n
#>    <chr>                    <fct>       <int>
#>  1 AT21352-A144173A-R03-C12 ss:<5Mb        21
#>  2 AT21352-A144173A-R03-C12 ss:5-10Mb      41
#>  3 AT21352-A144173A-R03-C12 ss:10-20Mb     39
#>  4 AT21352-A144173A-R03-C12 ss:20-50Mb     26
#>  5 AT21352-A144173A-R03-C12 ss:50-100Mb    12
#>  6 AT21352-A144173A-R03-C12 ss:100+Mb       3
#>  7 AT21352-A144173A-R03-C32 ss:<5Mb        31
#>  8 AT21352-A144173A-R03-C32 ss:5-10Mb      31
#>  9 AT21352-A144173A-R03-C32 ss:10-20Mb     21
#> 10 AT21352-A144173A-R03-C32 ss:20-50Mb     23
#> # ℹ 110 more rows
```

categories are: \< 5mb, 5-10 Mb, 20-50 Mb, 50-100 Mb, and 100+Mb.

  

2.  Copy number change points

These are the values of the change in copy number from segment to
segment. E.g., 2-5 would yield a 3.

See
[`?dlptools::extract_changepoint`](https://molonc.github.io/dlptools/reference/extract_changepoint.md)
for an important discussion of the `first_seg_correction` parameter.

``` r

dlptools::extract_changepoint(
  segs_df = segs_df
  # or returning the counts
  # return = 'counts'
) |>
  dplyr::slice_head(n = 10)
#> # A tibble: 10 × 2
#>    cell_id                  cn_change
#>    <chr>                        <dbl>
#>  1 AT21352-A144173A-R03-C12         1
#>  2 AT21352-A144173A-R03-C12         1
#>  3 AT21352-A144173A-R03-C12         1
#>  4 AT21352-A144173A-R03-C12         1
#>  5 AT21352-A144173A-R03-C12         1
#>  6 AT21352-A144173A-R03-C12         2
#>  7 AT21352-A144173A-R03-C12         2
#>  8 AT21352-A144173A-R03-C12         1
#>  9 AT21352-A144173A-R03-C12         1
#> 10 AT21352-A144173A-R03-C12         2
```

Predefined count categories are changepoint = 1, 2, 3, 4, 5+

  

3.  Extracting number of breakpoints per chromosome arm

A breakpoint is a change in the copy number.

``` r

dlptools::extract_breakpoints(
  segs_df = segs_df,
  scope = "arms"
  # or returning the counts
  # return = 'counts'
) |>
  dplyr::slice_head(n = 10)
#> # A tibble: 10 × 4
#>    cell_id                  chr   arm   n_bps
#>    <chr>                    <chr> <chr> <int>
#>  1 AT21352-A144173A-R03-C12 1     p        10
#>  2 AT21352-A144173A-R03-C12 1     q         5
#>  3 AT21352-A144173A-R03-C12 10    p         2
#>  4 AT21352-A144173A-R03-C12 10    q         1
#>  5 AT21352-A144173A-R03-C12 11    p         0
#>  6 AT21352-A144173A-R03-C12 11    q         4
#>  7 AT21352-A144173A-R03-C12 12    p         3
#>  8 AT21352-A144173A-R03-C12 12    q         2
#>  9 AT21352-A144173A-R03-C12 13    p         0
#> 10 AT21352-A144173A-R03-C12 13    q         1

# there is an alias of: dlptools::extract_bp_per_arm()
```

Predefined count categories are: breakpoint = 0, 1, 2, 3, 4, 5, 6+

  

4.  Extract breakpoints per MB window

Typically, publications have used 10Mb.

``` r

dlptools::extract_breakpoints(
  segs_df = segs_df,
  scope = "windows"
  # window_size = 10e6 # default
  # or returning the counts
  # return = 'counts'
) |>
  dplyr::slice_head(n = 10)
#> # A tibble: 10 × 6
#>    cell_id                  chr   window_name n_bps     start       end
#>    <chr>                    <chr> <chr>       <int>     <int>     <int>
#>  1 AT21352-A144173A-R03-C12 1     w1              0         1  10000000
#>  2 AT21352-A144173A-R03-C12 1     w10             1  90000001 100000000
#>  3 AT21352-A144173A-R03-C12 1     w11             1 100000001 110000000
#>  4 AT21352-A144173A-R03-C12 1     w12             0 110000001 120000000
#>  5 AT21352-A144173A-R03-C12 1     w13             0 120000001 130000000
#>  6 AT21352-A144173A-R03-C12 1     w14             0 130000001 140000000
#>  7 AT21352-A144173A-R03-C12 1     w15             0 140000001 150000000
#>  8 AT21352-A144173A-R03-C12 1     w16             0 150000001 160000000
#>  9 AT21352-A144173A-R03-C12 1     w17             0 160000001 170000000
#> 10 AT21352-A144173A-R03-C12 1     w18             1 170000001 180000000

# there is an alias of: dlptools::extract_bp_per_window()
```

Predefined count categories are: breakpoint = 0, 1, 2, 3+

  

5.  lengths of copy number oscillations

Overall, this is a bit of a funny feature. But basically, chromosomes
are assessed in 3-segment sets to determine if it counts as an
oscillation. If there are more than one oscillations in a row, they are
added up. If the 3-segment set is not an oscillation, it gets a 0.

  

- 1-2-3 is not an oscillation, length 0
- 2-3-2 is an oscillation of length 1
- 2-3-2-3 is an oscillation of length 2
- 2-3-2-3-5 is an oscillation of length 2 followed by a 0 for the
  non-oscillating component
- 1-2-3-4-5 has no oscillations, and returns 0 0 0 for the 3 possible
  3-set segments

  

I have dissected the code for various publications and packages that
count this feature and they all work differently, some have outright
bugs, and others lack some internal consistency in how they assess 0s.
So, please read
[`?dlptools::extract_oscillations`](https://molonc.github.io/dlptools/reference/extract_oscillations.md)
for a long discussion of how exactly this is working.

  

Overall, 0s are typically abundant, and should maybe be considered for
removal…? I don’t know, they just tend to overwhelm. Again, all
publications using this feature have done so differently, so it’s the
wild west. Do what feels right.

``` r

dlptools::extract_oscillations(
  segs_df = segs_df
  # or returning the counts
  # return = 'counts'

  # other important arguments are:
  # middle_bound = 2 # default
  # ends_bound = 0 # default
  # see function documentation for what these are
) |>
  dplyr::slice_head(n = 10)
#> # A tibble: 10 × 2
#>    cell_id                    OSc
#>    <chr>                    <dbl>
#>  1 AT21352-A144173A-R03-C12     4
#>  2 AT21352-A144173A-R03-C12     0
#>  3 AT21352-A144173A-R03-C12     1
#>  4 AT21352-A144173A-R03-C12     0
#>  5 AT21352-A144173A-R03-C12     1
#>  6 AT21352-A144173A-R03-C12     0
#>  7 AT21352-A144173A-R03-C12     0
#>  8 AT21352-A144173A-R03-C12     0
#>  9 AT21352-A144173A-R03-C12     0
#> 10 AT21352-A144173A-R03-C12     2
```

Predefined count categories are: Oscillation chains of length: 0, 1-3,
4-9, 10+.

E.g.,

``` r

dlptools::extract_oscillations(
  segs_df = segs_df,
  return = "counts"
) |>
  dplyr::slice_head(n = 10)
#> # A tibble: 80 × 3
#> # Groups:   cell_id, feat_cat [80]
#>    cell_id                  feat_cat     n
#>    <chr>                    <fct>    <int>
#>  1 AT21352-A144173A-R03-C12 osc:0       70
#>  2 AT21352-A144173A-R03-C12 osc:1-3     22
#>  3 AT21352-A144173A-R03-C12 osc:4-9      1
#>  4 AT21352-A144173A-R03-C12 osc:10+      0
#>  5 AT21352-A144173A-R03-C32 osc:0       60
#>  6 AT21352-A144173A-R03-C32 osc:1-3     18
#>  7 AT21352-A144173A-R03-C32 osc:4-9      0
#>  8 AT21352-A144173A-R03-C32 osc:10+      0
#>  9 AT21352-A144173A-R03-C33 osc:0       75
#> 10 AT21352-A144173A-R03-C33 osc:1-3     13
#> # ℹ 70 more rows
```

  

You can also play directly with
[`dlptools::count_oscillations`](https://molonc.github.io/dlptools/reference/count_oscillations.md)
on vectors to see how the inernals of this function are working:

``` r

dlptools::count_oscillations(1:3)
#> [1] 0

dlptools::count_oscillations(c(2, 3, 2))
#> [1] 1

dlptools::count_oscillations(c(2, 3, 2, 3))
#> [1] 2

dlptools::count_oscillations(c(2, 3, 2, 3, 8, 5, 4, 5))
#> [1] 2 0 0 0 1

# most other publications/package do not care about the middle value
dlptools::count_oscillations(
  c(5, 1000, 5),
  middle_bound = Inf
)
#> [1] 1

# some publications would allow this as counting as an oscillation. Where the
# first and last segments are allowed to be within a certain distance of one
# another, and the middle value doesn't matter.
dlptools::count_oscillations(
  c(5, 1000, 6),
  middle_bound = Inf,
  ends_bound = 1
)
#> [1] 1
```

  

6.  CN States

And for sake of completeness, also offered is CN. By default, function
just pulls the CN state column from your dataframe.

Really only slightly useful if you want to use the pre-defined
categories of:

CN = 1, 2, 3, 4, 5, 6+

``` r

dlptools::extract_cn(
  segs_df,
  return = "counts"
) |>
  dplyr::slice_head(n = 10)
#> # A tibble: 160 × 3
#> # Groups:   cell_id, feat_cat [160]
#>    cell_id                  feat_cat     n
#>    <chr>                    <fct>    <int>
#>  1 AT21352-A144173A-R03-C12 CN:0         0
#>  2 AT21352-A144173A-R03-C12 CN:1         4
#>  3 AT21352-A144173A-R03-C12 CN:2        19
#>  4 AT21352-A144173A-R03-C12 CN:3        35
#>  5 AT21352-A144173A-R03-C12 CN:4        38
#>  6 AT21352-A144173A-R03-C12 CN:5        25
#>  7 AT21352-A144173A-R03-C12 CN:6         4
#>  8 AT21352-A144173A-R03-C12 CN:7+       17
#>  9 AT21352-A144173A-R03-C32 CN:0         0
#> 10 AT21352-A144173A-R03-C32 CN:1         3
#> # ℹ 150 more rows
```

  
  

If you are happy with all of the default categories for counting, you
can call the wrapper to get them all at once:

``` r

# if you set up a future::plan(), this will run in parallel. See ?dlptools::count_default_process_feats

# could be something like
# future::plan(future::multicore, workers = 10)

dlptools::count_default_process_feats(
  segs_df
  # CN not included by default, use this to include it
  # with_cn=TRUE
) |>
  dplyr::slice_head(n = 10)
#> # A tibble: 10 × 3
#>    cell_id                  feat_cat        n
#>    <chr>                    <fct>       <int>
#>  1 AT21352-A144173A-R03-C12 ss:<5Mb        21
#>  2 AT21352-A144173A-R03-C12 ss:5-10Mb      41
#>  3 AT21352-A144173A-R03-C12 ss:10-20Mb     39
#>  4 AT21352-A144173A-R03-C12 ss:20-50Mb     26
#>  5 AT21352-A144173A-R03-C12 ss:50-100Mb    12
#>  6 AT21352-A144173A-R03-C12 ss:100+Mb       3
#>  7 AT21352-A144173A-R03-C12 CP:0            2
#>  8 AT21352-A144173A-R03-C12 CP:1           72
#>  9 AT21352-A144173A-R03-C12 CP:2           33
#> 10 AT21352-A144173A-R03-C12 CP:3            9
```

  
  

### Wu et al.

For extracting Wu et al. features:

``` r

dlptools::extract_wu_features(
  segs_df = segs_df
) |>
  dplyr::slice_head(n = 5)
#> # A tibble: 5 × 7
#>   cell_id                  cn_bin seg_bin seg_change seg_shape     n feat_cat   
#>   <chr>                    <fct>  <fct>   <fct>      <fct>     <int> <chr>      
#> 1 AT21352-A144173A-R03-C12 S0-1   S       AA         LL            0 S:LL:S0-1:…
#> 2 AT21352-A144173A-R03-C12 S0-1   S       AA         HH            0 S:HH:S0-1:…
#> 3 AT21352-A144173A-R03-C12 S0-1   S       AA         OT            0 S:OT:S0-1:…
#> 4 AT21352-A144173A-R03-C12 S0-1   S       BB         LL            1 S:LL:S0-1:…
#> 5 AT21352-A144173A-R03-C12 S0-1   S       BB         HH            0 S:HH:S0-1:…
```

  
  

There are other options, like annotating the input dataframe, or
returning a matrix for subsequent signature extraction:

``` r

feat_mtx <- dlptools::extract_wu_features(
  segs_df = segs_df,
  return_matrix = TRUE
)

feat_mtx[1:4, 1:4]
#>                          S:LL:S0-1:AA S:HH:S0-1:AA S:OT:S0-1:AA S:LL:S0-1:BB
#> AT21352-A144173A-R03-C12            0            0            0            1
#> AT21352-A144173A-R03-C32            0            0            0            1
#> AT21352-A144173A-R03-C33            0            0            0            0
#> AT21352-A144173A-R03-C35            0            0            0            0

# alternatively, there is a generic dlptools::make_cellid_matrix() than can
# pivot out and make the matrix if you already have the count dataframe and
# don't want to redo the feature extraction
```

  

see
[`?dlptools::extract_wu_features`](https://molonc.github.io/dlptools/reference/extract_wu_features.md)
for more.

  
  
  

### Wang et al., aka Sigminer

Alternatively, we can get Wang et al. features. Though if going this
route, it would be a good idea to look more at
[sigminer](https://github.com/ShixiangWang/sigminer?tab=readme-ov-file)
itself:

``` r

sig_feats <- extract_sigminer_wang_features(segs_df)

sig_feats[1:4, 1:4]
#>                          BP10MB[0] BP10MB[1] BP10MB[2] BP10MB[3]
#> AT21352-A144173A-R03-C12       198        87        12         3
#> AT21352-A144173A-R03-C32       215        69        13         3
#> AT21352-A144173A-R03-C33       201        85        13         0
#> AT21352-A144173A-R03-C35       224        65         8         3
```

  
  
  

## Other vaguely feature like things

These could be considered features, maybe only within certain contexts.
But, maybe they are just more useful as descriptions of the data as a
whole.

### ploidy informed CN change

CN change relative to ploidy might be interesting to have a feature.
Here we can mark 3 types of events for segments:

1.  CN is a gain relative to mode ploidy of sample/cell
2.  CN is a loss relative to mode ploidy of sample/cell
3.  CN is a match to mode ploidy

Employs the function dlptools::mode_ploidy() to calculate mode ploidy of
each sample/cell.

``` r

dlptools::extract_ploidy_cn_feature(
  segs_df,
  # with option just to annotate the input df
  # annotate_input = TRUE
  # or, alternatively, return a count matrix
  # return_matrix = TRUE
) |>
  dplyr::slice_head(n = 6)
#> # A tibble: 6 × 3
#>   cell_id                  cn_v_ploidy      n
#>   <chr>                    <fct>        <int>
#> 1 AT21352-A144173A-R03-C12 ploidy-gain     84
#> 2 AT21352-A144173A-R03-C12 ploidy-match    35
#> 3 AT21352-A144173A-R03-C12 ploidy-loss     23
#> 4 AT21352-A144173A-R03-C32 ploidy-gain     84
#> 5 AT21352-A144173A-R03-C32 ploidy-match    26
#> 6 AT21352-A144173A-R03-C32 ploidy-loss     16
```

  

### Segment Span Types

It might be useful to know where segments sit on chromosomes and how
much of a chromosome they cover. We can label segments based on 5
categories:

1.  telomere bound (`telo-bound`) - segment touches a telomere
2.  centromere bound (`centro-bound`) - segment touches or crosses the
    centromere
3.  arm (`arm`) - segment spans a whole are (\*with conditions)
4.  whole chromosome (`whole-chrom`) - segment spans the entire
    chromosome (\*with conditions)
5.  intersitial (`inter`) - occuring within the chromosome, not touching
    the centromere, telomeres, and not big enough to be an entire arm.

See
[`?dlptools::mark_segs_chromosome_span`](https://molonc.github.io/dlptools/reference/mark_segs_chromosome_span.md)
for more details about this marking, as there are various options to
tweak these definitions, such as:

- how much of an arm a CN needs to cover to be considered an “arm-level”
  event
- how much of a chromosome a segment needs to cover to be a
  “whole-chromosome” event
- how close a seg needs to be to a telo/centromere to be considered
  “touching”

We can mark a dataframe with the information:

``` r

dlptools::mark_segs_chromosome_span(
  segs_df,
  # see ?dlptools::mark_segs_chromosome_span for more args
) |>
  dplyr::relocate(
    # key added columns
    seg_span_event, telo_bound, centro_bound, spans_chrom, spans_arm
  ) |>
  dplyr::slice_head(n = 5)
#> # A tibble: 5 × 32
#>   seg_span_event telo_bound centro_bound spans_chrom spans_arm cell_id     chr  
#>   <fct>          <lgl>      <lgl>        <lgl>       <lgl>     <chr>       <chr>
#> 1 telo-bound     TRUE       FALSE        FALSE       FALSE     AT21352-A1… 1    
#> 2 inter          FALSE      FALSE        FALSE       FALSE     AT21352-A1… 1    
#> 3 inter          FALSE      FALSE        FALSE       FALSE     AT21352-A1… 1    
#> 4 inter          FALSE      FALSE        FALSE       FALSE     AT21352-A1… 1    
#> 5 inter          FALSE      FALSE        FALSE       FALSE     AT21352-A1… 1    
#> # ℹ 25 more variables: start <dbl>, end <dbl>, state <dbl>, seg_width <dbl>,
#> #   total_length <dbl>, centro_start <dbl>, centro_end <dbl>,
#> #   centro_span <dbl>, start_p <dbl>, start_q <dbl>, end_p <dbl>, end_q <dbl>,
#> #   telostart_p <dbl>, teloend_p <dbl>, telostart_q <dbl>, teloend_q <dbl>,
#> #   telo_p_dist <dbl>, telo_q_dist <dbl>, telo_dist <dbl>, centro_p_dist <dbl>,
#> #   centro_q_dist <dbl>, centro_dist <dbl>, spans_centro <lgl>,
#> #   p_arm_length <dbl>, q_arm_length <dbl>
```

  
  

More useful in the “feature” context, is to directly extract and count
them up. This function bascially just calls the one above, with counting
added:

``` r

dlptools::extract_segment_position_feature(
  segs_df
  # with option just to annotate the input df
  # annotate_input = TRUE
  # or, alternatively, return a count matrix
  # return_matrix = TRUE
) |>
  dplyr::slice_head(n = 6)
#> # A tibble: 6 × 3
#>   cell_id                  seg_span_event     n
#>   <chr>                    <fct>          <int>
#> 1 AT21352-A144173A-R03-C12 arm                1
#> 2 AT21352-A144173A-R03-C12 whole-chrom        2
#> 3 AT21352-A144173A-R03-C12 telo-bound        34
#> 4 AT21352-A144173A-R03-C12 centro-bound      20
#> 5 AT21352-A144173A-R03-C12 inter             85
#> 6 AT21352-A144173A-R03-C32 arm                0
```

  
  
  
