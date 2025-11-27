# Mark if CN states are gains or losses relative to cell ploidy

CNs \> 2 are not necessarily amplifications, and CNs \< 2 are not
necessarily the only states of losses. Ploidy of samples may not be
diploid. This function will give you an idea if the CN you see is a gain
or loss relative to the ploidy of a sample. For example, if the sample
has a ploidy of 4, then a CN of 3 is a loss.

## Usage

``` r
mark_cn_relative_to_ploidy(
  in_df,
  df_type = c("reads", "segs"),
  sample_col = "cell_id",
  ...
)
```

## Arguments

- in_df:

  dataframe of CN states.

- df_type:

  string. "reads" (default) or "segs" for CN segments, which will
  internally converted to bin based for mode calculation. reads in order
  to infer mode ploidy.

- sample_col:

  string. Name of the column with cell_id/other sample name

## Value

input dataframe, with new columns of information.

## Details

Ploidy is inferred by mode CN state using
[`mode_ploidy()`](https://molonc.github.io/dlptools/reference/mode_ploidy.md)
and states \> ploidy are marked as gains, states \< ploidy losses, and
states matching ploidy as matched.
