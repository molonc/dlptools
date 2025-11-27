# wrapper around classic process based features

Extracts counts using the default set categories of the process based
features, including:

## Usage

``` r
count_default_process_feats(segs_df, sample_col = "cell_id", with_cn = FALSE)
```

## Arguments

- segs_df:

  dataframe. Segmented CN for samples.

- sample_col:

  string. Column name of sample ID column.

- with_cn:

  bool. default FALSE. FALSE means do not include CN in the feature
  counts, TRUE will include it.

## Details

1.  [extract_segment_sizes](https://molonc.github.io/dlptools/reference/extract_segment_sizes.md)

2.  [extract_changepoint](https://molonc.github.io/dlptools/reference/extract_changepoint.md)

3.  [extract_bp_per_arm](https://molonc.github.io/dlptools/reference/extract_bp_per_arm.md)

4.  [extract_bp_per_window](https://molonc.github.io/dlptools/reference/extract_bp_per_window.md)

5.  [extract_oscillations](https://molonc.github.io/dlptools/reference/extract_oscillations.md)

Not included is CN, which can be added by specifying with_cn=TRUE

Internally using furrr::future_map, so if you set up a future::plan(),
this will execute in parallel. Recommended. See Example.

## Examples

``` r
if (FALSE) { # \dontrun{
future::plan(future::multicore, workers = 10)
dlptools::count_default_process_feats(segs_df, sample_id = "cell_id")
} # }
```
