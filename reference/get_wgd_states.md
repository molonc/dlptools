# Calculate WGD states based on major allele CN

Simple metric. If \>50% of the genome has 2 or more copies of the
"major" allele ("A"; referred to as "FM2" in the paper), it's at least
1x WGD. If there are 3 ("FM3") or more copies of A for \>50% of the
genome, then it's at least 2x WGD. Otherwise, there is no support for
WGD.

## Usage

``` r
get_wgd_states(
  signals_df,
  cell_col = "cell_id",
  start_col = "start",
  end_col = "end",
  chrom_col = "chr",
  equalize_chromosomes = FALSE
)
```

## Arguments

- signals_df:

  dataframe/tibble from signals, with A and B allele columns

- cell_col:

  string column name identifying cell IDs

- start_col:

  string column name identifying start of segments

- end_col:

  string column name identifying end of segments

- chrom_col:

  string column name for chromosomes

- equalize_chromosomes:

  boolean. See description.

## Value

tibble of per-sample ploidy estimates

## Details

Requires a [signals](https://shahcompbio.github.io/signals/) input, and
is based on the paper by [McPherson et al.
2025](https://www.nature.com/articles/s41586-025-09240-3).

Functions returns dataframe with "wgd_state" column, along with the
fraction estimates. Also calls on
[`weighted_ploidy()`](https://molonc.github.io/dlptools/reference/weighted_ploidy.md)
and adds those results.

A potential issue with this approach is that larger chromosomes
disproportionately contribute to the result (e.g., sum of chr 1 and 2 is
more than many of the smaller chromosomes combined). Setting
`equalize_chromosomes` to TRUE will perform the calculations within each
chromosome, then assess if \>50% of the chromosomes show FM2 \>50%, or
FM3 \>50% . This way, all chromosomes contribute equally to the
inference of WGD.
