# thin wrapper around sigminer feature tally

Reshapes/renames our typical dataframes into one sigminer is happy with
and runs sigminer::sig_tally with the "W" method.

## Usage

``` r
extract_sigminer_wang_features(segs_df)
```

## Arguments

- segs_df:

  dataframe of CN segments.

## Value

matrix of feature counts.

## Details

Features are based on the paper:

Wang et al. Copy number signature analysis tool and its application in
prostate cancer reveals distinct mutational processes and clinical
outcomes. PLOS Genetics. 2021.

https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009557

Some issues, I, Ben Furman, have with it are:

- any oscillations of CN count as an oscillation, including something
  like 2 - 500 - 2. Oscillations are supposed to reflect potential
  chromothripsis, which a pattern like that is not likely to be from.

- segment sizes are binned based on log10 values. This means there are a
  lot of CN size bins \< 1 Mb, and DLP starts at 500Kb (and segments
  that small should probably be filtered anyway). Thus, many of the
  segment size bins are not used with DLP. This isn't a sigminer issue,
  but an issue with the scale of DLP data.
