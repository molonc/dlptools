# count values that map to given categories

count values that map to given categories

## Usage

``` r
cut_categories_and_count(df, targ_col, sample_col, breaks, labels)
```

## Arguments

- df:

  dataframe. Contains a column of values you want to count per sample

- targ_col:

  string. Column you want to put in categories and counted.

- sample_col:

  string. Name of the column with sample ids.

- breaks:

  vector of doubles. Bounds for the categories.

- labels:

  vector of strings. What to call the categories.

## Value

tibble
