# internal for checking the names of chromosome columns in a dataframe.

used for error catching.

## Usage

``` r
chr_name_check(df, exp_chr_name)
```

## Arguments

- df:

  a dataframe being manipulated in some function

- exp_chr_name:

  a string of the chromosome column name the function is expecting to
  see.

## Value

bool TRUE if all good, FALSE with a message otherwise.
