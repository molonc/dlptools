# For a given length, create bins of a given size

e.g., length = 10, bin = 5, will get bins: 1-5, 6-10

## Usage

``` r
expand_length_to_bins(total_len, bin_size = 5e+05)
```

## Arguments

- total_len:

  int to create bins for

- bin_size:

  int of size of bin

## Value

tibble of bins with columns of start and end
