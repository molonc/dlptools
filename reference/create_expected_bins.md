# create bins of given size for a genome

create bins of given size for a genome

## Usage

``` r
create_expected_bins(version = c("hg19", "hg38"), bin_size = 5e+05)
```

## Arguments

- version:

  default "hg19", or can select "hg38"

- bins_size:

  int of how big to make the bins

## Value

tibble of bins for each chromosome. Columns: chr start end
