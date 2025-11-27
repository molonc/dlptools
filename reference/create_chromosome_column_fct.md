# create a sorted factor vector of chromosomes

This is mostly just for splitting the the states heatmap by chromosome.
Will be naturally sorted.

## Usage

``` r
create_chromosome_column_fct(states_mat)
```

## Arguments

- states_mat:

  matrix of cell_id named rows and bin (chr_start_end) columns

## Value

factor of chromosomes with sorted levels
