# extract sample ID from the typically formatted cell_ids

expecting cell IDs as AT21350-A143952A-R10-C37 with the first position
being the sample ID.

## Usage

``` r
sample_from_cell(cell_id)
```

## Arguments

- cell_id:

  string of a cell_id or vector of cell IDs

## Value

vector of the sample ID(s) contained within
