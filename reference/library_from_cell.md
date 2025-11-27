# extract library ID from the typically formatted cell_ids

expecting cell IDs as AT21350-A143952A-R10-C37 with the second position
being the library ID.

## Usage

``` r
library_from_cell(cell_id)
```

## Arguments

- cell_id:

  string of a cell_id or vector of cell IDs

## Value

vector of the library ID(s) contained within
