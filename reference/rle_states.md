# convert states to run length encoding

takes a vector of numbers (e.g., states) and returns a numeric group
value indicating the a group for each run of the same value.
c(5,5,5,6,6,5,5,5,2) -\> 1 1 1 2 2 3 3 3 4

## Usage

``` r
rle_states(states)
```

## Arguments

- states:

  really and vector of values.

## Value

vector of integers
