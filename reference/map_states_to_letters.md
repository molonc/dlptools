# convert a vector of states to letters

The point is to then use measure string distance between cells. With raw
states double digit values would count as 2 characters and throw things
off, with standard R-functions to measure string differences.

## Usage

``` r
map_states_to_letters(states)
```

## Arguments

- states:

  vector of state values

## Value

vector of letters corresponding to those states
