# get tips labels that will avoid duplicate sibling comparisons

I.e., in a tree: (A, (B, C)) we don't want to compare B to C and then C
to B, we only need to do one of those comparisons.

## Usage

``` r
get_tips_that_avoid_redundant_comps(tree)
```

## Arguments

- tree:

  a phylo object

## Value

vector of tips labels that will lead to non-redundant
