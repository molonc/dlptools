# clean tree tip labels and drop any locus tips from sitka trees

'Locus tips' are from sitka and are locus values that end up on the tip
of trees. Also removes the "cell\_" prefix from tip labels, which is
also a consequence of sitka.

## Usage

``` r
format_sitka_tree(tree)
```

## Arguments

- tree:

  phylo object as read by ape::read.tree

## Value

phylo object cleaned of "cell\_" notation
