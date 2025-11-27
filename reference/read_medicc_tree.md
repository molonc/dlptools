# read medicc tree to phylo object.

Medicc forces the inclusion of a diploid ancestor for all cells and
calls that tip "diploid". This optionally drops that tip.

## Usage

``` r
read_medicc_tree(tree_file, drop_diploid = TRUE)
```

## Arguments

- tree_file:

  file path to newick tree created by medicc

- drop_dipliod:

  boolean. Default TRUE. Remove diploid tip from tree.

## Value

phylo object
