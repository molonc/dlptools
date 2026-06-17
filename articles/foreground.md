# Copy Number 'Foreground'

``` r

library(dlptools)
```

“Foreground” *in the world of DLP+, in our lab*, and **at the current
moment**, is the concept of labelling copy number alterations that are
specific to a given cell.

These are events that are happening in single cells during division
generating variation among cells in a sequenced sample. The idea is that
this is the variation that selection ultimately functions on to fix
changes in tumours, driving tumour evolution.

How exactly to go about doing this is still a work in progress. The code
available here is under development and will change. Included are some
functions and plotting capabilities to make current efforts easier to
apply across DLP libraries.

## General Helpers

Retrieve node labels of immediate parents to tree tips:

``` r

ex_tree <- ape::read.tree("data/pkg_tree.newick")

tip_parents <- dlptools::get_tip_parents(ex_tree)

dplyr::slice_head(tip_parents, n = 10)
#> # A tibble: 10 × 2
#>    cell_id                  parent_node
#>    <chr>                    <chr>      
#>  1 AT23998-A138956A-R09-C10 internal583
#>  2 AT28335-A143820B-R64-C14 internal494
#>  3 AT28335-A143820B-R55-C55 internal465
#>  4 AT28335-A143820B-R57-C22 internal465
#>  5 AT23998-A138956A-R08-C27 internal499
#>  6 AT28335-A143820B-R65-C13 internal469
#>  7 AT28335-A143820B-R66-C49 internal469
#>  8 AT28335-A143820B-R53-C16 internal492
#>  9 AT23998-A138956A-R03-C61 internal507
#> 10 AT23998-A138956A-R19-C49 internal502
```

the result is a dataframe with the tip and it’s immediate parent node in
a tree.

  
  

Typically, we’d also have dataframe of states for bins in cells, and
this simple wrapper exists that will add this parent node information to
the dataframe:

``` r

states <- vroom::vroom("data/ex_state_dat.tsv.gz", show_col_types = FALSE)

states_w_node_info <- dlptools::add_tip_ancestors_to_df(states, ex_tree)

states_w_node_info |>
  dplyr::select(cell_id, chr, start, end, state, parent_node) |>
  dplyr::slice_head(n = 10)
#> # A tibble: 10 × 6
#>    cell_id                    chr   start     end state parent_node
#>    <chr>                    <dbl>   <dbl>   <dbl> <dbl> <chr>      
#>  1 AT23998-A138956A-R03-C59     1 2000001 2500000     4 internal582
#>  2 AT23998-A138956A-R03-C59     1 3000001 3500000     4 internal582
#>  3 AT23998-A138956A-R03-C59     1 4000001 4500000     4 internal582
#>  4 AT23998-A138956A-R03-C59     1 4500001 5000000     5 internal582
#>  5 AT23998-A138956A-R03-C59     1 5000001 5500000     5 internal582
#>  6 AT23998-A138956A-R03-C59     1 5500001 6000000     5 internal582
#>  7 AT23998-A138956A-R03-C59     1 6000001 6500000     5 internal582
#>  8 AT23998-A138956A-R03-C59     1 6500001 7000000     5 internal582
#>  9 AT23998-A138956A-R03-C59     1 7000001 7500000     5 internal582
#> 10 AT23998-A138956A-R03-C59     1 7500001 8000000     5 internal582
```

## Medicc

At the moment, one method being explored for “foreground” labelling is
with Medicc. Medicc, which we use to build phylogenetic trees, as part
of it’s methods builds CN profiles of the internal nodes of a tree that
reflect the input data.

So, if you input a set of bins with CN estimates for cells, it will
output the CN calls for that set of bins for each internal tree node.
Assessing the changes in those CN states between tips and parents is one
method to look at foreground changes.

There are a ton of caveats to this and limits to the model, best
discussed elsewhere.

Here, we’ll using some example data, which is a reduced set of data from
a full library. First loading the data:

``` r

# we can ingest the medicc tree. This is a simple wrapper to drop the diploid
# branch, which medicc invents and includes
med_tree <- dlptools::read_medicc_tree("data/ex_med_fg_tree.nwk")

# and we can import the profiles file output by medicc
# this function does some simple changes to the file to make things easier
med_profiles <- dlptools::read_medicc_profiles("data/ex_profiles.tsv.gz")
```

  

and now we can infer the changes between tips and immediate parents:

``` r

states_df <- dlptools::medicc_profiles_to_foreground(
  med_profiles, med_tree,
  cn_type = "total"
)
```

The result is a dataframe of the state data given to medicc for the
tips, and added columns of the states of the immediate parent node of
those tips.

``` r

states_df |>
  dplyr::select(
    cell_id, chr, start, end, state, parent_state, parent_node,
    fg_change, fg_type, foreground, background
  ) |>
  dplyr::slice_sample(n = 10)
#> # A tibble: 10 × 11
#>    cell_id  chr    start    end state parent_state parent_node fg_change fg_type
#>    <chr>    <chr>  <dbl>  <dbl> <dbl>        <dbl> <chr>           <dbl> <chr>  
#>  1 AT21352… 8     1.31e8 1.31e8     8            7 internal_2…         1 fg-gain
#>  2 AT21352… 5     7.30e7 7.35e7     3            2 internal_2…         1 fg-gain
#>  3 AT21352… 14    9.05e7 9.1 e7     4            4 internal_2…         0 NA     
#>  4 AT21352… 4     8.15e7 8.2 e7     4            4 internal_2…         0 NA     
#>  5 AT21352… 17    1.75e7 1.80e7     8            7 internal_2…         1 fg-gain
#>  6 AT21352… 8     5.15e7 5.20e7     7            6 internal_2…         1 fg-gain
#>  7 AT21352… 15    4.90e7 4.95e7     4            3 internal_2…         1 fg-gain
#>  8 AT21352… 13    6.70e7 6.75e7     3            2 internal_64         1 fg-gain
#>  9 AT21352… 17    7.45e7 7.5 e7     4            5 internal_2…        -1 fg-loss
#> 10 AT21352… 3     2.00e7 2.05e7     3            4 internal_2…        -1 fg-loss
#> # ℹ 2 more variables: foreground <lgl>, background <lgl>
```

## Plotting Foreground

We can visualize these inference in our classic heatmaps, highlighting
differnt features.

First, as a reference point the raw states:

``` r

# subsetting chromosomes for plotting purposes
states_sub <- dplyr::filter(states_df, chr %in% c(3, 4, 7))
```

``` r

dlptools::plot_state_hm(
  states_sub,
  state_col = "state",
  phylogeny = med_tree,
  file_name = "imgs/ex_med_raw_states.png"
)
```

![](imgs/ex_med_raw_states.png)

The raw states of the inferred foreground events:

``` r

dlptools::plot_fg_state_highlight(
  states_df = states_sub,
  phylogeny = med_tree,
  file_name = "imgs/medicc_fg_highlight.png"
)
```

![](imgs/medicc_fg_highlight.png)

The raw states of the implied background:

``` r

dlptools::plot_bg_state_highlight(
  states_df = states_sub,
  phylogeny = med_tree,
  file_name = "imgs/medicc_bg_highlight.png"
)
```

![](imgs/medicc_bg_highlight.png)

And the states of foreground events themselves, i.e., the change from
the parent node to tip:

``` r

dlptools::plot_heatmap_of_tip_changes(
  states_df = states_sub,
  phylogeny = med_tree,
  file_name = "imgs/medicc_fg_change.png"
  # can specify the name of the column, if you didn't use the above functions
  # to generate. Basically any column that is an integer spanning - to + changes.
  # changes_col = "fg_change"
)
```

Currently there is a limit to +/- 8 state changes between parent node
and tip, simply because I couldn’t find more colors that seemed to work
well in the palette. So anything bigger than that is capped at +/- 8.
Jumps beyond this are likely rare.

![](imgs/medicc_fg_change.png)

  
  
  
  
