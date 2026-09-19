# alter existing segments to fit new breakpoints

alter existing segments to fit new breakpoints

## Usage

``` r
snap_segs_to_new_breakpoints(jfit, input_mtx, input_cellids)
```

## Arguments

- jfit:

  a
  [`jointseg::jointSeg()`](https://rdrr.io/pkg/jointseg/man/jointSeg.html)
  output object

- input_mtx:

  matrix. cell by bin state matrix

- input_cellids:

  vector. Cell ids of cells in the matrix.

## Value

tibble. Cell segments that align with new breakpoints
