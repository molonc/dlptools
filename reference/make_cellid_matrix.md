# pivot a cell_id dataframe to a wide matrix.

Different functions may produce a dataframe of counts or something per
cell_id for various attributes. This converts it to a matrix with
rownames of cell_id and columns of some feature of interest from the
dataframe that can be pivoted out.

## Usage

``` r
make_cellid_matrix(cell_df, name_col, val_col)
```

## Arguments

- cell_df:

  some dataframe with a column of cell_id and values for each of those
  cells.

- name_col:

  the column to be pivoted out

- val_col:

  the column to use as values in the matrix

## Value

matrix

## Examples

``` r
cell_df <- data.frame(
  cell_id = c("cell_1", "cell_1", "cell_1"),
  some_col = c("A", "B", "C"),
  some_val = 1:3
)

make_cellid_matrix(cell_df, name_col = "some_col", val_col = "some_val")
#>        A B C
#> cell_1 1 2 3

#       A B C
# cell_1 1 2 3
```
