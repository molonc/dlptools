# builds the left-side annotations of the cells

builds the left-side annotations of the cells

## Usage

``` r
build_left_annot(
  anno_df = NULL,
  anno_cols_list = list(),
  annotation_legend_param = list(),
  annotation_name_gp = 16,
  show_annotation_legend = TRUE,
  clones_df = NULL,
  only_largest_clone_group = FALSE,
  labels_fontsize = 8,
  clone_palette = NULL
)
```

## Arguments

- anno_df:

  annotations dataframe with cell_id column and annotation for each cell
  id to be added to a heatmap.

- anno_cols_list:

  list of named vectors specifying colors for annotations example:
  list(passage=c(`1`='#2872bc', `19`='#d23e3e'))

- annotation_name_gp:

  Int: Fontsize for label names

- clones_df:

  dataframe of clone ideas (clone_id) for each cell_id. Both columns
  required.

- labels_fontsize:

  how large to make text labels

- clone_palette:

  named vector of colors to give to clones. E.g., c( A='#12345',
  B='#67890'). See dlptools::make_clone_palette()
