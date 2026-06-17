# main hm building function

Plot a heat map of states across chromosomes. Put annotations next to
it, a tree next to it, mark clone groups, or all of that at once!

## Usage

``` r
plot_state_hm(
  states_df,
  state_col,
  phylogeny = NULL,
  file_name = NULL,
  anno_df = NULL,
  anno_colors_list = list(),
  anno_columns = NULL,
  clone_column = NULL,
  clones_df = NULL,
  clone_palette = NULL,
  only_largest_clone_group = TRUE,
  labels_fontsize = 8,
  is_continuous = FALSE,
  continuous_hm_colours = FALSE,
  custom_continuous_colors = NULL,
  custom_continuous_range = NULL,
  hm_discrete_colors = NULL,
  hm_legend_title = NULL,
  legend_11plus = FALSE,
  color_tree_clones = FALSE,
  scale_y = FALSE,
  cell_height = 10,
  default_png_height = 1600,
  annotation_legend_param = list(),
  show_annotation_legend = TRUE,
  annotation_name_gp = 16,
  heatmap_legend_param = list(),
  show_heatmap_legend = TRUE,
  custom_legend = list(),
  hm_title = NULL,
  ...
)
```

## Arguments

- states_df:

  long format read bin data to be plotted

- state_col:

  string of column name to target for plotting in the heatmap. Examples
  include: state, BAF, state_AS, state_phase

- phylogeny:

  optional. phylo class object to be plotted.

- file_name:

  name of the file to save the plot to. Recommended for most cases as
  plots are big-ish.

- anno_df:

  optional. annotations dataframe with cell_id column and annotation for
  each cell id to be added to a heatmap.

- anno_colors_list:

  list of named vectors specifying colors for annotations example:
  `list(passage=c(p1='#2872bc', p19='#d23e3e'))`

- anno_columns:

  optional. Columns containing the annotation data to plot.

- clone_column:

  optional. Column of clone id labels for cells.

- clones_df:

  optional. dataframe of clone ideas (clone_id) for each cell_id. Both
  columns required.

- clone_palette:

  optional. ideally a named vector on how to color clones. Vector names
  are clond ID, values are hex codes. If vector is unnamed, names will
  be assigned, but without control on which get's which color.

- only_largest_clone_group:

  boolean. optional. Only put a letter label on the largest group of any
  given clone id.

- labels_fontsize:

  how large to make text labels

- continuous_hm_colours:

  plot heatmap colors on a continous scale.

- custom_continuous_colors:

  a vector of 3 colors to use as the scale for plotting a continuous
  variable, specified as hexcodes. E.g., c("#3182BD", "#CCCCCC",
  "#FDCC8A")

- custom_continuous_range:

  a vector of values to specify as the low, mid, and high bounds for the
  continuous color scale, e.g., c(1, 5, 10)

- hm_discrete_colors:

  specified colors for the values being plotted. Named vector:
  c(1="#3182BD", 2="#FDCC8A"). Need to specify for every value

- hm_legend_title:

  string. What to label the legend for what is plotted in the heatmap.
  Will default to the name of the column being plotted.

- legend_11plus:

  bool. Default FALSE. For HMMCopy state values, state 11 is really 11+,
  so we replace 11s with 11+ for the plot in the \#' legend.

- color_tree_clones:

  boolean. optional. Whether to color the tree with the same colors as
  the clone labels.

- scale_y:

  bool. Default FALSE. Scale the y size of the saved image accoring to
  the number of cells in the heatmap

- cell_height, :

  int. Default = 10 Heigh in pixels per cell if scale_y is True

- default_png_height, :

  int: Default = 1600. If scale_y not used, height of final png

- annotation_legend_param, :

  list(). Parameters topass to the annotation legend ex:
  annotation_legend_params \<- list( Small_inter_few_telo = list(
  title_gp = grid::gpar(fontsize = 20), \#fontsize of title labels_gp =
  grid::gpar(fontsize = 16), \# fontsize of labels gp =
  grid::gpar(fontsize = 16) \# siez of color bar or whatever legend is
  showing I think )

- show_annotation_legend, :

  bool.

- annotation_name_gp:

  Int: Fontsize for label names

- show_heatmap_legend, :

  bool

- custom_legend, :

  List(). custom legend for use when annotation and heatmap legends are
  FALSE For examples, see
