# Internal function. Creates a complex heatmap of a given matrix of states.

See
[`plot_state_hm()`](https://molonc.github.io/dlptools/reference/plot_state_hm.md)
for arguments description.

## Usage

``` r
generate_state_hm(
  states_mat,
  phylo,
  labels_fontsize = 8,
  plot_cols = STATE_COLORS,
  left_annot = NULL,
  hm_legend_title = NULL,
  legend_11plus = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list()
)
```
