#' Plot metrics across the chip
#'
#' Map of a metric as observed in the layout of the DLP chip.
#'
#' @param metrics_df dataframe of metrics for a run.
#' @param targ_val the value of interest for plotting. E.g.,
#' experimental_condition, quality, total_reads.
#' @export
chip_plot <- function(metrics_df, targ_val) {
  # boundaries of experimental conditions for outlining on plot
  exp_bounds <- metrics_df |>
    dplyr::filter(!is_control) |>
    dplyr::group_by(experimental_condition) |>
    dplyr::summarise(
      min_row = min(row),
      max_row = max(row),
      min_column = min(column),
      max_column = max(column)
    ) |>
    dplyr::ungroup()

  control_bounds <- metrics_df |>
    # hTERTs are standard and in multiple places, so just drop from labelling
    dplyr::filter(is_control & experimental_condition != "hTERT") |>
    dplyr::group_by(experimental_condition, column) |>
    dplyr::summarise(
      max_row = max(row)
    )

  # set up some colors for the experimental condition outlines
  non_control_exps <- unique(exp_bounds$experimental_condition)
  control_exps <- dplyr::filter(metrics_df, is_control) |>
    dplyr::pull(experimental_condition) |>
    unique()

  exp_colors <- GEN_PLOT_COLS
  names(exp_colors) <- c(non_control_exps, control_exps)

  chip_p <- ggplot2::ggplot(metrics_df) +
    ggplot2::geom_tile(
      ggplot2::aes(x = column, y = row, fill = .data[[targ_val]])
    ) +
    ggplot2::geom_rect(
      data = exp_bounds,
      ggplot2::aes(
        xmin = min_column - 0.25, xmax = max_column + 0.25,
        ymin = min_row - 0.25, ymax = max_row + 0.25,
        color = experimental_condition
      ),
      fill = NA,
      linewidth = 1.5
    ) +
    ggplot2::geom_label(
      data = exp_bounds,
      ggplot2::aes(
        x = (max_column + min_column) / 2,
        y = max_row - 1,
        label = experimental_condition
      )
    ) +
    ggrepel::geom_label_repel(
      data = control_bounds,
      ggplot2::aes(x = column, y = max_row, label = experimental_condition),
      min.segment.length = 0
    ) +
    BASE_PLOT_THEME +
    ggplot2::scale_color_manual(values = exp_colors) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank()
    ) +
    ggplot2::ggtitle("Chip Layout")

  if (targ_val == "experimental_condition") {
    # so the legend doens't double label
    chip_p <- chip_p +
      ggplot2::guides(color = "none")
  }

  if (is.numeric(metrics_df[[targ_val]])) {
    chip_p <- chip_p +
      ggsci::scale_fill_tw3("purple")
  } else {
    chip_p <- chip_p +
      ggplot2::scale_fill_manual(values = exp_colors)
  }
  return(chip_p)
}
