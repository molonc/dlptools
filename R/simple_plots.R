#' create a tile plot of read state calls.
#'
#' builds a basic ggplot with geom_tile on a reads df.
#' Only expects columns of cell_id, start, state, chr
#'
#' @param reads_df a table of reads data (e.g., could load with
#' [import_dlp_files()])
#' @return ggplot object
#' @export
#' @importFrom rlang .data
basic_tile_plot <- function(reads_df) {
  tile_p <- ggplot2::ggplot(
    reads_df,
    ggplot2::aes(
      .data$start, .data$cell_id,
      fill = base::as.factor(.data$state)
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(
      values = dlptools::CNV_COLOURS, "CNV", na.value = "white"
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(.data$chr),
      scales = "free", space = "free", switch = "x"
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL, "Chromosome") +
    ggplot2::theme(
      panel.spacing = ggplot2::unit(0.1, "lines"),
      strip.background = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
    ) +
    ggplot2::ylab("Cell ID") +
    ggplot2::xlab("Chromosome")

  return(tile_p)
}

#' Plot metrics across the chip
#'
#' Map of a metric as observed in the layout of the DLP chip.
#'
#' @param metrics_df dataframe of metrics for a run.
#' @param targ_val the value of interest for plotting. E.g.,
#' experimental_condition, quality, total_reads.
#' @return ggplot2 object
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

#' visualize GC correction for a cell
#'
#' Seeing how GC correction worked for a cell can give an indication if the
#' multiplier chosen by HMMCopy was appropriate to form integer copy numbers.
#' Also can just be informative to see what's going on behind the float copy
#' number to integer CN value.
#'
#' @param reads_df dataframe of reads like information.
#' @param cellid string of cell to plot
#' @param plot_choice string of plots to show. Default is both, alternatively "
#' 'raw' for just the raw GC values or 'corrected' for just the corrected
#' values.
#' @return ggplot2 object
#' @export
gc_plot <- function(
    reads_df, cellid, plot_choice = c("both", "raw", "corrected")) {
  plot_choice <- match.arg(plot_choice)

  cell_reads <- dplyr::filter(reads_df, cell_id == cellid)
  cell_reads <- dplyr::mutate(cell_reads, gc = dplyr::if_else(gc == -1, NA, gc))

  cell_reads <- dplyr::mutate(
    cell_reads,
    corrected_reads = reads / modal_curve * multiplier
  )

  reads_quant <- quantile(cell_reads$reads, 0.98)
  cor_reads_quant <- quantile(cell_reads$corrected_reads, 0.98, na.rm = TRUE)

  plot_stuff <- list(
    ggplot2::geom_point(size = 0.5),
    ggplot2::scale_color_manual(values = CNV_COLOURS),
    ggplot2::guides(color = ggplot2::guide_legend(title = "State")),
    ggplot2::xlab("GC"),
    ggplot2::theme_classic()
  )

  gc_p <- ggplot2::ggplot(
    cell_reads, ggplot2::aes(x = gc, y = reads, color = factor(state))
  ) +
    ggplot2::coord_cartesian(ylim = c(0, reads_quant)) +
    plot_stuff +
    ggplot2::geom_line(ggplot2::aes(y = modal_curve), color = "red") +
    ggplot2::ylab("Reads")

  gc_corr_p <- ggplot2::ggplot(
    cell_reads, ggplot2::aes(x = gc, y = corrected_reads, color = factor(state))
  ) +
    ggplot2::geom_hline(
      data = tibble::tibble(cn_ints = 1:round(max(cor_reads_quant))),
      ggplot2::aes(yintercept = cn_ints),
      color = "green"
    ) +
    plot_stuff +
    ggplot2::coord_cartesian(ylim = c(0, cor_reads_quant)) +
    ggplot2::scale_y_continuous(breaks = 1:round(max(cor_reads_quant))) +
    ggplot2::ylab("Corrected Reads")

  if (plot_choice == "both") {
    final_p <- patchwork::wrap_plots(gc_p, gc_corr_p, nrow = 1) +
      patchwork::plot_layout(guides = "collect")
  } else if (plot_choice == "raw") {
    final_p <- gc_p
  } else if (plot_choice == "corrected") {
    final_p <- gc_corr_p
  }

  final_p <- final_p + patchwork::plot_annotation(title = cellid)
  return(final_p)
}


#' plot copy number profile plot of a cell.
#'
#' Basically a geom_point of the float copy number, colored by the integer
#' state number. singals::plotCNprofile() is similar, possibly with more
#' features. But this function works in a pinch.
#'
#' @param reads_df dataframe of read bins like information.
#' @param cell_ids string/vector of cell ids to plot. Remember, too many and
#' the plot will be useless. If blank, will plot all cells in dataframe, but
#' only if less than 10, because honestly at that point, you probably need a
#' different plot.
#' @param pseudo_log_y boolean. Pseudo-log transform the y axis, as sometimes
#' high copy numbers obscure the plot. Or you can always set a
#' ggplot2::coord_cartesian() yourself on the plot returned by this function.
#' @return ggplot2 object
#' @export
cell_cn_profile <- function(reads_df, cell_ids = c(), pseduo_log_y = FALSE) {
  if (length(cell_ids) == 0) {
    n_cell_ids <- unique(reads_df$cell_id)
    if (length(n_cell_ids) > 10) {
      stop(paste0(
        "There are more than 10 cell ids in this dataframe. You need a",
        " different plot, 'casuse I won't do it!"
      ))
    }
    cells_df <- reads_df
  } else {
    cells_df <- dplyr::filter(reads_df, cell_id %in% cell_ids)
  }

  cells_df$chr <- factor_column_mixedsort(cells_df, "chr")

  cell_p <- ggplot2::ggplot(
    cells_df, ggplot2::aes(start, copy, col = as.factor(state))
  ) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::facet_grid(
      cell_id ~ chr,
      scales = "free",
      space = "free_x",
      switch = "x",
      axes = "all_x",
    ) +
    ggplot2::scale_colour_manual(values = CNV_COLOURS, "CNV") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(0.1, "lines"),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(),
      strip.placement = "outside",
      strip.text.y = ggplot2::element_text(angle = 0),
      legend.position = "bottom"
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
    ) +
    ggplot2::xlab("Chromosome")

  if (pseduo_log_y) {
    cell_p <- cell_p +
      ggplot2::scale_y_continuous(transform = scales::pseudo_log_trans())
  }

  return(cell_p)
}
