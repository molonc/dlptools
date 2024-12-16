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
