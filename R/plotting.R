#' standard colors used in dlp plots
#' @export
CNV_COLOURS <- structure(
  c(
    "#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
    "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA",
    "#D4B9DA"
  ),
  names = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+", "11")
)


#' create a tile plot of read state calls.
#'
#' builds a basic ggplot with geom_tile on a reads df.
#' Only expects columns of cell_id, start, state, chr
#'
#' @param reads_df a table of reads data (e.g., could load with
#' [import_dlp_files()])
#' @return ggplot object
#' @export
basic_tile_plot <- function(reads_df) {
  tile_p <- ggplot2::ggplot(
    reads_df,
    ggplot2::aes(start, cell_id, fill = base::as.factor(state))
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(
      values = dlptools::CNV_COLOURS, "CNV", na.value = "white"
    ) +
    ggplot2::facet_grid(~chr, scales = "free", space = "free", switch = "x") +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL, "Chromosome") +
    ggplot2::theme(panel.spacing = ggplot2::unit(0.1, "lines"))

  return(tile_p)
}
