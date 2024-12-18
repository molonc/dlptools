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

#' alias of standard colors used in dlp plots
#' @export
STATE_COLORS <- CNV_COLOURS

#' min, median, max for a continuous color range
DEFAULT_CONTINUOUS_COLOR_RANGE <- c(
  "#3182BD",
  "#CCCCCC",
  "#FDCC8A"
)

CONTINUOUS_COLOR_RANGE_ALT <- c("#000000", "#ffffff", "#5F9EA0")

#' colors for signals results of allele balances
#' @export
ASCN_COLORS <- c(
  `0|0` = "#3182BD",
  `1|0` = "#9ECAE1",
  `1|1` = "#CCCCCC",
  `2|0` = "#666666",
  `2|1` = "#FDCC8A",
  `3|0` = "#FEE2BC",
  `2|2` = "#FC8D59",
  `3|1` = "#FDC1A4",
  `4|0` = "#FB590E",
  `5` = "#E34A33",
  `6` = "#B30000",
  `7` = "#980043",
  `8` = "#DD1C77",
  `9` = "#DF65B0",
  `10` = "#C994C7",
  `11` = "#D4B9DA",
  `11+` = "#D4B9DA"
)


#' ASCN phase colors
#' @export
ASCN_PHASE_COLORS <- c(
  `A-Hom` = "#56941E",
  `B-Hom` = "#471871",
  `A-Gained` = "#94C773",
  `B-Gained` = "#7B52AE",
  `Balanced` = "#d5d5d4"
)

#' @export
BAF_COLORS <- circlize::colorRamp2(
  c(0, 0.5, 1),
  c(
    ASCN_PHASE_COLORS["A-Hom"],
    ASCN_PHASE_COLORS["Balanced"],
    ASCN_PHASE_COLORS["B-Hom"]
  )
)

#' picks beginning, middle, and end of a vector
#' to handle when vectors that are too long are
#' passed
bme_vec <- function(vec, vec_name) {
  warning(paste0(
    "more than 3 values given for ", vec_name,
    " Only using the first, middle, and last value."
  ))
  red_vec <- vec[c(
    1,
    ceiling(length(vec) / 2),
    length(vec)
  )]
  return(red_vec)
}


#' internal function for setting up heatmap continuous range colors
#' chooses defaults, unless overwritten by user.
fetch_continuous_color_ramp <- function(
    plotting_values, custom_continuous_colors = NULL, custom_continuous_range = NULL) {
  if (is.null(custom_continuous_range)) {
    metrics <- get_column_metrics(plotting_values)
  } else {
    if (length(custom_continuous_range) > 3) {
      custom_continuous_range <- bme_vec(custom_continuous_range, "custom_continuous_range")
    }
    metrics <- custom_continuous_range
  }

  if (!is.null(custom_continuous_colors)) {
    if (length(custom_continuous_colors) > 3) {
      custom_continuous_colors <- bme_vec(custom_continuous_colors, "custom_continuous_colors")
    } else if (length(custom_continuous_colors) < 3) {
      stop(paste0(
        "need to specify 3 colors (bottom, middle, top) for custom",
        " continuous color palette."
      ))
    }

    continuous_color_bounds <- custom_continuous_colors
  } else {
    continuous_color_bounds <- DEFAULT_CONTINUOUS_COLOR_RANGE
  }

  continuous_color_palette <- circlize::colorRamp2(
    sort(metrics),
    continuous_color_bounds
  )
}
