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
