#' metric of pairwise differences between two cells.
#'
#' Inspired by MSKCC SPECTRUM paper. Bins are aligned between two cells and
#' marked for if they have the same state Then segments of matching and
#' non-matching runs of bins are found (filtering those smaller than a
#' specified minimum). These segments are then re-split into 500kb bins and the
#' difference becomes the number of matching bins divided by the number of
#' considered bins.
#'
#' This function is slow, and the number of pairwise comparisons grows quickly.
#' Dramatic speed improvements can be had by setting up a parallel plan for
#' furrr like so:
#'
#' future::plan(future::multicore, workers=N_CORES_YOU_WANT)
#'
#' Example for 100 cells, which is 4950 pairs, this function will take 2 minutes
#' with 4 cores.
#'
#' @param bin_df a dataframe of read bins with states. Expected columns of:
#' cell_id, chr, start, end, state
#' @param cells optional vector specifying cells to compare. If it's blank, all
#' cells are compared. If it's 1 cell, then that one cell is compared to all
#' others. If it's 2 or more, then just the specified cells are compared to each
#' other.
#' @param min_seg_length double. This is the minium length of matching segment
#' bins to use when measuring similarity.
#' @return tibble of cell pairs and metrics about their differences.
#' @export
pairwise_bin_difference <- function(
    bin_df, cells = c(), min_seg_length = 2.5e6) {
  required_cols <- c("cell_id", "chr", "start", "end", "state")
  if (!all(required_cols %in% bin_df)) {
    stop(paste(
      "Need columns of: ",
      paste(required_cols, collapse = ", "),
      "\n Or someone needs to make this function more flexible."
    ))
  }

  if (is(future::plan(), "sequential")) {
    warning(paste(
      "This function will improve in speed dramatically if you set up a",
      "parallel plan for furrr. Use:\n\n",
      "future::plan(future::multicore, workers=N_CORES_YOU_WANT)",
      "\n\n",
      "ctrl-C to kill this and start over.",
      "\n\n",
      sep = " "
    ), immediate. = TRUE)
  }

  targ_cell <- NULL
  if (length(cells) == 0) {
    print("comparing all cells. Gonna take some time.")
    cells <- unique(bin_df$cell_id)
  } else if (length(cells) == 1) {
    print("comparing one cell to all others")
    targ_cell <- cells[1]
    print(targ_cell)
    cells <- unique(bin_df$cell_id)
  }

  cell_pairs <- t(combn(unique(cells), 2))

  if (!is.null(targ_cell)) {
    cell_pairs <- cell_pairs[
      c(
        cell_pairs[, 1] == targ_cell |
          cell_pairs[, 2] == targ_cell
      ), ,
      drop = FALSE
    ]
  }

  print(stringr::str_c(
    "processing:", nrow(cell_pairs), "pairs",
    sep = " "
  ))

  targ_cells_df <- dplyr::filter(bin_df, cell_id %in% cells)
  targ_cells_df <- split(targ_cells_df, f = as.factor(targ_cells_df$cell_id))

  all_cell_comps <- furrr::future_map2(
    cell_pairs[, 1], cell_pairs[, 2],
    \(c1, c2) {
      compare_two_cells(
        targ_cells_df[[c1]], targ_cells_df[[c2]],
        min_seg_length = min_seg_length
      )
    }
  ) |>
    purrr::list_rbind()

  return(all_cell_comps)
}

#' internal workhorse of pairwise_bin_difference function
#'
#' @param cell_1_df dataframe of bin state data for one cell
#' @param cell_2_df dataframe of bin state data for the other cell
#' @param min_seg_length int/double. minimum length of matching/non-matching
#' runs of bins to consider in the calculation of the differences between the
#' cells.
#' @return 1 row tibble summarizing differences.
compare_two_cells <- function(cell_1_df, cell_2_df, min_seg_length) {
  # first align bins
  targ_cols <- c("cell_id", "chr", "start", "end", "state")
  comp <- dplyr::inner_join(
    dplyr::select(cell_1_df, dplyr::all_of(targ_cols)),
    dplyr::select(cell_2_df, dplyr::all_of(targ_cols)),
    by = c("chr", "start", "end"),
    suffix = c("_cell1", "_cell2")
  )

  # now basically want to make segments of matches and not matches
  diff_res <- comp |>
    dplyr::group_by(.data$chr) |>
    dplyr::mutate(
      bin_match = state_cell1 == state_cell2,
      rle_matches = rle_states(bin_match)
    ) |>
    dplyr::group_by(.data$chr, .data$rle_matches) |>
    dplyr::summarise(
      start = base::min(.data$start),
      end = base::max(.data$end),
      span = end - start,
      match = unique(bin_match),
      .groups = "drop"
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(
      span >= min_seg_length
    ) |>
    segs_to_reads() |>
    dplyr::summarise(
      n_diff = sum(match == FALSE),
      tot_bins = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      prop_diff = n_diff / tot_bins,
      cell_one = unique(cell_1_df$cell_id)[1],
      cell_two = unique(cell_2_df$cell_id)[1]
    )

  return(diff_res)
}
