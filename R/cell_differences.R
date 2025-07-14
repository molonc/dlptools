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
#' The returned DF is organized by each cell and the distances to each other
#' cell (so there are some redundant comparisons, like cell 1 vs cell 2 and
#' cell 2 vs cell 1). There is also a column "nearest_neighbour" which is a
#' boolean identifying which comparison is the minimum distance for each cell.
#'
#' @param bin_df a dataframe of read bins with states. Expected columns of:
#' cell_id, chr, start, end, state
#' @param cells optional vector specifying cells to compare. If it's blank, all
#' cells are compared. If it's 1 cell, then that one cell is compared to all
#' others. If it's 2 or more, then just the specified cells are compared to each
#' other.
#' @param min_seg_length double. This is the minium length of matching segment
#' bins to use when measuring similarity.
#' @param return_pairs_matrix boolean. If TRUE, returns a pairwise matrix object of distances. This is useful to then pass to functions like hclust() and so forth. Can also do afterwards with dlptools::convert_dists_to_pairwise()
#' @return tibble of cell pairs and metrics about their differences.
#' @export
pairwise_bin_difference <- function(
    bin_df, cells = c(), min_seg_length = 2.5e6, return_pairs_matrix = FALSE) {
  required_cols <- c("cell_id", "chr", "start", "end", "state")
  if (!all(required_cols %in% colnames(bin_df))) {
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

  # this comparison was done in a non-redundant way (i.e., only doing 1 vs 2
  # and not also the converse of 2 vs 1). But want to reorg to see each cell
  # against all of the others.
  cell_based_comps <- furrr::future_map(
    cells,
    \(cell) {
      all_cell_comps |>
        dplyr::filter(cell_one == cell | cell_two == cell) |>
        dplyr::mutate(
          index_cell = cell,
          comp_cell = dplyr::if_else(cell_one == cell, cell_two, cell_one),
          nearest_neighbour = prop_diff == min(prop_diff)
        ) |>
        dplyr::select(-c(cell_one, cell_two))
    }
  ) |>
    purrr::list_rbind()

  if (return_pairs_matrix) {
    return(convert_dists_to_pairwise(cell_based_comps))
  } else {
    return(cell_based_comps)
  }
}


#' convert cell distances to a pairwise matrix
#'
#' Takes the output of dlptools::pairwise_bin_difference() and returns a
#' pairwise matrix of the distances. Useful to then pass to functions like
#' stats::hclust() and related ideas.
#'
#' @param cell_dists dataframe from dlptools::pairwise_bin_difference()
#' @return matrix of index vs comp cell distances.
#' @export
convert_dists_to_pairwise <- function(cell_dists) {
  p_mtx <- xtabs(
    prop_diff ~ .,
    dplyr::select(cell_dists, cell_one, cell_two, prop_diff)
  ) |> t()

  return(p_mtx)
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


#' Find outlier cells using a beta distribution
#'
#' Also inspired by MSKCC SPECTRUM paper. Using the measured pairwise bin
#' distances between cells (dlptools::pairwise_bin_differences()), this
#' function takes the nearest neighbour of each cell and fits a beta
#' distribution. Then using this distribution, it finds outlier cells based on
#' a selected percentile of the distribution (default 99th).
#'
#' @param cell_diffs the dataframe of differences from
#' dlptools::pairwise_bin_differences()
#' @param outlier_percentile double. Default 0.99. What percentile of the
#' distribution to consider an outlier cell.
#' @return tibble of information on cells considered outliers.
#' @export
find_outlier_cells <- function(cell_diffs, outlier_percentile = 0.99) {
  if (!"nearest_neighbour" %in% colnames(cell_diffs)) {
    stop(paste(
      "Requries a boolean column with the name: nearest_neighbour. First need",
      "to run dlptools::pairwise_bin_difference() to get the nearest_neighbour",
      "values for the cells.",
      sep = " "
    ))
  }

  # maybe do the average, I am capturing some interesting cells with a low
  # distance to one another, but high distance to everyone else.
  nn_cells <- dplyr::filter(cell_diffs, nearest_neighbour)
  diff_beta_dist <- fitdistrplus::fitdist(nn_cells$prop_diff, "beta")

  min_diff <- qbeta(
    p = outlier_percentile,
    shape1 = diff_beta_dist$estimate[["shape1"]],
    shape2 = diff_beta_dist$estimate[["shape2"]]
  )

  outlier_cells <- dplyr::filter(nn_cells, prop_diff >= min_diff)

  outlier_cell_info <- cell_diffs |>
    dplyr::filter(index_cell %in% outlier_cells$index_cell) |>
    dplyr::rename(outlier_cell = index_cell) |>
    dplyr::group_by(outlier_cell) |>
    dplyr::summarise(
      mean_diff_to_all_cells = mean(prop_diff),
      nn_dist = unique(prop_diff[nearest_neighbour]),
      nn_cell = unique(comp_cell[nearest_neighbour])
    )

  return(outlier_cell_info)
}

#' visualize cell nearest neighbour distance values
#'
#' Generates a simple plot of each cell and highlights the min NN-distance for
#' each. Also highlights cells that are considered outliers.
#'
#' @param cell_diffs output dataframe of dlptools::pairwise_bin_difference()
#' @param outlier_Cells output dataframe of dlptools::find_outlier_cells()
#' @return ggplot object
#' @export
plot_nnd_outlier_cells <- function(cell_diffs, outlier_cells) {
  p_dat <- cell_diffs |>
    dplyr::group_by(index_cell) |>
    dplyr::mutate(
      outlier = index_cell %in% outlier_cells,
      p_lab = dplyr::case_when(
        prop_diff == min(prop_diff) & !outlier ~ "min-NN",
        prop_diff == min(prop_diff) & outlier ~ "min-NN-outlier",
        .default = "other"
      )
    ) |>
    dplyr::ungroup()

  cell_ord <- p_dat |>
    dplyr::filter(nearest_neighbour) |>
    dplyr::arrange(prop_diff) |>
    dplyr::pull(index_cell) |>
    unique()

  p_dat$index_cell <- factor(p_dat$index_cell, levels = cell_ord)

  ggplot(p_dat, aes(x = index_cell, y = prop_diff, color = p_lab)) +
    geom_point() +
    # facet_grid(~outlier, scales='free_x') +
    theme_classic() +
    ggplot2::theme(
      plot.background = element_rect(color = "white"),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    scale_color_manual(
      values = c(
        `min-NN` = GEN_PLOT_COLS[2],
        `min-NN-outlier` = GEN_PLOT_COLS[3],
        `other` = "grey"
      )
    ) +
    xlab("Index Cell") +
    ylab("Nearest Neighbour Distance") +
    guides(color = guide_legend(title = ""))
}
