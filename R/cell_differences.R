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
#' future::plan(future::multisession, workers=N_CORES_YOU_WANT)
#'
#' Example for 100 cells, which is 4950 pairs, this function will take 4 minutes
#' with 4 cores.
#'
#' These comparisons are unique pairs! So if 3 input cells: A, B, C,
#' comparisons made are A-B, A-C, B-C. So one input cell will always be
#' "missing" from the index column. See
#' [the vignette](https://molonc.github.io/dlptools/articles/pairwise-differences.html)
#' for how you can turn this around, but it can create a very large DF with
#' redundant comparisons to rearrange to have each cell against all others.
#'
#' @param bin_df a dataframe of read bins with states. Expected columns of:
#' cell_id, chr, start, end, state
#' @param targ_cells optional vector specifying which cells to compare to all
#' other cells.
#' @param min_seg_length double. This is the minium length of matching segment
#' bins to use when measuring similarity.
#' @return tibble of cell pairs and metrics about their differences.
#' @export
pairwise_bin_difference <- function(
  bin_df,
  targ_cells = c(),
  min_seg_length = 2.5e6,
  return_pairs_matrix = FALSE
) {
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

  cells <- unique(bin_df$cell_id)
  cell_pairs <- t(combn(unique(cells), 2))

  if (length(targ_cells) > 0) {
    # remove unneeded comparisons
    cell_pairs <- cell_pairs[
      c(
        cell_pairs[, 1] %in% targ_cells |
          cell_pairs[, 2] %in% targ_cells
      ), ,
      drop = FALSE
    ]
  }

  warning(stringr::str_c(
    "processing:", nrow(cell_pairs), "pairs",
    sep = " "
  ), immediate. = TRUE)


  bind_df_split <- split(bin_df, f = as.factor(bin_df$cell_id))

  all_cell_comps <- furrr::future_map2(
    cell_pairs[, 1], cell_pairs[, 2],
    \(c1, c2) {
      compare_two_cells(
        bind_df_split[[c1]], bind_df_split[[c2]],
        min_seg_length = min_seg_length
      )
    }
  ) |>
    purrr::list_rbind()


  all_cell_comps <- all_cell_comps |>
    dplyr::rename(
      index_cell = cell_one,
      comp_cell = cell_two,
    )


  return(all_cell_comps)
}


#' mostly internal for rearranging the pairwise DF to focus on a cell.
#'
#' This is useful because with non-redundant pair comparisons, the "last cell"
#' in the list is never the "index_cell", as it has already been compared ot all
#' others. But sometimes useful to still see it.
#'
#' e.g.,
#' cells: A, B, C
#' comparisons would be: A-B, A-C, B-C
#' generate C focus with this function: C-A, C-B
#' or B focus: B-A, B-C
#' @param in_df dataframe most likely from [dlptools::pairwise_bin_difference()]
#' @param cell string. Target cell ID to put as index cell
#' @return tibble
#' @export
make_cell_focused_matrix <- function(in_df, cell) {
  in_df |>
    dplyr::filter(index_cell == cell | comp_cell == cell) |>
    dplyr::mutate(
      comp_cell = dplyr::if_else(
        index_cell == cell, comp_cell, index_cell
      ),
      index_cell = cell
    )
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
    dplyr::mutate(
      fake_cell_col = "placeholder"
    ) |>
    segs_to_reads(
      sample_id_col = "fake_cell_col",
      other_meta_cols = c("match"),
      bin_size = 5e5
    ) |>
    dplyr::select(-fake_cell_col) |>
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

#' find the nearest neighbour of each cell
#'
#' After running [dlptools::pairwise_bin_difference()], we can feed that output
#' dataframe to this function to find the nearest neighbour of each cell. The
#' pairwise function does the comparisons in a non-redundant way (A-B, A-C,
#' B-C), leaving out some redundant comparisons (i.e., C-A, C-B). If A's
#' nearest neighbour is C, it's not guaranteed that C's nearest neighbour is A.
#' This function takes a cell specific focus and finds the nearest neighbour of
#' each.
#'
#' As with [dlptools::pairwise_bin_difference()], using [furrr::future_map()]
#' internally, so speed improvements with a [future::plan()] being set.
#'
#' @param pairwise_diffs tibble output of [dlptools::pairwise_bin_difference()]
#' @return tibble/dataframe
#' @export
find_nearest_neighbours <- function(pairwise_diffs) {
  cells <- unique(c(pairwise_diffs$index_cell, pairwise_diffs$comp_cell))

  nearest_neighs <- furrr::future_map_dfr(
    cells, \(cell) {
      cell_df <- make_cell_focused_matrix(pairwise_diffs, cell)

      cell_df |>
        dplyr::slice_min(prop_diff) |>
        dplyr::rename(nn_diff = prop_diff) |>
        dplyr::mutate(
          max_diff_to_all = max(cell_df$prop_diff),
          mean_diff_to_all = mean(cell_df$prop_diff)
        )
    }
  )

  return(nearest_neighs)
}


#' Find outlier cells using a beta distribution
#'
#' Also inspired by MSKCC SPECTRUM paper. Using the measured pairwise bin
#' distances between cells [dlptools::pairwise_bin_differences()], this
#' function takes the nearest neighbour of each cell and fits a beta
#' distribution. Then using this distribution, it finds outlier cells based on
#' a selected percentile of the distribution (default 99th).
#'
#' @param cell_diffs dataframe of differences from
#' [dlptools::pairwise_bin_differences()]
#' @param nn_cells dataframe of [dlptools::find_nearest_neighbours()]
#' @param outlier_percentile double. Default 0.99. What percentile of the
#' distribution to consider an outlier cell.
#' @return NA or tibble of information on cells considered outliers. NA if no outliers found.
#' @export
find_outlier_cells <- function(
  pairwise_diffs,
  nn_cells,
  outlier_percentile = 0.99
) {
  diff_beta_dist <- fitdistrplus::fitdist(nn_cells$nn_diff, "beta")

  min_diff <- qbeta(
    p = outlier_percentile,
    shape1 = diff_beta_dist$estimate[["shape1"]],
    shape2 = diff_beta_dist$estimate[["shape2"]]
  )

  outlier_cells <- dplyr::filter(nn_cells, nn_diff >= min_diff) |>
    dplyr::rename(outlier_cell = index_cell, nn_cell = comp_cell)

  if (nrow(outlier_cells) == 0) {
    print("no outlier cells found!")
    return(NA)
  }

  return(outlier_cells)
}

#' visualize cell nearest neighbour distance values
#'
#' Generates a simple plot of the nearest neighbour distance for each cell,
#' which is highlighted with the point. Outlier cells are indicated by colour.
#' A violin plot is used to show distance to all other cells.
#'
#' @param pairwise_diffs dataframe of [dlptools::pairwise_bin_difference()]
#' @param nn_cells dataframe of [dlptools::find_nearest_neighbours()]
#' @param outlier_Cells dataframe of [dlptools::find_outlier_cells()]
#' @return ggplot object
#' @export
plot_nnd_outlier_cells <- function(
  pairwise_diffs,
  nn_cells,
  outlier_cells
) {
  nn_p_dat <- dplyr::mutate(
    nn_cells,
    `outlier cell` = index_cell %in% outlier_cells$outlier_cell
  )

  cell_dfs <- furrr::future_map_dfr(
    unique(nn_cells$index_cell),
    \(cell) {
      cell_df <- make_cell_focused_matrix(pairwise_diffs, cell)
    }
  )

  cell_ord <- nn_cells |>
    dplyr::arrange(nn_diff) |>
    dplyr::pull(index_cell) |>
    unique()

  nn_p_dat$index_cell <- factor(nn_p_dat$index_cell, levels = cell_ord)

  ggplot2::ggplot(
    nn_p_dat,
    ggplot2::aes(x = index_cell, y = nn_diff)
  ) +
    ggplot2::geom_violin(
      data = cell_dfs,
      bw = 0.02,
      ggplot2::aes(x = index_cell, y = prop_diff)
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = `outlier cell`),
      size = 3
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(color = "white"),
      axis.text.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::scale_color_manual(
      values = c(
        `FALSE` = GEN_PLOT_COLS[2],
        `TRUE` = GEN_PLOT_COLS[3]
      )
    ) +
    ggplot2::xlab("Index Cell") +
    ggplot2::ylab("Nearest Neighbour Distance") +
    ggplot2::guides(color = ggplot2::guide_legend(title = "Outlier Cell"))
}
