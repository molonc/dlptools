#' make a plot of just the foreground events
#'
#' Colors highlight the change in state between the tip and parent node.
#'
#' @param states_df bin based dataframe of state changes. Typically the output
#' of characterize_foreground_total_state_changes() or from
#' medicc_profiles_to_foreground(). But really, any bin dataframe with a column
#' indicating change between tip and parent will work.
#' @param phylogeny phylo object of the tree being used for analysis.
#' @param file_name what to call the resulting plot png.
#' @param fg_change the column to plot with the changes between parent and tips.
#' @export
plot_heatmap_of_tip_changes <- function(
    states_df, phylogeny, file_name, changes_col = "fg_change") {
  min_state <- min(states_df[[changes_col]], na.rm = TRUE)
  max_state <- max(states_df[[changes_col]], na.rm = TRUE)

  # temp until I get more colors going...should probably have up to 11 change
  # colors in both directions
  if (min_state < -8 || max_state > 8) {
    warning("Found state changes beyond 8 in one direction. Don't have colors for that. Will cap values at 8+ and -8+")
    states_df[[changes_col]][states_df[[changes_col]] >= 8] <- "8+"
    states_df[[changes_col]][states_df[[changes_col]] <= -8] <- "-8+"
  }

  loss_colors <- c(
    "#77bcf0",
    "#3182BD",
    "#0C46A0FF",
    "#470ca0",
    "#30b939",
    "#23842a",
    "#a8a8a8",
    "#000000",
    "#000000"
  )
  names(loss_colors) <- c(-1:-8, "-8+")

  gain_colors <- c(
    CNV_COLOURS[4:length(CNV_COLOURS)][1:8],
    CNV_COLOURS[4:length(CNV_COLOURS)][8]
  )
  names(gain_colors) <- c(1:8, "8+")

  fg_change_palette <- c("0" = "#ffffff", gain_colors, loss_colors)

  plot_state_hm(
    states_df,
    state_col = changes_col,
    phylogeny = phylogeny,
    file_name = file_name,
    hm_discrete_colors = fg_change_palette
  )
}

#' generate a plot of just the foreground states, making the background white.
#'
#' These are the actual state values of the bins in foreground events.
#'
#' @param states_df bin based dataframe that has foreground events marked.
#' @param file_name name for the resulting png.
#' @export
plot_fg_state_highlight <- function(states_df, file_name) {
  states_df |>
    dplyr::mutate(
      fg_highlight = dplyr::if_else(
        foreground, as.character(state), "background"
      )
    ) |>
    plot_state_hm(
      state_col = "fg_highlight",
      phylogeny = os052_p9_tree,
      file_name = file_name,
      hm_discrete_colors = c(dlptools::CNV_COLOURS, `background` = "#f4f8f6")
    )
}

#' generate a plot of just the background states, making the foreground white.
#'
#' These are the actual state values of the bins in background bins.
#'
#' @param states_df bin based dataframe that has foreground events marked.
#' @param file_name name for the resulting png.
#' @export
plot_bg_state_highlight <- function(states_df, file_name) {
  states_df |>
    dplyr::mutate(
      bg_highlight = dplyr::if_else(
        background, as.character(state), "foreground"
      )
    ) |>
    plot_state_hm(
      state_col = "bg_highlight",
      phylogeny = os052_p9_tree,
      file_name = file_name,
      hm_discrete_colors = c(dlptools::CNV_COLOURS, `foreground` = "#f4f8f6")
    )
}
