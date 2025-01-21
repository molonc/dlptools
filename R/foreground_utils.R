#' mark the types of change between a parent node state and the tip.
#'
#' This compares the tip and parent state columns to evaluate if a foreground
#' change has occurred, what type of change, and how much of a change.
#'
#' e.g., if parent state is 5, tip is 4, it will get marked as a 'foreground'
#' event, fg-gain, with a change of 1.
#'
#' 4 coulmns are added:
#'  foreground: boolean, if parent and tip states are different
#'  fg_type: character, fg-gain, fg-loss, or NA
#'  fg_change: numeric, tip - parent state. Positive indicates gain in tip.
#'  background: boolean, inverse of foreground
#'
#' @param bin_states_df bin based dataframe with parent and tip state columns.
#' @param state_col name of tip state column
#' @param parent_state_col name of parent state column.
#' @return input dataframe with new columns
#' @export
characterize_foreground_total_state_changes <- function(
    bin_states_df, state_col = "state", parent_state_col = "parent_state") {
  # confirm that we are working with numbers or ordered factors
  check_if_state_columns_are_valid(bin_states_df, state_col, parent_state_col)

  bin_states_df <- bin_states_df |>
    dplyr::mutate(
      foreground = .data[[state_col]] != .data[[parent_state_col]],
      fg_type = dplyr::case_when(
        .data[[state_col]] > .data[[parent_state_col]] ~ "fg-gain",
        .data[[state_col]] < .data[[parent_state_col]] ~ "fg-loss",
        .default = NA
      ),
      fg_change = .data[[state_col]] - .data[[parent_state_col]],
      background = !foreground
    )

  return(bin_states_df)
}

#' mark types of allele specific foreground state changes
characterize_foreground_allele_state_changes <- function(bin_states_df) {
  # TODO: implement this.
  warning(
    "Haven't implemented characterizing allele state foreground changes yet."
  )
  return(bin_states_df)
}

#' confirms a vector is not numeric and is unorderd
#'
#' returns TRUE if the vector is not numeric and is missing a defined order.
#' Used as a check of suitability for comparing events between tip and parent
#' nodes.
is_not_numeric_lacks_order <- function(vec) {
  return(!is.numeric(vec[1]) && !is.ordered(vec))
}

#' internal check for using the characterization function
#'
#' To characterize forground changes, need wither numeric columns or factors
#' with a defined order to them.
check_if_state_columns_are_valid <- function(
    states_df, state_col = "state", parent_state_col = "parent_state") {
  tip_check <- is_not_numeric_lacks_order(states_df[[state_col]])
  parent_check <- is_not_numeric_lacks_order(states_df[[parent_state_col]])

  if (tip_check || parent_check) {
    stop(paste0(
      "Passed columns to summarize are not numeric, but some sort of strings.",
      " These need to be ordered factors that explain the relationships",
      " between the discrete states. E.g., ",
      " factor(",
      "c('a', 'a', 'b', 'c', 'c'),",
      " levels=c('c', 'a', 'b'),",
      " ordered=TRUE",
      "),",
      " which would define c < a < b."
    ))
  }
}
