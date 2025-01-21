characterize_foreground_total_state_changes <- function(bin_states_df, state_col = "state", parent_state_col = "parent_state") {
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

characterize_foreground_allele_state_changes <- function(bin_states_df) {
  warning(
    "Haven't implemented characterizing allele state foreground changes yet."
  )
  return(bin_states_df)
}


is_not_numeric_lacks_order <- function(vec) {
  return(!is.numeric(vec[1]) && !is.ordered(vec))
}

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
