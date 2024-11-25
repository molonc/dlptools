#' convert reads table to segments table
#'
#' Often the segments file given by the DLP pipeline is not what you want. What
#' you likely want is to do various filtering at the reads level, and then make
#' a segments file from that. This function will take a reads table and covert
#' it to a segments table.
#'
#' @param reads_df a table with standard reads data (e.g., created with
#' [import_dlp_files(file_type='reads')])
#' @return tibble/dataframe with read bins organized into segment blocks.
#' @export
#' @importFrom rlang .data
reads_to_segs <- function(reads_df) {
  new_segs <- reads_df %>%
    dplyr::select("cell_id", "chr", "start", "end", "state", "copy") %>%
    dplyr::group_by(.data$cell_id, .data$chr) %>%
    dplyr::mutate(
      rle_group = rle_states(.data$state)
    ) %>%
    dplyr::group_by(.data$cell_id, .data$chr, .data$rle_group) %>%
    dplyr::summarise(
      start = base::min(.data$start),
      end = base::max(.data$end),
      state = base::unique(.data$state), # will only be one state
      median_copy = stats::median(.data$copy)
    ) %>%
    dplyr::select(-c(rle_group))

  return(new_segs)
}

#' convert states to run length encoding
#'
#' takes a vector of numbers (e.g., states) and returns a numeric group
#' value indicating the a group for each run of the same value.
#' c(5,5,5,6,6,5,5,5,2) -> 1 1 1 2 2 3 3 3 4
#'
#' @param states really and vector of values.
#' @return vector of integers
rle_states <- function(states) {
  states_rle <- base::rle(states)
  rles <- base::rep(base::seq_along(states_rle$lengths), states_rle$lengths)
  return(rles)
}
