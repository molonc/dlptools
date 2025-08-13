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
  if (any(is.na(reads_df$state))) {
    stop("Error: 'state' column in reads_df contains NA values.")
  }
  new_segs <- reads_df |>
    dplyr::select("cell_id", "chr", "start", "end", "state") |>
    dplyr::group_by(.data$cell_id, .data$chr) %>%
    dplyr::mutate(
      rle_group = rle_states(.data$state)
    ) |>
    dplyr::group_by(.data$cell_id, .data$chr, .data$rle_group) |>
    dplyr::summarise(
      start = base::min(.data$start),
      end = base::max(.data$end),
      state = base::unique(.data$state), # will only be one state
    ) |>
    dplyr::select(-c(rle_group)) |>
    dplyr::mutate(
      seg_width = end - start
    ) |>
    dplyr::ungroup()

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
#' @export
rle_states <- function(states) {
  states_rle <- base::rle(states)
  rles <- base::rep(base::seq_along(states_rle$lengths), states_rle$lengths)
  return(rles)
}


#' split segments into bins of a requested size.
#'
#' Takes a dataframe with segment start and end columns and returns a dataframe
#' with those segments split into individual bins. All segment information
#' applied to the generated bins.
#'
#' Warnings issues if requested bin size is bigger than some segments or if
#' segments can't be split evenly into the bins.
#'
#' Bin end will not exceed a segment end. Depending on what is input and
#' requested, this can lead to one smaller bin at the end of the segment.
#'
#' @param segs_df dataframe (or similar) of segment dat
#' @param bin_size width of bins to split into. Default is standard 500kb
#' @param seg_start_col name of the column that indicates the start of a segment
#' @param seg_end_col name of the column that indicates the end of a segment
#' @importFrom rlang .data
#' @return input frame split into bins
#' @export
segs_to_reads <- function(
    segs_df, bin_size = 5e5, seg_start_col = "start", seg_end_col = "end") {
  binned_segs <- segs_df |>
    dplyr::mutate(
      seg_width = .data[[seg_end_col]] - .data[[seg_start_col]] + 1,
      bins = purrr::map(seg_width, \(width) {
        expand_length_to_bins(width, bin_size = bin_size) |>
          dplyr::rename(bin_start = start, bin_end = end)
      })
    ) |>
    tidyr::unnest(bins) |>
    dplyr::mutate(
      bin_start = bin_start + .data[[seg_start_col]] - 1,
      bin_end = bin_end + .data[[seg_start_col]] - 1,
      short_seg = seg_width < bin_size,
      bin_end = dplyr::if_else(
        # reset any bins that exceed end of segment
        bin_end > .data[[seg_end_col]], .data[[seg_end_col]], bin_end
      ),
      uneven_bin = bin_end - bin_start + 1 < bin_size
    ) |>
    dplyr::rename(
      seg_start = .data[[seg_start_col]],
      seg_end = .data[[seg_end_col]],
      start = bin_start,
      end = bin_end
    )

  if (TRUE %in% unique(binned_segs$uneven_bin)) {
    warning(paste0(
      "Some segments unable to be evenly split into requested bin sizes.",
      " Use df$uneven_bin column to find which ones."
    ))
  }

  if (TRUE %in% unique(binned_segs$short_seg)) {
    warning(paste0(
      "Some segments shorter than the requested bin sizes.",
      " Use df$short_seg column to find which ones."
    ))
  }

  return(binned_segs)
}
