# idea is to use fill up/down
# would like to work in the concept of perhaps splitting large gaps, not sure
# how to handle that.
# but could specify splitting over the centromere gap. Basically, split by
# cell + chrom for filling, but can then optionally further split by being
# within the centromere buffer and then can do a downup fill to split the gap
# ...or can I? How would I tell it where to split? Would need to calculate the
# mid point of the centromere.

#' loading UCSC chromosome length files
#' @param version default "hg19", can also load "hg38"
#' @return tibgble of chromosome, total length, etc.
#' @export
load_chrom_info_file <- function(version = c("hg19", "hg38")) {
  version <- match.arg(version)

  chrom_files <- list(
    hg19 = "hg19_chromInfo.txt.gz",
    hg38 = "hg19_chromInfo.txt.gz"
  )

  chrom_info <- get_package_file_path(chrom_files[version]) |>
    vroom::vroom(
      col_names = c("chr", "total_length", "misc"),
      show_col_types = FALSE
    ) |>
    dplyr::filter(
      # remove the unnecessary chromosomes
      stringr::str_detect(chr, "_|M", negate = TRUE)
    ) |>
    dplyr::mutate(
      chr = stringr::str_replace(chr, "chr", "")
    )

  return(chrom_info)
}

#' For a given length, create bins of a given size
#'
#' e.g., length = 10, bin = 5, will get bins: 1-5, 6-10
#'
#' @param total_len int to create bins for
#' @param bin_size int of size of bin
#' @return tibble of bins with columns of start and end
#' @export
expand_length_to_bins <- function(total_len, bin_size = 5e5) {
  starts <- seq(1, total_len, by = bin_size)
  ends <- seq(bin_size, total_len + bin_size, by = bin_size)
  if ((total_len + bin_size) %% bin_size == 0) {
    # TODO: this seems unnecessary. Better seq to fix?
    ends <- ends[seq_along(ends) - 1]
  }
  bins <- tibble::tibble(start = starts, end = ends)
  return(bins)
}

#' create bins of given size for a genome
#'
#' @param version default "hg19", or can select "hg38"
#' @param bins_size int of how big to make the bins
#' @return tibble of bins for each chromosome. Columns: chr start end
#' @export
create_expected_bins <- function(version = c("hg19", "hg38"), bin_size = 5e5) {
  version <- match.arg(version)

  chroms <- load_chrom_info_file(version)

  chrom_bins <- purrr::pmap(
    chroms,
    \(chr, total_length, ...) {
      bins <- expand_length_to_bins(total_length)
      bins |>
        dplyr::mutate(chr = chr) |>
        dplyr::relocate(chr)
    }
  ) |>
    purrr::list_rbind()

  return(chrom_bins)
}


#' Add back chr, start, end of bins missing from cells
#'
#' Often we filter read bins based on various criteria, or tools we use might
#' drop bins of read/state data. This function will put missing bins back into
#' a dataframe for cells.
#'
#' Basic operation will result in NAs for every column except cell_id, chr,
#' start, end for the newly added bins. But you can also carry over cell
#' metadata by specifying which columns to keep for the inferred bins. This
#' should be cell level data.
#'
#' @param state_df the dataframe/tibble to insert the missing bins into
#' @param version which genome version you are using. Default "hg19", or select
#' "hg38"
#' @param bin_size size of bins that should be there.
#' @param cell_metadata_cols vector of columns to carry through to inferred bins
#' @param input_chroms_only bool. Only insert bins for chromosomes in the input
#' state dataframe. I.e., if you've filtered X or Y, don't add bins for X or Y.
#' Setting to FALSE adds bins for entire genome. Default: TRUE
#' @return input dataframe with extra bins for each cells that were missing.
#' @export
add_missing_bins_for_cells <- function(
    state_df, version = c("hg19", "hg38"), bin_size = 500000,
    cell_metadata_cols = c(),
    input_chroms_only = TRUE) {
  version <- match.arg(version)
  exp_bins <- create_expected_bins(version, bin_size)

  if (input_chroms_only) {
    exp_bins <- dplyr::filter(exp_bins, chr %in% unique(state_df$chr))
  }

  cell_bins <- dplyr::cross_join(
    # could optionally carry in other columns too
    dplyr::distinct(state_df, dplyr::across(c("cell_id", cell_metadata_cols))),
    exp_bins
  )

  state_df <- dplyr::full_join(
    state_df, cell_bins,
    by = c(cell_metadata_cols, "cell_id", "chr", "start", "end")
  ) |>
    dplyr::arrange(cell_id, chr, start)

  return(state_df)
}

#' infer missing data using the upstream neighbour bin
#'
#' This function fills in NAs of columns with the information from an upstream
#' bin. The exception are bins missing from the start of chromosomes, which are
#' filled by their downstream neighbour if there is nothing upstream.
#'
#' Primarily used to fill in "state" data for filtered bins. But can also pass
#' other states to fill in.
#'
#' @param state_df dataframe of state data with NAs for some bins that need
#' to be filled in. E.g., after dlptools::running add_missing_bins_for_cells()
#' @export
#' @importFrom rlang .data
fill_state_from_neighbours <- function(state_df, cols_to_fill = c("state")) {
  filled_df <- state_df |>
    dplyr::group_by(.data$cell_id, .data$chr) |>
    tidyr::fill(cols_to_fill, .direction = "downup") |>
    dplyr::ungroup()

  return(filled_df)
}
