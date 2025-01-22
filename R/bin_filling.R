# idea is to use fill up/down
# would like to work in the concept of perhaps splitting large gaps, not sure
# how to handle that.
# but could specify splitting over the centromere gap. Basically, split by
# cell + chrom for filling, but can then optionally further split by being
# within the centromere buffer and then can do a downup fill to split the gap
# ...or can I? How would I tell it where to split? Would need to calculate the
# mid point of the centromere.

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

#' @export
expand_length_to_bins <- function(total_len, bin_size = 5e5) {
  starts <- seq(1, total_len, by = bin_size)
  ends <- seq(bin_size, total_len + bin_size, by = bin_size)
  bins <- tibble::tibble(start = starts, end = ends)
  return(bins)
}

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


#' @export
add_missing_bins_for_cells <- function(
    state_df, version = c("hg19", "hg38"), bin_size = 500000,
    cell_metadata_cols = c()) {
  version <- match.arg(version)
  exp_bins <- create_expected_bins(version, bin_size)

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
