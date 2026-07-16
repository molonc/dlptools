#' convert reads table to segments table
#'
#' Often the segments file given by the DLP pipeline is not what you want. What
#' you likely want is to do various filtering at the reads level, and then make
#' a segments file from that. This function will take a reads table and covert
#' it to a segments table.
#'
#' @param reads_df a table with standard reads data
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


#' original function to split segments back to read bins
#'
#' kept for posterity. This splits into exact 500kb bins, which may extend
#' beyond end of chromosome, as original DLP did. New function is wildly faster
#' and more simple.
#' @export
old_segs_to_reads <- function(
  segs_df, bin_size = 5e5, seg_start_col = "start", seg_end_col = "end"
) {
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

#' Chop Genomic Segments into Fixed-Size Genomic Bins
#'
#' @description
#' Splits variable-length copy number segments into fixed-size genomic bins
#' (windows) of a specified width (e.g., 500 Kb).
#'
#' @param segs_df A data.frame or tibble containing genomic segments with
#'   start, end, and chromosome name columns.
#' @param bin_size Numeric. The width of the genomic windows in base pairs.
#'   Defaults to \code{5e5} (500 Kb).
#' @param genome_version Character. The assembly build to use for creating
#'   genomic windows. Must be one of \code{"hg19"} or \code{"hg38"}.
#' @param return_type Character. The format of the returned object. Must be
#'   one of \code{"granges"} or \code{"tibble"}. Defaults to \code{"tibble"}.
#' @param seg_start_col Character. The name of the column in \code{segs_df}
#'   representing segment start coordinates. Defaults to \code{"start"}.
#' @param seg_end_col Character. The name of the column in \code{segs_df}
#'   representing segment end coordinates. Defaults to \code{"end"}.
#' @param sample_id_col Character. The name of the column in \code{segs_df}
#'   identifying the sample or cell. Defaults to \code{"cell_id"}.
#' @param state_col Character. The name of the column in \code{segs_df}
#'   representing the copy number state. Defaults to \code{"state"}.
#' @param other_meta_cols Character vector. Optional additional metadata column
#'   names in \code{segs_df} to carry over to the output. Defaults to an empty
#' vector.
#' @param chrom_col Character. The name of the column in \code{segs_df}
#'   representing chromosomes. Defaults to \code{"chr"}.
#'
#' @return Depending on \code{return_type}:
#'   \itemize{
#'     \item \code{"granges"}: A \code{\link[GenomicRanges]{GRanges}} object
#'       containing the chopped bins with associated metadata columns.
#'     \item \code{"tibble"}: A \code{\link[tibble]{tibble}} version of the
#'       GRanges object, with the chromosome column renamed back to
#'        \code{chrom_col}.
#'   }
#'   In both formats, the output retains the original segment start and end
#'   positions in the \code{seg_start} and \code{seg_end} metadata columns.
#'
#' @export
segs_to_reads <- function(
  segs_df,
  bin_size = 5e5,
  genome_version = c("hg19", "hg38"),
  return_type = c("tibble", "granges"),
  seg_start_col = "start",
  seg_end_col = "end",
  sample_id_col = "cell_id",
  state_col = "state",
  other_meta_cols = c(),
  chrom_col = "chr"
) {
  return_type <- match.arg(return_type)

  bins_tile <- create_chrom_window_intervals(
    window_size = bin_size,
    genome_version = genome_version,
    return_type = "granges"
  )

  segs_grs <- cap_dlp_to_chrom_lengths(
    segs_df,
    chrom_col = chrom_col,
    genome_version = genome_version
  ) |>
    GenomicRanges::GRanges()

  hits <- GenomicRanges::findOverlaps(bins_tile, segs_grs)

  chopped_segments <- GenomicRanges::pintersect(
    bins_tile[S4Vectors::queryHits(hits)],
    segs_grs[S4Vectors::subjectHits(hits)]
  )
  S4Vectors::mcols(chopped_segments)$hit <- NULL

  # update with meta data
  purrr::walk(
    c(state_col, sample_id_col, other_meta_cols),
    \(targ_col) {
      (
        S4Vectors::mcols(chopped_segments)[[targ_col]] <<-
          S4Vectors::mcols(segs_grs)[[targ_col]][S4Vectors::subjectHits(hits)]
      )
    }
  )

  # re-insert segment info
  S4Vectors::mcols(chopped_segments)$seg_start <-
    S4Vectors::start(segs_grs)[S4Vectors::subjectHits(hits)]
  S4Vectors::mcols(chopped_segments)$seg_end <-
    S4Vectors::end(segs_grs)[S4Vectors::subjectHits(hits)]

  if (return_type == "granges") {
    chopped_segments
  } else {
    tibble::as_tibble(chopped_segments) |>
      dplyr::rename(!!chrom_col := seqnames) |>
      dplyr::relocate(
        dplyr::all_of(c(sample_id_col, chrom_col)),
        start, end
      )
  }
}


#' cap DLP data at the lengths of chromosomes
#'
#' DLP can output bins/segments longer than the actual chromosomes as it chops
#' to 500kb bins. This function sets the \code{end} column to be no longer than
#' the chromosome it occurs on.
#'
#' @param rs_df tibble/dataframe. Reads or segments dataframe
#' @param chrom_col string. Name of chromosome column
#' @param genome_version string. Name of genome version to use. Default: hg19
#'
#' @return input dataframe with final bins/segments not exceeding total
#' chromosome lengths
#' @export
cap_dlp_to_chrom_lengths <- function(
  rs_df,
  chrom_col = "chr",
  genome_version = c("hg19", "hg38")
) {
  chr_info <- suppressWarnings(load_chrom_info_file(version = genome_version))

  if (class(rs_df[[chrom_col]]) != "character") {
    rs_df[[chrom_col]] <- as.character(rs_df[[chrom_col]])
  }

  # DLP can output bins longer than the chromosome to maintain it's 500Kb
  # bin size
  rs_df_cap <- dplyr::left_join(
    rs_df, chr_info,
    by = setNames(chrom_col, "chr")
  ) |>
    dplyr::mutate(
      end = dplyr::if_else(end > total_length, total_length, end)
    )

  rs_df_cap
}
