# code to perform a light joint segmentation on read bin data
# Mostly useful for aligning hmmcopy segments. Primary effect is dropping
# 500kb variability in segment cutoffs.

#' convert a dataframe to a wide matrix of bin states
#' @return list of Y, the matrix, and cell_ord, the cell IDs in the matrix
convert_to_mtx <- function(chr_df) {
  wide_res <- dlptools::convert_long_reads_to_wide(chr_df)
  mtx <- wide_res |>
    dplyr::select(-cell_id) |>
    t()
  return(list(
    Y = mtx,
    cell_ord = wide_res$cell_id
  ))
}

#' alter existing segments to fit new breakpoints
#' @param jfit a [jointseg::jointSeg()] output object
#' @param input_mtx matrix. cell by bin state matrix
#' @param input_cellids vector. Cell ids of cells in the matrix.
#' @return tibble. Cell segments that align with new breakpoints
snap_segs_to_new_breakpoints <- function(
  jfit,
  input_mtx,
  input_cellids
) {
  starts <- c(1, jfit$bestBkp + 1)
  ends <- c(jfit$bestBkp, nrow(input_mtx))

  snapped_profile_mat <- matrix(
    NA_integer_,
    nrow = nrow(input_mtx),
    ncol = ncol(input_mtx),
    dimnames = dimnames(input_mtx)
  )

  snapped_segment_matrix <- matrix(
    NA_integer_,
    nrow = length(starts),
    ncol = ncol(input_mtx),
    dimnames = list(paste0("Seg_", seq_along(starts)), colnames(input_mtx))
  )

  for (i in seq_along(starts)) {
    seg_range <- starts[i]:ends[i]

    if (length(seg_range) == 1) {
      # Single-bin segment: take values directly
      majority_cn <- input_mtx[seg_range, ]
    } else {
      # Multi-bin segment: calculate majority integer call per cell
      majority_cn <- apply(
        input_mtx[seg_range, , drop = FALSE], 2, dlptools::cust_mode
      )
    }

    # Assign discrete segment state
    snapped_segment_matrix[i, ] <- majority_cn

    # Expand back into the full profile matrix
    snapped_profile_mat[seg_range, ] <- matrix(
      rep(majority_cn, length(seg_range)),
      nrow = length(seg_range),
      byrow = TRUE
    )
  }

  colnames(snapped_profile_mat) <- input_cellids
  new_segs <- tibble::as_tibble(snapped_profile_mat) |>
    dplyr::mutate(
      bin = rownames(snapped_profile_mat)
    ) |>
    tidyr::separate_wider_delim(
      bin,
      delim = "_",
      names = c("chr", "start", "end")
    ) |>
    dplyr::relocate(chr, start, end) |>
    dplyr::mutate(
      dplyr::across(c(start, end), as.numeric)
    ) |>
    tidyr::pivot_longer(
      cols = -c(chr, start, end),
      names_to = "cell_id",
      values_to = "state"
    ) |>
    dlptools::reads_to_segs()

  return(new_segs)
}


#' Perform joint segmentation on a single chromosome
#'
#' @details see [dlptools::joint_seg_reads] for explanation. But effectively,
#' joint segmentation needs to happen per chromosome. So this function does the
#' segmentation, and the other is just a wrapper for all chromosomes.
#'
#' @param chr_df dataframe. State calls of the DLP bins data.
#' @param chrom_col string. Name of the chromosome column
#' @param max_bps int. Can specify the maximum number of breakpoints to find
#' @param bin_fraction float. What fraction of the bins to target as the maximum number of segments.
#'
#' @return dataframe. Cell segments that match the jointly found segments.
#' @export
joint_seg_chromosome <- function(
  chr_df,
  chrom_col = "chr",
  max_bps = NULL,
  bin_fraction = 0.5
) {
  working_chr <- unique(chr_df[[chrom_col]])

  stopifnot(
    "more than one chrom passed for joint seg" = length(working_chr) == 1
  )

  mtx_info <- convert_to_mtx(chr_df)
  mtx <- mtx_info$Y
  mtx_cellids <- mtx_info$cell_ord

  if (is.null(max_bps)) {
    max_bps <- round(bin_fraction * nrow(mtx))
  }

  jfit <- jointseg::jointSeg(
    mtx,
    method = "RBS",
    K = max_bps
  )

  found_bk_prop <- length(jfit$bestBkp) / max_bps
  if (found_bk_prop >= 0.95) {
    warning(paste0(
      "num breakpoints found for chrom ", working_chr,
      " is within 5% of max. Might want to try higher number"
    ))
  } else if (found_bk_prop == 1) {
    stop(paste0(
      "Number of breakpoits for ", working_chr,
      "hit the max possible. Try again with a larger search space"
    ))
  }

  new_segs <- snap_segs_to_new_breakpoints(jfit, mtx, mtx_cellids)

  return(new_segs)
}


#' Perform joint segmentation of read state calls
#'
#' @description
#' This function takes a cell bin state calls and uses recursive binary
#' segmentation to find as many breakpoints as it can. By default, it looks for
#' a maximum number of breakpoints equal to half the number of input bins.
#' E.g., 100 bins of data, will search for at most 50 breakpoints.
#'
#' Joint segmentation happens across all cells for a single chromosome at a
#' time. This function is a wrapper for all chromosomes, calling
#' [dlptools::joint_seg_chromosome] for each to do the joint segmentation.
#'
#' This function was intended for use with the already segmented profiles of
#' hmmcopy, by leveraging the state calls. And the goal is to modify them as
#' little as possible to correct for slight variances among cells in where the
#' segment breakpoints are inferred.
#'
#' @param reads_df dataframe. Binned read data.
#' @param chrom_col string. Name of the chromosome column
#' @param max_bps NULL/int. Can specify the maximum number of breakpoints to
#' find
#' @param bin_fraction float. What fraction of the bins to target as the
#' maximum number of segments.
#'
#' @return dataframe. Cell segments that match the jointly found segments.
#' @export
joint_seg_reads <- function(
  reads_df,
  chrom_col = "chr",
  max_bps = NULL,
  bin_fraction = 0.5
) {
  # safely map incase a chromosome fails, then can redo just that chromosome
  safely_joint_seg_chr <- purrr::safely(joint_seg_chromosome)

  joint_seg_res <- reads_df |>
    dplyr::group_split(.data[[chrom_col]]) |>
    furrr::future_map(
      safely_joint_seg_chr,
      chrom_col = chrom_col,
      max_bps = max_bps,
      bin_fraction = bin_fraction
    )

  joint_segs <- purrr::map_dfr(joint_seg_res, \(r) r$result)

  chrom_errs <- purrr::map(joint_seg_res, "error") |>
    purrr::compact() |>
    paste()

  if (length(chrom_errs) > 0) {
    print(chrom_errs)
    print(paste0(
      "Some chromosomes failed. Re-rurn chromosomes directly with",
      " dlptools::joint_seg_chromosome(chrom_reads_df)"
    ))
  }

  return(joint_segs)
}
