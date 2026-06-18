#' cut segments to equal widths across cells/samples.
#'
#' Segments are continuous stretches of a chromosome with the same copy number.
#' But, generally, these will be different across cells/samples. This function
#' will create segments in the smallest consistent intervals possible across
#' the samples. See example, but it's basically equalizing the segment widths
#' among all of the cells, so copy number is measured in the same "bin" across
#' cells.
#'
#' The main value of this is that some tools (medicc, or plotting) want the same
#' intervals across all samples (i.e., bins of equal width). But going all the
#' way back to 500 kb bins (i.e., typical DLP) can be excessive when there are
#' longer segments within cells. Cutting segments into their minimal spans will
#' generally be a smaller set of data than a full 500kb bin breakdown (unless
#' CN data is highly noisy among samples).
#'
#' @examples
#' sgs <- tibble::tibble(
#'   cell_id = c("A", "A", "B", "B"),
#'   chrom = rep("chr1", 4),
#'   start = c(1, 11, 1, 8),
#'   end = c(10, 25, 7, 25),
#'   state = c(2, 4, 3, 8)
#' )
#'
#' sgs
#'
#' align_segments(sgs)
#'
#' @param segs_df dataframe. Copy number segments of cells/samples. Need
#' chromosome, start, end.
#' @param chrom_col string. name of column containing chromosome information
#' @param cell_col string. name of column with cell/sample IDs
#' @return tibble of broken down segments with same input columns + seg_width
#' of the new segments.
#' @export
align_segments <- function(
  segs_df,
  chrom_col = "chr",
  cell_col = "cell_id"
) {
  ply1 <- plyranges::as_granges(
    segs_df,
    seqnames = !!rlang::sym(chrom_col)
  )

  dis1 <- GenomicRanges::disjoin(ply1)

  align1 <- plyranges::join_overlap_intersect(dis1, ply1)

  align1 |>
    as.data.frame() |>
    tibble::tibble() |>
    dplyr::rename(!!chrom_col := seqnames) |>
    dplyr::select(dplyr::all_of(names(segs_df))) |>
    dplyr::mutate(
      seg_width = end - start
    ) |>
    dplyr::arrange(.data[[cell_col]])
}
