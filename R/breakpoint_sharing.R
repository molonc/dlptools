#' calculate jaccard similarity using breakpoint matrices
#'
#' Jaccard is the number of shared breakpoints (intersection), over the total
#' number of unique breakpoints (union) between two cells.
#'
#' @param bps_mtx sparse matrix. Breakpoint presence in cells
#' @param shared_bps sparse matrix. the numbers of shared breakpoints between
#' cells.
#' @return sparse matrix of jaccard similarity
calc_jaccard_similarity <- function(bps_mtx, shared_bps) {
  cell_bp_totals <- Matrix::rowSums(bps_mtx)
  intersection_bps <- as.matrix(shared_bps)
  card <- outer(cell_bp_totals, cell_bp_totals, FUN = "+")
  union_bps <- card - intersection_bps
  jaccard <- intersection_bps / union_bps
  jaccard[is.nan(jaccard)] <- 0 # for empty cells
  Matrix::Matrix(jaccard, sparse = TRUE)
}


#' summarize pairwise breakpoint sharing across cells
#'
#' @description
#' Will summarize the number of shared breakpoints (changes in copy number) for
#' all pairs of cells. Alternatively, will return jaccard similarity between
#' all pairs (intersection of shared breakpoints over the union of breakpoints).
#'
#' Returns a pairwise matrix that can be used for various downstream things,
#' like clustering.
#'
#' @param seg_df dataframe. Segmented copy number profiles per cell.
#' @param jaccard bool. True will calculate and return jaccard similarity.
#'
#' @return pairwise matrix of the number of shared breakpoints between each
#' pair of cells, or the jaccard similarity.
#'
#' @export
breakpoint_sharing <- function(
  seg_df,
  jaccard = FALSE
) {
  req_cols <- c("cell_id", "start", "chr")

  if (!all(req_cols %in% colnames(seg_df))) {
    stop(paste0("require columns of: ", paste(req_cols, sep = ", ")))
  }

  bps_df <- seg_df |>
    dplyr::filter(start != 1) |>
    dplyr::select(cell_id, chr, start) |>
    dplyr::mutate(
      chr_bp = stringr::str_c(chr, start, sep = ":"),
      bp = 1
    ) |>
    dplyr::select(cell_id, chr_bp)

  cell_factor <- factor(bps_df$cell_id)
  cell_ids <- levels(cell_factor)
  row_idx <- as.integer(cell_factor)
  breakpoints_fac <- factor(bps_df$chr_bp)
  breakpoints <- levels(breakpoints_fac)
  col_idx <- as.integer(breakpoints_fac)

  sparce_bp <- Matrix::sparseMatrix(
    i = row_idx,
    j = col_idx,
    x = 1,
    dims = c(
      length(cell_ids), length(breakpoints)
    ),
    dimnames = list(cell_ids, breakpoints)
  )

  shared_bps <- sparce_bp %*% Matrix::t(sparce_bp)

  if (jaccard) {
    return(calc_jaccard_similarity(sparce_bp, shared_bps))
  }

  return(shared_bps)
}
