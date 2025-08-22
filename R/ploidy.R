#' calculate sample ploidy with a weighted CN mean
#'
#' CN is first weighted by the length of its segment, a weighted mean of these
#' segments for each chromosome is made, then the mean across all chrosomomes
#' is made for each sample.
#'
#' In this way, the CN contribution to the ploidy is based on it's length, and
#' each chromosome contributes equally to the overall sample ploidy.
#'
#' @param cn_df dataframe input of segmented (or not) sample copy numbers
#' @param sample_col string column name identifying samples
#' @param start_col string column name identifying start of segments
#' @param end_col string column name identifying end of segments
#' @param cn_col string column name for copy number states
#' @param chrom_col string column name for chromosomes
#' @return tibble of per-sample ploidy estimates
#' @export
weighted_ploidy <- function(
    cn_df,
    sample_col = "cell_id",
    start_col = "start",
    end_col = "end",
    cn_col = "state",
    chrom_col = "chr") {
  sample_ploidy <- cn_df |>
    dplyr::ungroup() |>
    dplyr::mutate(
      span = cn_df[[end_col]] - cn_df[[start_col]],
      seg_cn_w = cn_df[[cn_col]] * span
    ) |>
    dplyr::group_by(.data[[sample_col]], .data[[chrom_col]]) |>
    dplyr::summarise(
      w_chr_cn = sum(seg_cn_w) / sum(span)
    ) |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::summarise(
      ploidy = mean(w_chr_cn)
    ) |>
    dplyr::ungroup()

  return(sample_ploidy)
}


#' find the mode CN per chromosome, then mode across the chromosomes
#'
#' Done per chromosome first so that large chromosomes don't dominate the
#' result. It is critical that the input be even bin level data, or this makes
#' no sense to do. Measures of the copy number states need to be done in evenly
#' sized bins.
#'
#' @param bin_df dataframe of bin level data
#' @param sample_col string column name identifying samples
#' @param cn_col string column name for copy number states
#' @param chrom_col string column name for chromosomes
#' @return tibble/dataframe of results by cell_id/sample
mode_ploidy <- function(
    bin_df,
    sample_col = "cell_id",
    cn_col = "state",
    chrom_col = "chr") {
  bin_df |>
    dplyr::group_by(.data[[sample_col]], .data[[chrom_col]]) |>
    dplyr::summarise(
      chr_mode = cust_mode(.data[[cn_col]])
    ) |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::summarise(
      mode_ploidy = cust_mode(chr_mode)
    ) |>
    dplyr::ungroup()
}
