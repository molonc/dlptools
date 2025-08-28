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
#' @export
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

#' Mark if CN states are gains or losses relative to cell ploidy
#'
#' CNs > 2 are not necessarily amplifications, and CNs < 2 are not necessarily
#' the only states of losses. Ploidy of samples may not be diploid. This
#' function will give you an idea if the CN you see is a gain or loss relative
#' to the ploidy of a sample. For example, if the sample has a ploidy of 4,
#' then a CN of 3 is a loss.
#'
#' Ploidy is inferred by mode CN state using [dlptools::mode_ploidy()] and
#' states > ploidy are marked as gains, states < ploidy losses, and states
#' matching ploidy as matched.
#'
#' @param in_df dataframe of CN states.
#' @param df_type string. "reads" (default) or "segs" for CN segments, which
#' will internally converted to bin based for mode calculation.
#' reads in order to infer mode ploidy.
#' @param sample_col string. Name of the column with cell_id/other sample name
#' @return input dataframe, with new columns of information.
#' @export
mark_cn_relative_to_ploidy <- function(
    in_df,
    df_type = c("reads", "segs"),
    sample_col = "cell_id",
    ...) {
  df_type <- match.arg(df_type)

  event_labels <- c(
    gain = "ploidy-gain",
    match = "ploidy-match",
    loss = "ploidy-loss"
  )

  if (df_type == "segs") {
    reads_df <- segs_to_reads(in_df, ...)
  } else if (df_type == "reads") {
    reads_df <- in_df
  }

  sample_mode_ploidy <- mode_ploidy(reads_df, sample_col = sample_col, ...)

  in_df <- dplyr::left_join(
    in_df,
    sample_mode_ploidy,
    by = sample_col
  )

  in_df <- in_df |>
    dplyr::mutate(
      cn_v_ploidy = dplyr::case_when(
        state < mode_ploidy ~ event_labels["loss"],
        state > mode_ploidy ~ event_labels["gain"],
        state == mode_ploidy ~ event_labels["match"]
      ),
      cn_v_ploidy = factor(cn_v_ploidy, levels = unname(event_labels))
    )

  return(in_df)
}
